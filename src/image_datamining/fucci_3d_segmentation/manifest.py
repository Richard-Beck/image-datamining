from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


CHANNELS = ("00", "01", "02")
Z_RE = re.compile(r"_z(?P<z>\d+)_ch(?P<channel>\d+)\.tif$", re.IGNORECASE)


@dataclass(frozen=True)
class FieldRecord:
    sample_id: str
    cell_line: str
    source_dir: Path
    prefix: str
    n_z: int
    z_min: int
    z_max: int
    ch00_files: str
    ch01_files: str
    ch02_files: str


def cell_line_from_path(path: Path) -> str:
    parts = {part.lower(): part for part in path.parts}
    if "nci-n87" in parts:
        return "NCI-N87"
    if "mdamb" in parts:
        return "MDAMB"
    return ""


def sample_id_from_dir(path: Path) -> str:
    cell_line = cell_line_from_path(path)
    raw = path.name
    raw = raw.replace(".nucleus", "_nucleus")
    raw = re.sub(r"[^A-Za-z0-9]+", "_", raw).strip("_")
    return f"{cell_line}_{raw}" if cell_line else raw


def field_dirs_under(root: Path, recursive: bool = False) -> Iterable[Path]:
    if recursive:
        yield from sorted(p for p in root.rglob("*") if p.is_dir())
    else:
        yield from sorted(p for p in root.iterdir() if p.is_dir())


def discover_z_stack_fields(roots: Iterable[Path], recursive: bool = False) -> list[FieldRecord]:
    records: list[FieldRecord] = []
    seen: set[Path] = set()
    for root in roots:
        if not root.exists():
            continue
        for field_dir in field_dirs_under(root, recursive=recursive):
            resolved = field_dir.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            by_channel: dict[str, list[tuple[int, Path]]] = {channel: [] for channel in CHANNELS}
            prefix = ""
            for tif in sorted(field_dir.glob("*.tif")):
                match = Z_RE.search(tif.name)
                if not match:
                    continue
                channel = match.group("channel")
                if channel not in by_channel:
                    continue
                z = int(match.group("z"))
                by_channel[channel].append((z, tif))
                if not prefix:
                    prefix = tif.name[: match.start()]

            if not all(by_channel.values()):
                continue

            z_sets = [{z for z, _path in by_channel[channel]} for channel in CHANNELS]
            common_z = sorted(set.intersection(*z_sets))
            if not common_z:
                continue

            def paths_for(channel: str) -> str:
                lookup = dict(by_channel[channel])
                return ";".join(str(lookup[z]) for z in common_z)

            records.append(
                FieldRecord(
                    sample_id=sample_id_from_dir(field_dir),
                    cell_line=cell_line_from_path(field_dir),
                    source_dir=field_dir,
                    prefix=prefix,
                    n_z=len(common_z),
                    z_min=common_z[0],
                    z_max=common_z[-1],
                    ch00_files=paths_for("00"),
                    ch01_files=paths_for("01"),
                    ch02_files=paths_for("02"),
                )
            )
    return records


def write_manifest(records: Iterable[FieldRecord], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(FieldRecord.__dataclass_fields__.keys())
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record in records:
            row = {
                "sample_id": record.sample_id,
                "cell_line": record.cell_line,
                "source_dir": str(record.source_dir),
                "prefix": record.prefix,
                "n_z": record.n_z,
                "z_min": record.z_min,
                "z_max": record.z_max,
                "ch00_files": record.ch00_files,
                "ch01_files": record.ch01_files,
                "ch02_files": record.ch02_files,
            }
            writer.writerow(row)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

#!/usr/bin/env python3
"""Bounded directory-tree walk for large shared microscopy folders.

Each directory is scanned with GNU find under timeout. Directories that time out
are marked as black holes, sampled cheaply, and not explored further.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any


def extension_for(name: str) -> str:
    suffix = Path(name).suffix.lower()
    if suffix:
        return suffix
    return "[no extension]"


def depth_for(path: Path, root: Path) -> int:
    rel = path.relative_to(root)
    if str(rel) == ".":
        return 0
    return len(rel.parts)


def run_find_scan(path: Path, timeout_seconds: float) -> tuple[str, bytes, str]:
    command = [
        "timeout",
        f"{timeout_seconds:g}s",
        "find",
        str(path),
        "-mindepth",
        "1",
        "-maxdepth",
        "1",
        "-printf",
        "%y\\0%f\\0%p\\0",
    ]
    result = subprocess.run(command, capture_output=True, check=False)
    if result.returncode == 0:
        return "ok", result.stdout, result.stderr.decode("utf-8", errors="replace")
    if result.returncode == 124:
        return "timeout_blackhole", result.stdout, result.stderr.decode("utf-8", errors="replace")
    return "read_error", result.stdout, result.stderr.decode("utf-8", errors="replace")


def run_name_sample(path: Path, timeout_seconds: float, sample_limit: int) -> list[str]:
    command = [
        "timeout",
        f"{timeout_seconds:g}s",
        "find",
        str(path),
        "-mindepth",
        "1",
        "-maxdepth",
        "1",
        "-printf",
        "%f\\0",
    ]
    result = subprocess.run(command, capture_output=True, check=False)
    names = result.stdout.split(b"\0")
    sample: list[str] = []
    for raw_name in names:
        if not raw_name:
            continue
        sample.append(raw_name.decode("utf-8", errors="replace"))
        if len(sample) >= sample_limit:
            break
    return sample


def parse_find_records(stdout: bytes) -> list[tuple[str, str, str]]:
    parts = stdout.split(b"\0")
    if parts and not parts[-1]:
        parts = parts[:-1]

    records: list[tuple[str, str, str]] = []
    for i in range(0, len(parts) - 2, 3):
        entry_type = parts[i].decode("utf-8", errors="replace")
        name = parts[i + 1].decode("utf-8", errors="replace")
        full_path = parts[i + 2].decode("utf-8", errors="replace")
        records.append((entry_type, name, full_path))
    return records


def summarize_ok_record(
    path: Path,
    root: Path,
    records: list[tuple[str, str, str]],
    sample_limit: int,
    stderr_text: str,
) -> tuple[dict[str, Any], list[Path]]:
    type_counts = Counter(entry_type for entry_type, _name, _full_path in records)
    extension_counts = Counter(extension_for(name) for _entry_type, name, _full_path in records)
    example_names = [name for _entry_type, name, _full_path in records[:sample_limit]]
    sample_extensions = Counter(extension_for(name) for name in example_names)
    child_dirs = [Path(full_path) for entry_type, _name, full_path in records if entry_type == "d"]

    record: dict[str, Any] = {
        "path": str(path),
        "relative_path": str(path.relative_to(root)),
        "depth": depth_for(path, root),
        "status": "ok",
        "Nentries": len(records),
        "Ndirs": type_counts.get("d", 0),
        "Nfiles": type_counts.get("f", 0),
        "type_counts": dict(sorted(type_counts.items())),
        "sampleNames": example_names,
        "sampleExtensions": dict(sorted(sample_extensions.items())),
        "extension_counts": dict(sorted(extension_counts.items())),
    }
    if stderr_text:
        record["stderr"] = stderr_text.strip()
    return record, sorted(child_dirs, reverse=True)


def summarize_non_ok_record(
    path: Path,
    root: Path,
    status: str,
    stderr_text: str,
    sample_names: list[str],
) -> dict[str, Any]:
    sample_extensions = Counter(extension_for(name) for name in sample_names)
    record: dict[str, Any] = {
        "path": str(path),
        "relative_path": str(path.relative_to(root)),
        "depth": depth_for(path, root),
        "status": status,
        "Nentries": None,
        "Ndirs": None,
        "Nfiles": None,
        "sampleNames": sample_names,
        "sampleExtensions": dict(sorted(sample_extensions.items())),
    }
    if stderr_text:
        record["stderr"] = stderr_text.strip()
    return record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Walk a directory tree with per-directory timeouts.")
    parser.add_argument("root", help="Root directory to walk.")
    parser.add_argument("--jsonl-output", required=True, help="Output JSONL path.")
    parser.add_argument("--tsv-output", required=True, help="Output TSV path.")
    parser.add_argument("--scan-timeout", type=float, default=2.0, help="Seconds allowed for each direct-child scan.")
    parser.add_argument("--sample-timeout", type=float, default=2.0, help="Seconds allowed for timeout directory sampling.")
    parser.add_argument("--sample-limit", type=int, default=50, help="Maximum sample names stored per directory.")
    parser.add_argument("--progress-every", type=int, default=100, help="Print progress after this many directories.")
    parser.add_argument("--progress-seconds", type=float, default=10.0, help="Print progress after this many seconds.")
    return parser.parse_args()


def write_tsv_row(writer: csv.DictWriter, record: dict[str, Any]) -> None:
    writer.writerow(
        {
            "path": record["path"],
            "depth": record["depth"],
            "status": record["status"],
            "Nentries": "" if record["Nentries"] is None else record["Nentries"],
            "Ndirs": "" if record["Ndirs"] is None else record["Ndirs"],
            "Nfiles": "" if record["Nfiles"] is None else record["Nfiles"],
            "sampleNames": "; ".join(record["sampleNames"]),
            "sampleExtensions": json.dumps(record["sampleExtensions"], sort_keys=True),
        }
    )


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"Root is not a directory: {root}", file=sys.stderr)
        return 2
    if args.scan_timeout <= 0 or args.sample_timeout <= 0:
        print("--scan-timeout and --sample-timeout must be > 0", file=sys.stderr)
        return 2
    if args.sample_limit < 0:
        print("--sample-limit must be >= 0", file=sys.stderr)
        return 2

    jsonl_output = Path(args.jsonl_output)
    tsv_output = Path(args.tsv_output)
    jsonl_output.parent.mkdir(parents=True, exist_ok=True)
    tsv_output.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    status_counts: Counter[str] = Counter()
    max_depth_seen = 0
    stack = [root]
    started = time.monotonic()
    last_progress = started

    with jsonl_output.open("w", encoding="utf-8") as jsonl_handle, tsv_output.open(
        "w", encoding="utf-8", newline=""
    ) as tsv_handle:
        tsv_writer: csv.DictWriter = csv.DictWriter(
            tsv_handle,
            fieldnames=["path", "depth", "status", "Nentries", "Ndirs", "Nfiles", "sampleNames", "sampleExtensions"],
            delimiter="\t",
        )
        tsv_writer.writeheader()

        while stack:
            path = stack.pop()
            status, stdout, stderr_text = run_find_scan(path, args.scan_timeout)
            if status == "ok":
                find_records = parse_find_records(stdout)
                record, child_dirs = summarize_ok_record(path, root, find_records, args.sample_limit, stderr_text)
                stack.extend(child_dirs)
            else:
                sample_names = run_name_sample(path, args.sample_timeout, args.sample_limit)
                record = summarize_non_ok_record(path, root, status, stderr_text, sample_names)

            jsonl_handle.write(json.dumps(record, sort_keys=True) + "\n")
            jsonl_handle.flush()
            write_tsv_row(tsv_writer, record)
            tsv_handle.flush()

            count += 1
            status_counts[record["status"]] += 1
            max_depth_seen = max(max_depth_seen, record["depth"])

            now = time.monotonic()
            if (
                args.progress_every > 0
                and count % args.progress_every == 0
                or args.progress_seconds > 0
                and now - last_progress >= args.progress_seconds
            ):
                elapsed = now - started
                print(
                    "progress "
                    f"dirs={count} queued={len(stack)} max_depth={max_depth_seen} "
                    f"statuses={dict(status_counts)} elapsed={elapsed:.1f}s",
                    file=sys.stderr,
                    flush=True,
                )
                last_progress = now

    print(
        f"wrote {count} records to {jsonl_output} and {tsv_output}; statuses={dict(status_counts)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

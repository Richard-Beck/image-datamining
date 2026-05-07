from __future__ import annotations

import csv
import shlex
from pathlib import Path


def write_cellpose_command_manifest(
    prepared_manifest: Path,
    run_dir: Path,
    output: Path,
    mode: str,
    anisotropy: float,
    min_size: int,
    gpu: bool = True,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    masks_dir = run_dir / "masks"
    logs_dir = run_dir / "logs"
    masks_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    with prepared_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    fieldnames = ["sample_id", "mode", "input_path", "mask_path", "log_path", "command"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            sample_id = row["sample_id"]
            input_path = Path(row["prepared_input"])
            mask_path = masks_dir / f"{sample_id}_{mode}_cpsam_masks.tif"
            log_path = logs_dir / f"{sample_id}_{mode}_cpsam.log"
            command = [
                "python",
                "pipelines/fucci_3d_segmentation/bin/run_cellpose_cpsam.py",
                "--input",
                str(input_path),
                "--output",
                str(mask_path),
                "--mode",
                mode,
                "--anisotropy",
                str(anisotropy),
                "--min-size",
                str(min_size),
            ]
            if gpu:
                command.append("--gpu")
            shell_command = " ".join(shlex.quote(part) for part in command)
            shell_command = f"{shell_command} > {shlex.quote(str(log_path))} 2>&1"
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "mode": mode,
                    "input_path": str(input_path),
                    "mask_path": str(mask_path),
                    "log_path": str(log_path),
                    "command": shell_command,
                }
            )

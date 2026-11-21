#!/usr/bin/env python3
"""Create side-by-side video from rendered left/right image sequences."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run_ffmpeg(ffmpeg: str, args: list[str]) -> None:
    cmd = [ffmpeg, "-y", *args]
    completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if completed.returncode != 0:
        print(completed.stdout, end="")
        raise RuntimeError(f"ffmpeg failed: {' '.join(cmd)}")


def ensure_frames_exist(pattern: str) -> None:
    # Convert printf-style `%03d` to a wildcard and check at least one frame.
    wildcard = pattern.replace("%03d", "*").replace("%d", "*")
    if not list(Path(".").glob(wildcard)):
        raise FileNotFoundError(f"No input frames found for pattern: {pattern}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create SBS MP4 from rendered frame sequences.")
    parser.add_argument("--left-pattern", default="output_L_%03d.png", help="FFmpeg input pattern for left-eye frames.")
    parser.add_argument("--right-pattern", default="output_R_%03d.png", help="FFmpeg input pattern for right-eye frames.")
    parser.add_argument("--fps", type=int, default=30, help="Input/output frame rate.")
    parser.add_argument("--loops", type=int, default=8, help="How many times to repeat forward+reverse clip.")
    parser.add_argument("--output", default="sbs.mp4", help="Output SBS MP4 path.")
    parser.add_argument("--ffmpeg", default="ffmpeg", help="Path to ffmpeg executable.")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    ffmpeg = shutil.which(args.ffmpeg) if args.ffmpeg == "ffmpeg" else args.ffmpeg
    if not ffmpeg:
        print("ffmpeg not found in PATH. Install ffmpeg or pass --ffmpeg /path/to/ffmpeg.", file=sys.stderr)
        return 1

    if args.fps <= 0:
        print("--fps must be positive.", file=sys.stderr)
        return 1
    if args.loops <= 0:
        print("--loops must be positive.", file=sys.stderr)
        return 1

    try:
        ensure_frames_exist(args.left_pattern)
        ensure_frames_exist(args.right_pattern)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="smallpt4d_video_") as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        left = temp_dir / "left.mp4"
        right = temp_dir / "right.mp4"
        left_rev = temp_dir / "left_reverse.mp4"
        right_rev = temp_dir / "right_reverse.mp4"
        left_loop = temp_dir / "left_looped.mp4"
        right_loop = temp_dir / "right_looped.mp4"

        run_ffmpeg(
            ffmpeg,
            ["-framerate", str(args.fps), "-i", args.left_pattern, "-an", "-c:v", "libx264", "-pix_fmt", "yuv420p", str(left)],
        )
        run_ffmpeg(
            ffmpeg,
            ["-framerate", str(args.fps), "-i", args.right_pattern, "-an", "-c:v", "libx264", "-pix_fmt", "yuv420p", str(right)],
        )

        reverse_concat = "[0:v]split[f][r];[r]reverse[rev];[f][rev]concat=n=2:v=1[v]"
        run_ffmpeg(ffmpeg, ["-i", str(left), "-filter_complex", reverse_concat, "-map", "[v]", str(left_rev)])
        run_ffmpeg(ffmpeg, ["-i", str(right), "-filter_complex", reverse_concat, "-map", "[v]", str(right_rev)])

        # Repeat clip without extra concat files.
        run_ffmpeg(ffmpeg, ["-stream_loop", str(args.loops - 1), "-i", str(left_rev), "-c", "copy", str(left_loop)])
        run_ffmpeg(ffmpeg, ["-stream_loop", str(args.loops - 1), "-i", str(right_rev), "-c", "copy", str(right_loop)])

        run_ffmpeg(
            ffmpeg,
            [
                "-i",
                str(left_loop),
                "-i",
                str(right_loop),
                "-filter_complex",
                "[0:v][1:v]hstack=inputs=2[v]",
                "-map",
                "[v]",
                "-an",
                "-c:v",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                str(output_path),
            ],
        )

    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)

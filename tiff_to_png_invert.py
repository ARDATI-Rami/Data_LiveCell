#!/usr/bin/env python3
"""
Batch TIFF → PNG conversion with intensity inversion.

- Walks files and/or directories (optionally recursive).
- Supports uint8/uint16 grayscale, RGB, LA, RGBA; preserves alpha (not inverted).
- Writes next to inputs OR mirrors the directory tree into an output root.

Usage examples
--------------
Convert a single directory (non-recursive), write PNGs alongside inputs:
    python tiff_to_png_invert.py Real_wing_tissue_figs

Convert recursively and mirror into a clean output root:
    python tiff_to_png_invert.py Real_wing_tissue_figs --recursive -o Real_wing_tissue_figs_png

Convert only specific files:
    python tiff_to_png_invert.py Real_wing_tissue_figs/wingframe_001.tif Real_wing_tissue_figs/wingframe_100.tif
"""

from __future__ import annotations
import argparse
import pathlib
import sys
from typing import Iterable, Iterator, List, Sequence, Tuple

import numpy as np
from PIL import Image


TIFF_EXTS = {".tif", ".tiff"}


def _split_rgb_alpha(arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray | None]:
    """
    Split an array into color and alpha planes if present.

    :param arr: Image array with shape (H, W) or (H, W, C); C∈{2,4} may contain alpha.
    :type  arr: np.ndarray
    :return: (color_plane, alpha_plane or None). Color excludes the alpha channel if present.
    :rtype: Tuple[np.ndarray, Optional[np.ndarray]]
    """
    if arr.ndim == 3 and arr.shape[-1] in (2, 4):
        return arr[..., :-1], arr[..., -1]
    return arr, None


def _rejoin_rgb_alpha(color: np.ndarray, alpha: np.ndarray | None) -> np.ndarray:
    """
    Rejoin color and alpha planes.

    :param color: Color plane, shape (H, W) or (H, W, C_color).
    :type  color: np.ndarray
    :param alpha: Alpha plane or None.
    :type  alpha: Optional[np.ndarray]
    :return: Combined array with alpha appended if provided.
    :rtype: np.ndarray
    """
    if alpha is None:
        return color
    if color.ndim == 2:
        color = color[..., None]
    return np.concatenate([color, alpha[..., None]], axis=-1)


def invert_ndarray(arr: np.ndarray) -> np.ndarray:
    """
    Invert pixel intensities for uint8/uint16 images; preserve alpha channels.

    :param arr: Input image array. dtypes supported: uint8, uint16.
    :type  arr: np.ndarray
    :return: Inverted image array, same dtype/shape as input.
    :rtype: np.ndarray
    """
    if arr.dtype == np.uint8:
        maxv = np.uint8(255)
    elif arr.dtype == np.uint16:
        maxv = np.uint16(65535)
    else:
        # Fallback: convert to uint8, invert, return uint8
        arr = arr.astype(np.uint8, copy=False)
        maxv = np.uint8(255)

    color, alpha = _split_rgb_alpha(arr)
    inv_color = (maxv - color).astype(arr.dtype, copy=False)
    return _rejoin_rgb_alpha(inv_color, alpha)


def _iter_tiffs(paths: Sequence[pathlib.Path], recursive: bool) -> Iterator[pathlib.Path]:
    """
    Yield TIFF files from a list of paths (files or directories).

    :param paths: Input paths (files and/or directories).
    :type  paths: Sequence[pathlib.Path]
    :param recursive: If True, walk directories recursively.
    :type  recursive: bool
    :return: Iterator of TIFF file paths.
    :rtype: Iterator[pathlib.Path]
    """
    for p in paths:
        if p.is_file() and p.suffix.lower() in TIFF_EXTS:
            yield p
        elif p.is_dir():
            if recursive:
                for f in p.rglob("*"):
                    if f.is_file() and f.suffix.lower() in TIFF_EXTS:
                        yield f
            else:
                for f in p.glob("*"):
                    if f.is_file() and f.suffix.lower() in TIFF_EXTS:
                        yield f


def _decide_mode(inv: np.ndarray) -> str | None:
    """
    Choose a PIL mode for a given ndarray, when possible.

    :param inv: Image array after inversion.
    :type  inv: np.ndarray
    :return: PIL mode string or None if needs special handling (e.g., LA 16-bit).
    :rtype: Optional[str]
    """
    if inv.ndim == 2:
        return "I;16" if inv.dtype == np.uint16 else "L"
    if inv.ndim == 3 and inv.shape[-1] == 3:
        return "RGB"
    if inv.ndim == 3 and inv.shape[-1] == 4:
        return "RGBA"
    if inv.ndim == 3 and inv.shape[-1] == 2:
        # LA may require manual merge; return None to trigger special path
        return None
    raise ValueError(f"Unsupported shape {inv.shape}")


def process_one(in_path: pathlib.Path, out_root: pathlib.Path | None, base_root: pathlib.Path | None) -> pathlib.Path:
    """
    Process a single TIFF: invert intensities and write PNG.

    :param in_path: Input TIFF path.
    :type  in_path: pathlib.Path
    :param out_root: Output root directory; if None, write next to input.
    :type  out_root: Optional[pathlib.Path]
    :param base_root: If set, path relative to this root is preserved under out_root.
    :type  base_root: Optional[pathlib.Path]
    :return: Path to the written PNG file.
    :rtype: pathlib.Path
    """
    with Image.open(in_path) as im:
        im.load()
        # Preserve true bit depth when possible
        if im.mode in ("I;16", "I;16B", "I;16L"):
            arr = np.array(im, dtype=np.uint16)
        else:
            arr = np.array(im)

    inv = invert_ndarray(arr)

    mode = _decide_mode(inv)
    if mode is None:
        # Handle LA explicitly (8- or 16-bit)
        L = Image.fromarray(inv[..., 0], "I;16" if inv.dtype == np.uint16 else "L")
        A = Image.fromarray(inv[..., 1], "I;16" if inv.dtype == np.uint16 else "L")
        out_img = Image.merge("LA", (L, A))
    else:
        out_img = Image.fromarray(inv, mode)

    # Decide output path
    if out_root is None:
        target = in_path.with_suffix(".png")
    else:
        if base_root and in_path.is_relative_to(base_root):
            rel = in_path.relative_to(base_root).with_suffix(".png")
            target = out_root / rel
        else:
            target = out_root / in_path.name.replace(in_path.suffix, ".png")
        target.parent.mkdir(parents=True, exist_ok=True)

    out_img.save(target)
    return target


def main(argv: Sequence[str] | None = None) -> int:
    """
    CLI entrypoint.

    :param argv: Optional command-line arguments (defaults to sys.argv[1:]).
    :type  argv: Optional[Sequence[str]]
    :return: Process exit code (0 on success, non-zero on failure).
    :rtype: int
    """
    ap = argparse.ArgumentParser(description="Convert TIFF(s) to PNG with inverted intensities.")
    ap.add_argument("inputs", nargs="+",
                    help="Files and/or directories. Use --recursive to walk subdirectories.")
    ap.add_argument("-o", "--out-dir", type=pathlib.Path, default=None,
                    help="Output root directory. If omitted, PNGs are written next to inputs.")
    ap.add_argument("-r", "--recursive", action="store_true",
                    help="Recurse into subdirectories of any directory inputs.")
    ap.add_argument("--mirror-structure", action="store_true",
                    help="When using --out-dir, mirror the input directory structure under it. "
                         "If not set, files go flat into --out-dir.")
    args = ap.parse_args(argv)

    in_paths: List[pathlib.Path] = [pathlib.Path(p) for p in args.inputs]
    if args.out_dir:
        args.out_dir.mkdir(parents=True, exist_ok=True)

    # If exactly one directory provided, we can treat it as the base to mirror
    base_root: pathlib.Path | None = None
    if args.out_dir and args.mirror_structure:
        # Use the first directory among inputs as base_root; if multiple, choose their
        # common parent so structure is preserved sensibly.
        dirs = [p for p in in_paths if p.is_dir()]
        if len(dirs) == 1:
            base_root = dirs[0].resolve()
        elif len(dirs) > 1:
            # Find common parent
            base_root = pathlib.Path(*pathlib.Path(*[d.resolve() for d in dirs][0].parts)).anchor
            # Fallback: use first dir to keep behavior deterministic
            base_root = dirs[0].resolve()

    tiff_files = sorted(set(_iter_tiffs(in_paths, recursive=args.recursive)))
    if not tiff_files:
        print("No TIFF files found.", file=sys.stderr)
        return 2

    n_ok = 0
    for f in tiff_files:
        try:
            written = process_one(f, args.out_dir, base_root if args.mirror_structure else None)
            print(f"[OK] {f} → {written}")
            n_ok += 1
        except Exception as e:
            print(f"[FAIL] {f}: {e}", file=sys.stderr)

    print(f"Done: {n_ok}/{len(tiff_files)} converted.")
    return 0 if n_ok == len(tiff_files) else 1


if __name__ == "__main__":
    raise SystemExit(main())

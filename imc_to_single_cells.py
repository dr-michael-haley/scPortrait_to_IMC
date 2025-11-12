#!/usr/bin/env python3
"""
IMC → scPortrait helper.

This script discovers Imaging Mass Cytometry ROIs, loads the per-channel TIFFs,
attaches a single whole-cell segmentation mask, and runs scPortrait's
HDF5CellExtraction to generate single-cell image collections (.h5sc).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
from tifffile import imread
from skimage.segmentation import expand_labels

from scportrait.pipeline._utils.constants import DEFAULT_SEGMENTATION_DTYPE
from scportrait.pipeline.extraction import HDF5CellExtraction
from scportrait.pipeline.project import Project


CHANNEL_EXT_DEFAULT = [".tif", ".tiff", ".ome.tif", ".ome.tiff"]


def _normalize_ext(exts: Iterable[str]) -> set[str]:
    normalized: set[str] = set()
    for ext in exts:
        ext = ext if ext.startswith(".") else f".{ext}"
        normalized.add(ext.lower())
        normalized.add(ext.upper())
    return normalized


def _list_roi_dirs(images_root: Path) -> list[Path]:
    if not images_root.is_dir():
        raise ValueError(f"Images directory {images_root} not found.")
    roi_dirs = sorted([p for p in images_root.iterdir() if p.is_dir()])
    if not roi_dirs:
        raise ValueError(f"No ROI subdirectories found in {images_root}.")
    return roi_dirs


def _stems_in_dir(directory: Path, extensions: set[str]) -> set[str]:
    stems: set[str] = set()
    for ext in extensions:
        for path in directory.glob(f"*{ext}"):
            stems.add(path.stem)
    return stems


def _collect_rois(images_root: Path, mask_dir: Path, extensions: set[str]) -> list[str]:
    if not mask_dir.is_dir():
        raise ValueError(f"Mask directory {mask_dir} not found.")
    roi_dirs = _list_roi_dirs(images_root)
    roi_names = {d.name for d in roi_dirs}
    mask_ids = _stems_in_dir(mask_dir, extensions)
    roi_ids = sorted(roi_names & mask_ids)
    if not roi_ids:
        raise ValueError("No overlapping ROI identifiers between ROI folders and mask directory.")
    return roi_ids


def _resolve_file(directory: Path, roi: str, extensions: set[str]) -> Path:
    for ext in extensions:
        candidate = directory / f"{roi}{ext}"
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No file for ROI '{roi}' in {directory} (extensions {sorted(extensions)}).")


def _spatial_shape(shape: tuple[int, ...]) -> tuple[int, int]:
    if len(shape) == 2:
        return shape[0], shape[1]
    if len(shape) == 3:
        if shape[0] == 1:
            return shape[1], shape[2]
        if shape[-1] == 1:
            return shape[0], shape[1]
    raise ValueError(f"Cannot infer 2D plane from shape {shape}.")


def _load_mask(mask_path: Path, expected_shape: tuple[int, int], expand_px: int = 0) -> np.ndarray:
    mask = imread(mask_path)
    spatial_shape = _spatial_shape(mask.shape)
    if spatial_shape != expected_shape:
        raise ValueError(
            f"Mask {mask_path} has shape {spatial_shape}, expected {expected_shape} based on channel data."
        )
    mask_data = np.squeeze(mask)
    if expand_px > 0:
        mask_data = expand_labels(mask_data.astype(np.int32, copy=False), distance=expand_px)
    return mask_data.astype(DEFAULT_SEGMENTATION_DTYPE, copy=False)


def _channel_files_and_names(roi_dir: Path, extensions: set[str]) -> tuple[list[str], list[str]]:
    tiffs = sorted(
        [
            path
            for ext in extensions
            for path in roi_dir.glob(f"*{ext}")
            if path.is_file()
        ]
    )
    if not tiffs:
        raise ValueError(f"No channel images found in {roi_dir}.")
    files = [str(path) for path in tiffs]
    names = [path.stem for path in tiffs]
    return files, names


def _run_single_roi(
    roi: str,
    project_root: Path,
    config_path: Path,
    roi_dir: Path,
    mask_dir: Path,
    extensions: set[str],
    overwrite: bool,
    debug: bool,
    mask_expand_px: int,
) -> Path:
    project_dir = project_root / roi
    project_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = project_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    channel_files, channel_names = _channel_files_and_names(roi_dir, extensions)
    example_image = imread(channel_files[0])
    spatial_shape = _spatial_shape(example_image.shape)
    mask_path = _resolve_file(mask_dir, roi, extensions)
    mask = _load_mask(mask_path, spatial_shape, expand_px=mask_expand_px)

    project = Project(
        project_location=str(project_dir),
        config_path=str(config_path),
        extraction_f=HDF5CellExtraction,
        overwrite=overwrite,
        debug=debug,
    )
    project.load_input_from_tif_files(
        channel_files,
        channel_names=channel_names,
        overwrite=True,
        cache=str(cache_dir),
    )

    project.filehandler._write_segmentation_sdata(
        mask,
        project.cyto_seg_name,
        chunks=project.DEFAULT_CHUNK_SIZE_2D,
        overwrite=True,
    )
    project.filehandler._add_centers(project.cyto_seg_name, overwrite=True)
    project.extraction_f.register_parameter("segmentation_mask", project.cyto_seg_name)

    project.extract(overwrite=overwrite)

    output_file = (
        Path(project.project_location)
        / project.DEFAULT_EXTRACTION_DIR_NAME
        / project.DEFAULT_DATA_DIR
        / project.DEFAULT_EXTRACTION_FILE
    )
    return output_file


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Batch IMC extraction with scPortrait (single-mask workflow).")
    parser.add_argument("--channels-dir", required=True, type=Path, help="Folder with one subdirectory per IMC channel.")
    parser.add_argument("--mask-dir", required=True, type=Path, help="Folder containing whole-cell mask TIFFs.")
    parser.add_argument(
        "--projects-root",
        required=True,
        type=Path,
        help="Directory where per-ROI scPortrait projects will be created.",
    )
    parser.add_argument("--config", required=True, type=Path, help="scPortrait config.yml containing HDF5CellExtraction.")
    parser.add_argument(
        "--roi",
        action="append",
        dest="roi_list",
        help="Specific ROI identifier to process (repeat for multiple). Default: run all ROIs discovered from masks.",
    )
    parser.add_argument(
        "--image-ext",
        nargs="+",
        default=CHANNEL_EXT_DEFAULT,
        help=f"Extensions for both channel images and masks (default: {CHANNEL_EXT_DEFAULT}).",
    )
    parser.add_argument(
        "--mask-expand-px",
        type=int,
        default=0,
        help="Expand each labelled cell mask by the given number of pixels (uses skimage.segmentation.expand_labels).",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite any existing project directories.")
    parser.add_argument("--debug", action="store_true", help="Enable verbose project logging.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    extensions = _normalize_ext(args.image_ext)
    roi_dirs = {roi_dir.name: roi_dir for roi_dir in _list_roi_dirs(args.channels_dir)}
    roi_ids = args.roi_list or _collect_rois(args.channels_dir, args.mask_dir, extensions)

    args.projects_root.mkdir(parents=True, exist_ok=True)
    for roi in roi_ids:
        try:
            if roi not in roi_dirs:
                raise ValueError(f"Roi folder '{roi}' not found below {args.channels_dir}.")
            output = _run_single_roi(
                roi=roi,
                project_root=args.projects_root,
                config_path=args.config,
                roi_dir=roi_dirs[roi],
                mask_dir=args.mask_dir,
                extensions=extensions,
                overwrite=args.overwrite,
                debug=args.debug,
                mask_expand_px=args.mask_expand_px,
            )
            print(f"[{roi}] single-cell collection written to {output}")
        except Exception as exc:  # pragma: no cover - user feedback path
            print(f"[{roi}] failed: {exc}", file=sys.stderr)
            if not args.debug:
                print("Re-run with --debug for more details.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

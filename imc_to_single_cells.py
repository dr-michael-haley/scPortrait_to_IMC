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
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
from uuid import uuid4

import numpy as np
from tifffile import imread, imwrite
from skimage.segmentation import expand_labels

from scportrait.io.h5sc import read_h5sc
from scportrait.pipeline._utils.constants import (
    DEFAULT_CELL_ID_NAME,
    DEFAULT_IMAGE_DTYPE,
    DEFAULT_SEGMENTATION_DTYPE,
)
from scportrait.pipeline.extraction import HDF5CellExtraction
from scportrait.pipeline.project import Project


CHANNEL_EXT_DEFAULT = [".tif", ".tiff", ".ome.tif", ".ome.tiff"]
ROI_VERTICAL_SPACING = 32


@dataclass
class ROIContext:
    name: str
    channel_files: list[Path]
    mask_path: Path
    spatial_shape: tuple[int, int]


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


def _derive_channel_name(path: Path) -> str:
    """
    Extract human-readable channel name from filenames formatted as
    {chan#}_{roi#}_{str}_{channel name}.tiff
    """
    stem = path.stem
    parts = stem.split("_", 3)
    if len(parts) == 4:
        return parts[3]
    return stem


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
    names = [_derive_channel_name(path) for path in tiffs]
    return files, names


def _safe_channel_filename(index: int, name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    if not slug:
        slug = f"channel_{index:03d}"
    return f"{index:03d}_{slug}.tif"


def _compute_canvas_layout(contexts: Sequence[ROIContext], spacing: int = ROI_VERTICAL_SPACING) -> tuple[int, int, dict[str, tuple[int, int]]]:
    if not contexts:
        raise ValueError("No ROI contexts available to build canvas layout.")
    max_width = max(ctx.spatial_shape[1] for ctx in contexts)
    total_height = sum(ctx.spatial_shape[0] for ctx in contexts)
    if len(contexts) > 1:
        total_height += spacing * (len(contexts) - 1)
    offsets: dict[str, tuple[int, int]] = {}
    cursor = 0
    for ctx in contexts:
        offsets[ctx.name] = (cursor, cursor + ctx.spatial_shape[0])
        cursor += ctx.spatial_shape[0] + spacing
    return total_height, max_width, offsets


def _write_combined_channel_images(
    contexts: Sequence[ROIContext],
    channel_names: list[str],
    offsets: dict[str, tuple[int, int]],
    assembled_dir: Path,
    canvas_height: int,
    canvas_width: int,
) -> list[str]:
    assembled_dir.mkdir(parents=True, exist_ok=True)
    channel_paths: list[str] = []
    for idx, channel_name in enumerate(channel_names):
        canvas = np.zeros((canvas_height, canvas_width), dtype=DEFAULT_IMAGE_DTYPE)
        for ctx in contexts:
            start_y, end_y = offsets[ctx.name]
            tile = imread(ctx.channel_files[idx])
            tile = np.squeeze(tile)
            tile = np.asarray(tile, dtype=DEFAULT_IMAGE_DTYPE, order="C")
            height, width = tile.shape
            if height != ctx.spatial_shape[0] or width != ctx.spatial_shape[1]:
                raise ValueError(
                    f"Channel '{channel_name}' for ROI '{ctx.name}' has unexpected shape {(height, width)}, expected {ctx.spatial_shape}."
                )
            canvas[start_y:end_y, 0:width] = tile
        filename = _safe_channel_filename(idx, channel_name)
        out_path = assembled_dir / filename
        imwrite(out_path, canvas)
        channel_paths.append(str(out_path))
    return channel_paths


def _assemble_combined_mask(
    contexts: Sequence[ROIContext],
    offsets: dict[str, tuple[int, int]],
    canvas_height: int,
    canvas_width: int,
    mask_expand_px: int,
) -> tuple[np.ndarray, dict[str, tuple[int, int]]]:
    combined = np.zeros((canvas_height, canvas_width), dtype=DEFAULT_SEGMENTATION_DTYPE)
    roi_cell_ranges: dict[str, tuple[int, int]] = {}
    next_offset = 0
    for ctx in contexts:
        start_y, end_y = offsets[ctx.name]
        mask = _load_mask(ctx.mask_path, ctx.spatial_shape, expand_px=mask_expand_px)
        mask = np.asarray(mask, dtype=DEFAULT_SEGMENTATION_DTYPE, order="C")
        non_zero = mask > 0
        if np.any(non_zero):
            roi_values = mask[non_zero]
            roi_min = int(roi_values.min())
            roi_max = int(roi_values.max())
            range_start = next_offset + roi_min
            range_end = next_offset + roi_max
            roi_cell_ranges[ctx.name] = (range_start, range_end)
            mask[non_zero] = roi_values + next_offset
            next_offset = range_end
        else:
            roi_cell_ranges[ctx.name] = (0, 0)
        combined[start_y:end_y, 0 : mask.shape[1]] = mask
    return combined, roi_cell_ranges


def _annotate_roi_membership(output_file: Path, roi_ranges: dict[str, tuple[int, int]]) -> None:
    if not roi_ranges:
        return
    try:
        adata = read_h5sc(str(output_file))
    except Exception as exc:  # pragma: no cover - runtime safeguard
        print(f"Failed to annotate ROI membership: {exc}", file=sys.stderr)
        return

    if DEFAULT_CELL_ID_NAME not in adata.obs:
        return

    cell_ids = adata.obs[DEFAULT_CELL_ID_NAME].to_numpy()
    roi_labels = np.full(cell_ids.shape, "unassigned", dtype=object)
    for roi, (start_id, end_id) in roi_ranges.items():
        if start_id == 0 and end_id == 0:
            continue
        mask = (cell_ids >= start_id) & (cell_ids <= end_id)
        roi_labels[mask] = roi

    adata.obs["roi"] = roi_labels
    adata.write_h5ad(output_file)


def _ensure_cache_directory(cache_value: str | None) -> None:
    if not cache_value:
        return
    cache_path = Path(cache_value)
    if cache_path.is_absolute():
        target = cache_path
    else:
        tmp_root = Path(tempfile.gettempdir())
        target = tmp_root / cache_path
    target.mkdir(parents=True, exist_ok=True)


def _run_combined_project(
    project_dir: Path,
    config_path: Path,
    channel_files: list[str],
    channel_names: list[str],
    mask: np.ndarray,
    overwrite: bool,
    debug: bool,
) -> Path:
    project_dir.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(tempfile.gettempdir())
    cache_token = f"scportrait_imc_cache_{uuid4().hex[:8]}"
    cache_dir = tmp_root / cache_token
    cache_dir.mkdir(parents=True, exist_ok=True)

    project = Project(
        project_location=str(project_dir),
        config_path=str(config_path),
        extraction_f=HDF5CellExtraction,
        overwrite=overwrite,
        debug=debug,
    )
    try:
        project.load_input_from_tif_files(
            channel_files,
            channel_names=channel_names,
            overwrite=True,
            cache=cache_token,
        )
        extraction_cfg = project.config.get("HDF5CellExtraction", {})
        _ensure_cache_directory(extraction_cfg.get("cache"))

        project.filehandler._write_segmentation_sdata(
            mask,
            project.cyto_seg_name,
            chunks=project.DEFAULT_CHUNK_SIZE_2D,
            overwrite=True,
        )
        project.filehandler._add_centers(project.cyto_seg_name, overwrite=True)
        project.extraction_f.register_parameter("segmentation_mask", project.cyto_seg_name)

        project.extract(overwrite=overwrite)
    finally:
        shutil.rmtree(cache_dir, ignore_errors=True)

    output_file = (
        Path(project.project_location)
        / project.DEFAULT_EXTRACTION_DIR_NAME
        / project.DEFAULT_DATA_DIR
        / project.DEFAULT_EXTRACTION_FILE
    )
    return output_file


def process_imc_rois(
    channels_dir: str | Path,
    mask_dir: str | Path,
    projects_root: str | Path,
    config_path: str | Path,
    roi_list: Sequence[str] | None = None,
    image_ext: Iterable[str] | None = None,
    mask_expand_px: int = 0,
    overwrite: bool = False,
    debug: bool = False,
) -> tuple[dict[str, Path], dict[str, Exception]]:
    """
    Process IMC ROIs and extract single-cell images via scPortrait using a single combined project.

    Returns:
        (successes, failures) where `successes` maps each successfully processed ROI
        to the shared output h5sc Path and `failures` maps ROI -> raised Exception.
    """
    channels_dir = Path(channels_dir)
    mask_dir = Path(mask_dir)
    projects_root = Path(projects_root)
    config_path = Path(config_path)

    extensions = _normalize_ext(image_ext or CHANNEL_EXT_DEFAULT)
    roi_dirs = {roi_dir.name: roi_dir for roi_dir in _list_roi_dirs(channels_dir)}
    roi_ids = roi_list or _collect_rois(channels_dir, mask_dir, extensions)

    projects_root.mkdir(parents=True, exist_ok=True)

    successes: dict[str, Path] = {}
    failures: dict[str, Exception] = {}

    contexts: list[ROIContext] = []
    channel_names_template: list[str] | None = None

    for roi in roi_ids:
        if roi not in roi_dirs:
            failures[roi] = KeyError(f"ROI folder '{roi}' not found under {channels_dir}.")
            continue
        try:
            channel_files, channel_names = _channel_files_and_names(roi_dirs[roi], extensions)
            example_image = imread(channel_files[0])
            spatial_shape = _spatial_shape(example_image.shape)
            mask_path = _resolve_file(mask_dir, roi, extensions)
        except Exception as exc:
            failures[roi] = exc
            continue

        if channel_names_template is None:
            channel_names_template = channel_names
        elif channel_names_template != channel_names:
            failures[roi] = ValueError(f"Channel list for ROI '{roi}' does not match the first ROI.")
            continue

        contexts.append(
            ROIContext(
                name=roi,
                channel_files=[Path(p) for p in channel_files],
                mask_path=mask_path,
                spatial_shape=spatial_shape,
            )
        )

    if not contexts:
        return successes, failures
    if channel_names_template is None:
        raise RuntimeError("Unable to determine channel names for the combined project.")

    assembled_dir = projects_root / "assembled_channels"
    if assembled_dir.exists():
        shutil.rmtree(assembled_dir, ignore_errors=True)

    canvas_height, canvas_width, offsets = _compute_canvas_layout(contexts)
    aggregated_channel_files = _write_combined_channel_images(
        contexts,
        channel_names_template or [],
        offsets,
        assembled_dir,
        canvas_height,
        canvas_width,
    )
    combined_mask, roi_ranges = _assemble_combined_mask(
        contexts,
        offsets,
        canvas_height,
        canvas_width,
        mask_expand_px,
    )

    try:
        output = _run_combined_project(
            project_dir=projects_root,
            config_path=config_path,
            channel_files=aggregated_channel_files,
            channel_names=channel_names_template or [],
            mask=combined_mask,
            overwrite=overwrite,
            debug=debug,
        )
    finally:
        shutil.rmtree(assembled_dir, ignore_errors=True)

    _annotate_roi_membership(output, roi_ranges)

    for ctx in contexts:
        successes[ctx.name] = output
    return successes, failures


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
    outputs, errors = process_imc_rois(
        channels_dir=args.channels_dir,
        mask_dir=args.mask_dir,
        projects_root=args.projects_root,
        config_path=args.config,
        roi_list=args.roi_list,
        image_ext=args.image_ext,
        mask_expand_px=args.mask_expand_px,
        overwrite=args.overwrite,
        debug=args.debug,
    )

    for roi, path in outputs.items():
        print(f"[{roi}] single-cell collection written to {path}")

    if errors:
        for roi, exc in errors.items():
            print(f"[{roi}] failed: {exc}", file=sys.stderr)
        if not args.debug:
            print("Re-run with --debug for more details.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import ConvexHull, QhullError


OUTPUT_COLUMNS = [
    "FILENAME",
    "PANEL_SIZE_TOTAL_TARGETS",
    "QC_FEATURE_NUMBER",
    "QC_MEAN_READS_PER_FEATURE",
    "QC_TOTAL_GENES_DETECTED",
    "QC_TOTAL_NUMBER_OF_READS",
    "REGION_AREA",
]

COORDINATE_COLUMN_CANDIDATES = [
    ("CenterX_global_px", "CenterY_global_px"),
    ("CenterX_global", "CenterY_global"),
    ("CenterX", "CenterY"),
    ("x", "y"),
    ("X", "Y"),
    ("x_coord", "y_coord"),
    ("x_coordinate", "y_coordinate"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize one or more AnnData count files ending in 'counts.h5ad' "
            "and write the requested QC metrics to a CSV table."
        )
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing one or more files ending in 'counts.h5ad'.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("output.csv"),
        help="Output CSV path. Default: output.csv",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Number of cells to process at a time when summing counts. Default: 10000",
    )
    parser.add_argument(
        "--region-area-method",
        choices=("convex_hull", "bounding_box"),
        default="convex_hull",
        help=(
            "How to estimate REGION_AREA from cell x/y coordinates. "
            "Default: convex_hull"
        ),
    )
    return parser.parse_args()


def iter_count_files(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    return sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.name.endswith("counts.h5ad")
    )


def chunk_sum(chunk: object) -> float:
    if sparse.issparse(chunk):
        return float(chunk.sum())
    return float(np.asarray(chunk).sum(dtype=np.float64))


def total_reads(adata: ad.AnnData, chunk_size: int) -> int:
    total = 0.0
    for chunk, _, _ in adata.chunked_X(chunk_size):
        total += chunk_sum(chunk)
    return int(round(total))


def unique_valid_points(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError("Expected a 2D array with at least two coordinate columns.")

    points = points[:, :2]
    points = points[~np.isnan(points).any(axis=1)]
    if len(points) == 0:
        return points
    return np.unique(points, axis=0)


def get_spatial_points(adata: ad.AnnData) -> np.ndarray:
    for x_col, y_col in COORDINATE_COLUMN_CANDIDATES:
        if x_col in adata.obs.columns and y_col in adata.obs.columns:
            return np.column_stack(
                [
                    adata.obs[x_col].to_numpy(dtype=float),
                    adata.obs[y_col].to_numpy(dtype=float),
                ]
            )

    if "spatial" in adata.obsm:
        return np.asarray(adata.obsm["spatial"], dtype=float)

    raise KeyError(
        "Could not find cell x/y coordinates in adata.obs or adata.obsm['spatial']."
    )


def region_area(adata: ad.AnnData, method: str) -> float:
    points = unique_valid_points(get_spatial_points(adata))
    if len(points) < 2:
        return 0.0

    x_span = float(points[:, 0].max() - points[:, 0].min())
    y_span = float(points[:, 1].max() - points[:, 1].min())
    bbox_area = x_span * y_span

    if method == "bounding_box" or len(points) < 3:
        return bbox_area

    try:
        hull = ConvexHull(points)
        return float(hull.volume)
    except QhullError:
        return bbox_area


def summarize_file(path: Path, chunk_size: int, area_method: str) -> dict[str, object]:
    adata = ad.read_h5ad(path, backed="r")
    try:
        n_genes = int(adata.n_vars)
        n_cells = int(adata.n_obs)
        read_total = total_reads(adata, chunk_size=chunk_size)
        mean_reads = float(read_total / n_cells) if n_cells else 0.0
        area = float(region_area(adata, method=area_method))

        return {
            "FILENAME": path.name,
            "PANEL_SIZE_TOTAL_TARGETS": n_genes,
            "QC_FEATURE_NUMBER": n_cells,
            "QC_MEAN_READS_PER_FEATURE": mean_reads,
            "QC_TOTAL_GENES_DETECTED": n_genes,
            "QC_TOTAL_NUMBER_OF_READS": read_total,
            "REGION_AREA": area,
        }
    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()


def summarize_files(
    paths: Iterable[Path], chunk_size: int, area_method: str
) -> pd.DataFrame:
    rows = [summarize_file(path, chunk_size, area_method) for path in paths]
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def main() -> None:
    args = parse_args()
    paths = iter_count_files(args.input_dir)
    summary_df = summarize_files(
        paths,
        chunk_size=args.chunk_size,
        area_method=args.region_area_method,
    )
    summary_df.to_csv(args.output, index=False)
    print(f"Wrote {len(summary_df)} row(s) to {args.output}")


if __name__ == "__main__":
    main()

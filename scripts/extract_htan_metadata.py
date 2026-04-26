#!/usr/bin/env python3
"""Fill HTAN Xenium Level 3 metadata from Xenium output bundles.

Usage:
    python extract_htan_metadata.py config.json

The script reads a metadata CSV, inspects Xenium output directories (or tarballs),
fills missing required Level 3 fields, normalizes a few spec-bound enum values,
and writes the updated metadata plus panel-support outputs into `outputs/`.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import tarfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REQUIRED_FIELDS_BASE = [
    "PLATFORM",
    "ASSAY_CHEMISTRY_VERSION",
    "RNA_MEASURED",
    "PROTEIN_MEASURED",
    "PANEL_SIZE_TOTAL_TARGETS",
    "REGION_AREA",
    "BUNDLE_CONTENTS",
    "HAS_CELL_SEGMENTATION",
    "HAS_CLUSTERING",
    "QC_SPATIAL_UNIT",
    "QC_FEATURE_NUMBER",
    "QC_MEAN_READS_PER_FEATURE",
    "QC_TOTAL_GENES_DETECTED",
    "QC_TOTAL_NUMBER_OF_READS",
    "FILENAME",
    "FILE_FORMAT",
    "HTAN_DATA_FILE_ID",
    "HTAN_PARENT_ID",
]

ORDERED_BUNDLE_CONTENTS = [
    "morphology.ome.tif",
    "analysis.zarr.zip",
    "cell_boundaries.parquet",
    "cell_feature_matrix.zarr.zip",
    "cells.parquet",
    "cells.zarr.zip",
    "experiment.xenium",
    "nucleus_boundaries.parquet",
    "transcripts.parquet",
    "transcripts.zarr.zip",
    "gene_panel.json",
]

AUTHORITATIVE_FIELDS = {
    "ASSAY_CHEMISTRY_VERSION",
    "BUNDLE_CONTENTS",
    "DIMENSIONALITY_REDUCTION_METHOD",
    "FILE_FORMAT",
    "FILENAME",
    "HAS_DIMENSIONALITY_REDUCTION",
    "NUMBER_OF_CLUSTERS",
    "NUMBER_OF_SEGMENTED_CELLS",
    "PLATFORM",
    "RUN_ID",
    "SLIDE_SERIAL_NUMBER",
    "SOFTWARE_AND_VERSION",
    "SPATIAL_ASSAY_TYPE",
}

BOOL_TRUE = {"true", "t", "1", "yes", "y"}
BOOL_FALSE = {"false", "f", "0", "no", "n"}
VALID_HTAN_GENE_ID = re.compile(r"(ENSG\d+|\d+)$")


def is_blank(value: object) -> bool:
    return value is None or str(value).strip() == ""


def as_clean_string(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def parse_bool(value: object) -> bool | None:
    text = as_clean_string(value).lower()
    if text in BOOL_TRUE:
        return True
    if text in BOOL_FALSE:
        return False
    return None


def normalize_panel_key(value: str) -> str:
    text = as_clean_string(value)
    text = text.replace("/", "_")
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text.upper()


def safe_panel_id(value: str) -> str:
    text = as_clean_string(value)
    text = text.replace("/", "_")
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "panel"


def csv_read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = []
        for row in reader:
            cleaned = {key: as_clean_string(value) for key, value in row.items()}
            rows.append(cleaned)
    return rows, fieldnames


def csv_write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def resolve_path(path_text: str, base_dir: Path) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


@dataclass
class PanelDefinition:
    display_name: str
    canonical_key: str
    panel_type: str
    predesigned_panel: str
    description: str
    n_gene: str
    n_gene_total: str


@dataclass
class PanelReference:
    definitions_by_key: dict[str, PanelDefinition]
    specimen_to_panel_key: dict[str, str]


def parse_panel_reference(path: Path) -> PanelReference:
    definitions: dict[str, PanelDefinition] = {}
    specimen_to_panel_key: dict[str, str] = {}
    panel_type_values = {"Predesigned", "Custom"}

    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        rows = [list(row) for row in reader]

    for row in rows:
        padded = row + [""] * (6 - len(row))
        c0, c1, c2, c3, c4, c5 = [as_clean_string(value) for value in padded[:6]]
        if c0 and c1 in panel_type_values:
            canonical_key = normalize_panel_key(c0)
            definitions[canonical_key] = PanelDefinition(
                display_name=c0,
                canonical_key=canonical_key,
                panel_type=c1,
                predesigned_panel=c2,
                description=c3,
                n_gene=c4,
                n_gene_total=c5,
            )

    for row in rows:
        padded = row + [""] * (6 - len(row))
        c0, c1, c2, c3 = [as_clean_string(value) for value in padded[:4]]
        panel_key = normalize_panel_key(c3)
        if (
            c0
            and c3
            and panel_key in definitions
            and c2.lower() in {"breast", "prostate"}
            and not c1.startswith("ENSG")
            and not c1.startswith("ENST")
        ):
            specimen_to_panel_key[c0] = panel_key

    return PanelReference(
        definitions_by_key=definitions,
        specimen_to_panel_key=specimen_to_panel_key,
    )


class BundleReader:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.is_dir = path.is_dir()
        if self.is_dir:
            self._tar = None
            self._members = []
        else:
            self._tar = tarfile.open(path, "r:*")
            self._members = [member for member in self._tar.getmembers() if member.isfile()]

    def close(self) -> None:
        if self._tar is not None:
            self._tar.close()

    def _find_member(self, relative_path: str):
        assert self._tar is not None
        wanted = relative_path.lstrip("./")
        matches = []
        for member in self._members:
            name = member.name.lstrip("./")
            if name == wanted or name.endswith(f"/{wanted}"):
                matches.append(member)
        if not matches:
            return None
        matches.sort(key=lambda member: len(member.name))
        return matches[0]

    def exists(self, relative_path: str) -> bool:
        if self.is_dir:
            return (self.path / relative_path).exists()
        return self._find_member(relative_path) is not None

    def read_text(self, relative_path: str) -> str:
        if self.is_dir:
            return (self.path / relative_path).read_text()
        member = self._find_member(relative_path)
        if member is None:
            raise FileNotFoundError(relative_path)
        extracted = self._tar.extractfile(member)
        if extracted is None:
            raise FileNotFoundError(relative_path)
        return extracted.read().decode()

    def list_root_entries(self) -> list[str]:
        if self.is_dir:
            return sorted(child.name for child in self.path.iterdir())

        roots = set()
        for member in self._members:
            parts = Path(member.name.lstrip("./")).parts
            if not parts:
                continue
            roots.add(parts[0])

        if len(roots) == 1:
            root_prefix = next(iter(roots))
            second_level = set()
            for member in self._members:
                parts = Path(member.name.lstrip("./")).parts
                if len(parts) >= 2 and parts[0] == root_prefix:
                    second_level.add(parts[1])
            return sorted(second_level)

        return sorted(roots)

    def find_root_html(self) -> str:
        if self.is_dir:
            html_paths = sorted(child.name for child in self.path.iterdir() if child.suffix.lower() == ".html")
            if "analysis_summary.html" in html_paths:
                return "analysis_summary.html"
            return html_paths[0] if html_paths else ""

        html_names: list[str] = []
        for member in self._members:
            name = member.name.lstrip("./")
            parts = Path(name).parts
            if not parts:
                continue
            rel_name = parts[-1] if len(self.list_root_entries()) != 1 else "/".join(parts[1:])
            if rel_name.endswith(".html") and "/" not in rel_name:
                html_names.append(rel_name)
        html_names = sorted(set(html_names))
        if "analysis_summary.html" in html_names:
            return "analysis_summary.html"
        return html_names[0] if html_names else ""


def read_json(bundle: BundleReader, relative_path: str) -> dict:
    return json.loads(bundle.read_text(relative_path))


def read_single_row_csv(bundle: BundleReader, relative_path: str) -> dict[str, str]:
    text = bundle.read_text(relative_path)
    reader = csv.DictReader(io.StringIO(text))
    row = next(reader, None)
    return {key: as_clean_string(value) for key, value in (row or {}).items()}


def determine_filename(source_path_text: str) -> str:
    source_name = Path(source_path_text).name
    if source_name.endswith(".tar.gz") or source_name.endswith(".gz"):
        return source_name
    return f"{source_name}.tar.gz"


def determine_file_format(filename: str) -> str:
    return "tar.gz" if filename.endswith(".tar.gz") else "gz"


def select_bundle_contents(bundle: BundleReader) -> str:
    entries = set(bundle.list_root_entries())
    selected = [name for name in ORDERED_BUNDLE_CONTENTS if name in entries]
    if selected:
        return ",".join(selected)
    return ",".join(sorted(entries))


def cluster_count(bundle: BundleReader) -> int | None:
    cluster_path = "analysis/clustering/gene_expression_graphclust/clusters.csv"
    if not bundle.exists(cluster_path):
        return None
    text = bundle.read_text(cluster_path)
    reader = csv.DictReader(io.StringIO(text))
    clusters = {as_clean_string(row.get("Cluster")) for row in reader if as_clean_string(row.get("Cluster"))}
    return len(clusters) or None


def has_pca(bundle: BundleReader) -> bool:
    return bundle.exists("analysis/pca/gene_expression_10_components/projection.csv")


def has_umap(bundle: BundleReader) -> bool:
    return bundle.exists("analysis/umap/gene_expression_2_components/projection.csv")


def normalize_platform(value: str) -> str:
    text = as_clean_string(value)
    if text.lower() == "xenium":
        return "10x Genomics Xenium"
    return text


def normalize_segmented_object_type(value: str) -> str:
    text = as_clean_string(value)
    lookup = {
        "whole cell": "Whole cell",
        "nucleus": "nucleus",
        "cytoplasm": "cytoplasm",
    }
    return lookup.get(text.lower(), text)


def required_fields_for_row(row: dict[str, str]) -> list[str]:
    required = list(REQUIRED_FIELDS_BASE)

    if parse_bool(row.get("RNA_MEASURED")):
        required.append("TRANSCRIPTOME_TYPE")
    if parse_bool(row.get("RNA_MEASURED")) and as_clean_string(row.get("TRANSCRIPTOME_TYPE")) == "Targeted":
        required.append("PANEL_NAME")
    if parse_bool(row.get("HAS_CELL_SEGMENTATION")):
        required.extend(
            [
                "CELL_SEGMENTATION_METHOD",
                "CELL_SEGMENTED_OBJECT_TYPE",
                "NUMBER_OF_SEGMENTED_CELLS",
            ]
        )
    if parse_bool(row.get("HAS_CLUSTERING")):
        required.extend(["CLUSTERING_METHOD", "NUMBER_OF_CLUSTERS"])
    if normalize_platform(row.get("PLATFORM")) == "10x Genomics Xenium":
        required.append("SLIDE_SERIAL_NUMBER")

    deduped: list[str] = []
    seen = set()
    for field in required:
        if field not in seen:
            deduped.append(field)
            seen.add(field)
    return deduped


def derive_panel_name(
    row: dict[str, str],
    panel_reference: PanelReference,
    gene_panel: dict | None,
) -> str:
    specimen = as_clean_string(row.get("WUSTL Specimen"))
    panel_key = panel_reference.specimen_to_panel_key.get(specimen)
    if panel_key and panel_key in panel_reference.definitions_by_key:
        return panel_reference.definitions_by_key[panel_key].display_name

    if gene_panel:
        panel_identity = (((gene_panel.get("payload") or {}).get("panel") or {}).get("identity") or {})
        if panel_identity.get("name"):
            return as_clean_string(panel_identity["name"])
        if panel_identity.get("design_id"):
            return as_clean_string(panel_identity["design_id"])

    return ""


def derive_panel_size_total_targets(gene_panel: dict | None) -> str:
    if not gene_panel:
        return ""
    payload = gene_panel.get("payload") or {}
    targets = payload.get("targets") or []
    gene_targets = [target for target in targets if ((target.get("type") or {}).get("descriptor") == "gene")]
    return str(len(gene_targets)) if gene_targets else ""


def derive_software_string(experiment: dict) -> str:
    analysis_sw = as_clean_string(experiment.get("analysis_sw_version"))
    instrument_sw = as_clean_string(experiment.get("instrument_sw_version"))
    parts = []
    if analysis_sw:
        parts.append(f"analysis_sw_version={analysis_sw}")
    if instrument_sw:
        parts.append(f"instrument_sw_version={instrument_sw}")
    return "; ".join(parts)


def derive_run_id(row: dict[str, str], experiment: dict) -> str:
    path_text = as_clean_string(row.get("WUSTL Path"))
    if path_text:
        parent_name = Path(path_text).parent.name
        if parent_name:
            return parent_name
    return as_clean_string(experiment.get("run_name"))


def derive_assay_chemistry_version(row: dict[str, str], experiment: dict, gene_panel: dict | None) -> str:
    chemistry = as_clean_string(experiment.get("chemistry_version"))
    if chemistry:
        return chemistry
    if gene_panel:
        chemistry = as_clean_string((((gene_panel.get("payload") or {}).get("chemistry") or {}).get("version")))
        if chemistry:
            return chemistry
    existing = as_clean_string(row.get("ASSAY_CHEMISTRY_VERSION"))
    if existing == "Xenium Prime":
        return "v2"
    return existing


def derive_number_of_segmented_cells(experiment: dict, metrics: dict[str, str]) -> str:
    if not is_blank(experiment.get("num_cells")):
        return as_clean_string(experiment.get("num_cells"))
    if not is_blank(metrics.get("num_cells_detected")):
        return as_clean_string(metrics.get("num_cells_detected"))
    return ""


def derive_metadata_fields(
    row: dict[str, str],
    bundle: BundleReader,
    panel_reference: PanelReference,
) -> tuple[dict[str, str], dict | None, list[str]]:
    warnings: list[str] = []
    experiment = read_json(bundle, "experiment.xenium")
    metrics = read_single_row_csv(bundle, "metrics_summary.csv") if bundle.exists("metrics_summary.csv") else {}
    gene_panel = read_json(bundle, "gene_panel.json") if bundle.exists("gene_panel.json") else None

    filename = determine_filename(as_clean_string(row.get("WUSTL Path")))
    derived = {
        "ASSAY_CHEMISTRY_VERSION": derive_assay_chemistry_version(row, experiment, gene_panel),
        "BUNDLE_CONTENTS": select_bundle_contents(bundle),
        "CLUSTERING_METHOD": "Graph-Based" if bundle.exists("analysis/clustering/gene_expression_graphclust/clusters.csv") else "",
        "DIMENSIONALITY_REDUCTION_METHOD": "PCA" if has_pca(bundle) else ("UMAP" if has_umap(bundle) else ""),
        "FILE_FORMAT": determine_file_format(filename),
        "FILENAME": filename,
        "HAS_DIMENSIONALITY_REDUCTION": "True" if (has_pca(bundle) or has_umap(bundle)) else "False",
        "NUMBER_OF_CLUSTERS": as_clean_string(cluster_count(bundle)),
        "NUMBER_OF_SEGMENTED_CELLS": derive_number_of_segmented_cells(experiment, metrics),
        "PANEL_NAME": derive_panel_name(row, panel_reference, gene_panel),
        "PANEL_SIZE_TOTAL_TARGETS": derive_panel_size_total_targets(gene_panel),
        "PLATFORM": "10x Genomics Xenium",
        "PORTAL_PREVIEW_FILE": bundle.find_root_html(),
        "RUN_ID": derive_run_id(row, experiment),
        "SLIDE_SERIAL_NUMBER": as_clean_string(experiment.get("slide_id")),
        "SOFTWARE_AND_VERSION": derive_software_string(experiment),
        "SPATIAL_ASSAY_TYPE": "In situ",
    }

    if not bundle.exists("gene_panel.json"):
        warnings.append("Missing gene_panel.json")
    if not derived["NUMBER_OF_CLUSTERS"] and parse_bool(row.get("HAS_CLUSTERING")):
        warnings.append("Clustering declared but graphclust cluster file was not found")

    return derived, gene_panel, warnings


def build_panel_outputs(
    output_dir: Path,
    panel_reference: PanelReference,
    rows: list[dict[str, str]],
    gene_panels_by_specimen: dict[str, dict],
    hgnc_version: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)

    panel_manifest: list[dict[str, str]] = []
    excluded_target_rows: list[dict[str, str]] = []
    written_panel_ids: set[str] = set()

    for row in rows:
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        gene_panel = gene_panels_by_specimen.get(specimen)
        if not gene_panel:
            continue

        panel_name = as_clean_string(row.get("PANEL_NAME"))
        if not panel_name:
            panel_name = derive_panel_name(row, panel_reference, gene_panel)

        panel_id = safe_panel_id(panel_name or specimen)
        if panel_id in written_panel_ids:
            continue
        written_panel_ids.add(panel_id)

        payload = gene_panel.get("payload") or {}
        targets = payload.get("targets") or []
        unique_valid: OrderedDict[tuple[str, str], dict[str, str]] = OrderedDict()
        invalid_targets: list[dict[str, str]] = []

        for target in targets:
            target_type = target.get("type") or {}
            if target_type.get("descriptor") != "gene":
                continue

            data = target_type.get("data") or {}
            gene_symbol = as_clean_string(data.get("name"))
            gene_id = as_clean_string(data.get("id"))
            key = (gene_symbol, gene_id)

            if VALID_HTAN_GENE_ID.fullmatch(gene_id):
                if key not in unique_valid:
                    unique_valid[key] = {
                        "HTAN_PANEL_ID": panel_id,
                        "GENE_SYMBOL": gene_symbol,
                        "HGNC_VERSION": hgnc_version,
                        "GENE_ID": gene_id,
                        "USER_GENE_NAME": "",
                    }
            else:
                invalid_targets.append(
                    {
                        "HTAN_PANEL_ID": panel_id,
                        "PANEL_NAME": panel_name,
                        "GENE_SYMBOL": gene_symbol,
                        "GENE_ID": gene_id,
                        "REASON": "GENE_ID is not ENSG* or numeric",
                    }
                )

        panel_rows = list(unique_valid.values())
        panel_output_path = panels_dir / f"{panel_id}.csv"
        csv_write_rows(
            panel_output_path,
            panel_rows,
            ["HTAN_PANEL_ID", "GENE_SYMBOL", "HGNC_VERSION", "GENE_ID", "USER_GENE_NAME"],
        )

        excluded_output_path = ""
        if invalid_targets:
            excluded_output = panels_dir / f"{panel_id}.excluded_targets.csv"
            csv_write_rows(
                excluded_output,
                invalid_targets,
                ["HTAN_PANEL_ID", "PANEL_NAME", "GENE_SYMBOL", "GENE_ID", "REASON"],
            )
            excluded_output_path = str(excluded_output)
            excluded_target_rows.extend(invalid_targets)

        ref_key = panel_reference.specimen_to_panel_key.get(specimen, "")
        ref_def = panel_reference.definitions_by_key.get(ref_key) if ref_key else None

        panel_manifest.append(
            {
                "HTAN_PANEL_ID": panel_id,
                "PANEL_NAME": panel_name,
                "SOURCE_SPECIMEN": specimen,
                "PANEL_TYPE": ref_def.panel_type if ref_def else "",
                "PREDESIGNED_PANEL": ref_def.predesigned_panel if ref_def else "",
                "EXPECTED_TARGETS_FROM_REFERENCE": ref_def.n_gene_total if ref_def else "",
                "GENE_TARGETS_FROM_GENE_PANEL_JSON": str(
                    sum(1 for target in targets if ((target.get("type") or {}).get("descriptor") == "gene"))
                ),
                "VALID_HTAN_PANEL_ROWS": str(len(panel_rows)),
                "EXCLUDED_NON_STANDARD_TARGETS": str(len(invalid_targets)),
                "PANEL_OUTPUT_FILE": str(panel_output_path),
                "EXCLUDED_TARGETS_FILE": excluded_output_path,
            }
        )

    return panel_manifest, excluded_target_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to the JSON configuration file.")
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    config_dir = config_path.parent
    config = json.loads(config_path.read_text())

    metadata_csv = resolve_path(config["metadata_csv"], config_dir)
    panel_reference_csv = resolve_path(config["panel_reference_csv"], config_dir)
    output_dir = resolve_path(config.get("output_dir", "outputs"), config_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_rows, fieldnames = csv_read_rows(metadata_csv)
    panel_reference = parse_panel_reference(panel_reference_csv)

    include_specimens = set(config.get("include_specimens") or [])
    path_overrides = {key: as_clean_string(value) for key, value in (config.get("path_overrides_by_specimen") or {}).items()}
    hgnc_version = as_clean_string(config.get("hgnc_version", ""))

    output_filename = config.get("metadata_output_filename") or f"{metadata_csv.stem}.filled.csv"
    metadata_output_path = output_dir / output_filename
    fill_log_path = output_dir / (Path(output_filename).stem + ".fill_log.csv")
    unresolved_path = output_dir / (Path(output_filename).stem + ".unresolved_required_fields.csv")
    panel_manifest_path = output_dir / (Path(output_filename).stem + ".panel_manifest.csv")

    fill_log_rows: list[dict[str, str]] = []
    unresolved_rows: list[dict[str, str]] = []
    processed_rows: list[dict[str, str]] = []
    gene_panels_by_specimen: dict[str, dict] = {}
    specimen_rows = [row for row in metadata_rows if not include_specimens or as_clean_string(row.get("WUSTL Specimen")) in include_specimens]

    log(f"Loaded metadata CSV: {metadata_csv}")
    log(f"Loaded panel reference CSV: {panel_reference_csv}")
    log(f"Output directory: {output_dir}")
    log(f"Specimens selected for processing: {len(specimen_rows)}")

    for index, row in enumerate(metadata_rows, start=1):
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        if include_specimens and specimen not in include_specimens:
            continue

        working_row = dict(row)

        row_path = as_clean_string(row.get("WUSTL Path"))
        source_path_text = path_overrides.get(specimen, row_path)
        source_path = resolve_path(source_path_text, config_dir)

        log(f"[{index}/{len(metadata_rows)}] Processing specimen: {specimen or '<blank>'}")
        log(f"  Source path: {source_path}")

        if not source_path.exists():
            log("  Source path not found; leaving row unresolved")
            unresolved_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "WUSTL Path",
                    "VALUE": source_path_text,
                    "DETAIL": "Source Xenium path does not exist",
                }
            )
            processed_rows.append(working_row)
            continue

        bundle = BundleReader(source_path)
        try:
            derived, gene_panel, warnings = derive_metadata_fields(working_row, bundle, panel_reference)
        finally:
            bundle.close()

        if gene_panel:
            gene_panels_by_specimen[specimen] = gene_panel

        for warning in warnings:
            fill_log_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "__warning__",
                    "OLD_VALUE": "",
                    "NEW_VALUE": "",
                    "ACTION": "warning",
                    "SOURCE": warning,
                }
            )

        for field, value in derived.items():
            if is_blank(value):
                continue

            old_value = as_clean_string(working_row.get(field, ""))
            if field in AUTHORITATIVE_FIELDS or is_blank(old_value):
                if old_value != value:
                    action = "filled" if is_blank(old_value) else "overrode"
                    fill_log_rows.append(
                        {
                            "WUSTL Specimen": specimen,
                            "FIELD": field,
                            "OLD_VALUE": old_value,
                            "NEW_VALUE": value,
                            "ACTION": action,
                            "SOURCE": "xenium_bundle",
                        }
                    )
                working_row[field] = value

        normalized_platform = normalize_platform(working_row.get("PLATFORM", ""))
        if normalized_platform and normalized_platform != working_row.get("PLATFORM", ""):
            fill_log_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "PLATFORM",
                    "OLD_VALUE": as_clean_string(working_row.get("PLATFORM")),
                    "NEW_VALUE": normalized_platform,
                    "ACTION": "normalized",
                    "SOURCE": "htan_enum",
                }
            )
            working_row["PLATFORM"] = normalized_platform

        normalized_object_type = normalize_segmented_object_type(working_row.get("CELL_SEGMENTED_OBJECT_TYPE", ""))
        if normalized_object_type and normalized_object_type != working_row.get("CELL_SEGMENTED_OBJECT_TYPE", ""):
            fill_log_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "CELL_SEGMENTED_OBJECT_TYPE",
                    "OLD_VALUE": as_clean_string(working_row.get("CELL_SEGMENTED_OBJECT_TYPE")),
                    "NEW_VALUE": normalized_object_type,
                    "ACTION": "normalized",
                    "SOURCE": "htan_enum",
                }
            )
            working_row["CELL_SEGMENTED_OBJECT_TYPE"] = normalized_object_type

        for field in required_fields_for_row(working_row):
            if is_blank(working_row.get(field)):
                unresolved_rows.append(
                    {
                        "WUSTL Specimen": specimen,
                        "FIELD": field,
                        "VALUE": "",
                        "DETAIL": "Required field still blank after automated filling",
                    }
                )

        log(
            "  Completed"
            f" | panel={as_clean_string(working_row.get('PANEL_NAME')) or '<blank>'}"
            f" | clusters={as_clean_string(working_row.get('NUMBER_OF_CLUSTERS')) or '<blank>'}"
            f" | segmented_cells={as_clean_string(working_row.get('NUMBER_OF_SEGMENTED_CELLS')) or '<blank>'}"
        )

        processed_rows.append(working_row)

    rows_to_write = []
    processed_by_specimen = {row.get("WUSTL Specimen", ""): row for row in processed_rows}
    for row in metadata_rows:
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        if specimen in processed_by_specimen:
            rows_to_write.append(processed_by_specimen[specimen])
        else:
            rows_to_write.append(row)

    panel_manifest_rows, excluded_target_rows = build_panel_outputs(
        output_dir=output_dir,
        panel_reference=panel_reference,
        rows=processed_rows,
        gene_panels_by_specimen=gene_panels_by_specimen,
        hgnc_version=hgnc_version,
    )

    csv_write_rows(metadata_output_path, rows_to_write, fieldnames)
    log(f"Wrote metadata output: {metadata_output_path}")
    csv_write_rows(
        fill_log_path,
        fill_log_rows,
        ["WUSTL Specimen", "FIELD", "OLD_VALUE", "NEW_VALUE", "ACTION", "SOURCE"],
    )
    log(f"Wrote fill log: {fill_log_path}")
    csv_write_rows(
        unresolved_path,
        unresolved_rows,
        ["WUSTL Specimen", "FIELD", "VALUE", "DETAIL"],
    )
    log(f"Wrote unresolved-fields report: {unresolved_path}")
    csv_write_rows(
        panel_manifest_path,
        panel_manifest_rows,
        [
            "HTAN_PANEL_ID",
            "PANEL_NAME",
            "SOURCE_SPECIMEN",
            "PANEL_TYPE",
            "PREDESIGNED_PANEL",
            "EXPECTED_TARGETS_FROM_REFERENCE",
            "GENE_TARGETS_FROM_GENE_PANEL_JSON",
            "VALID_HTAN_PANEL_ROWS",
            "EXCLUDED_NON_STANDARD_TARGETS",
            "PANEL_OUTPUT_FILE",
            "EXCLUDED_TARGETS_FILE",
        ],
    )
    log(f"Wrote panel manifest: {panel_manifest_path}")
    log(f"Generated panel CSVs: {len(panel_manifest_rows)}")

    summary = {
        "metadata_output": str(metadata_output_path),
        "fill_log_output": str(fill_log_path),
        "unresolved_output": str(unresolved_path),
        "panel_manifest_output": str(panel_manifest_path),
        "processed_specimen_count": len(processed_rows),
        "generated_panel_count": len(panel_manifest_rows),
        "excluded_nonstandard_panel_targets": len(excluded_target_rows),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

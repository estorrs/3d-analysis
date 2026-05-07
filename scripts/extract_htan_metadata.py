#!/usr/bin/env python3
"""Fill HTAN Spatial Omics Level 3 metadata from processed output bundles.

Usage:
    python extract_htan_metadata.py config.json

The config can contain one job or multiple jobs. Each job represents a modality
such as Xenium or Visium and produces:
  - a filled metadata CSV
  - a fill log CSV
  - an unresolved-required-fields CSV
  - optional panel CSVs and a panel manifest
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import sys
import tarfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import h5py


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

BOOL_TRUE = {"true", "t", "1", "yes", "y"}
BOOL_FALSE = {"false", "f", "0", "no", "n"}
VALID_HTAN_GENE_ID = re.compile(r"(ENSG\d+|\d+)$")

XENIUM_BUNDLE_CONTENTS = [
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

VISIUM_BUNDLE_CONTENTS = [
    "analysis",
    "binned_outputs",
    "cloupe.cloupe",
    "feature_slice.h5",
    "filtered_feature_bc_matrix",
    "filtered_feature_bc_matrix.h5",
    "metrics_summary.csv",
    "molecule_info.h5",
    "probe_set.csv",
    "raw_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
    "spatial",
    "web_summary.html",
]

OME_TIFF_EXTENSIONS = (".ome.tif", ".ome.tiff", ".ome.tf2", ".ome.tf8", ".ome.btf")

CODEX_METADATA_FIELDS = [
    "WUSTL Participant",
    "WUSTL Specimen",
    "WUSTL Path",
    "HTAN_DATA_FILE_ID",
    "HTAN_PARENT_ID",
    "FILE_FORMAT",
    "FILENAME",
    "WORKING_DISTANCE",
    "IMAGING_ASSAY_TYPE",
    "PYRAMID",
    "PHYSICAL_SIZE_X",
    "PHYSICAL_SIZE_Y",
    "PHYSICAL_SIZE_Z",
    "SIZE_C",
    "SIZE_T",
    "SIZE_X",
    "SIZE_Y",
    "SIZE_Z",
    "CHANNEL_METADATA_ID",
    "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES",
    "DE_IDENTIFICATION_METHOD_TYPE",
    "DE_IDENTIFICATION_METHOD_DESCRIPTION",
    "DE_IDENTIFICATION_SOFTWARE",
    "LICENSE",
    "IMAGE_MODALITY",
    "IMAGING_EQUIPMENT_MANUFACTURER",
    "IMAGING_EQUIPMENT_MODEL",
    "IMAGING_SOFTWARE",
    "CITATION_OR_DOI",
    "IMAGING_PROTOCOL",
    "STAINING_METHOD",
    "OBJECTIVE",
    "NOMINAL_MAGNIFICATION",
    "IMMERSION",
    "LENS_NUMERICAL_APERTURE",
    "PASSED_QC",
    "QC_COMMENT",
    "SPECIES",
    "HAS_SLIDE_LABEL",
    "SLIDE_LABEL_REDACTED",
    "DE_IDENTIFIED",
]

CODEX_REQUIRED_FIELDS = [
    "HTAN_DATA_FILE_ID",
    "HTAN_PARENT_ID",
    "FILE_FORMAT",
    "FILENAME",
    "IMAGING_ASSAY_TYPE",
    "PHYSICAL_SIZE_X",
    "PHYSICAL_SIZE_Y",
    "PHYSICAL_SIZE_Z",
    "SIZE_C",
    "SIZE_T",
    "SIZE_X",
    "SIZE_Y",
    "SIZE_Z",
    "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES",
    "DE_IDENTIFICATION_METHOD_TYPE",
    "LICENSE",
    "IMAGE_MODALITY",
    "IMAGING_EQUIPMENT_MANUFACTURER",
    "CITATION_OR_DOI",
    "STAINING_METHOD",
    "PASSED_QC",
    "SPECIES",
    "HAS_SLIDE_LABEL",
    "DE_IDENTIFIED",
]

CODEX_CHANNEL_FIELDS = [
    "CHANNEL_ID",
    "CHANNEL_NAME",
    "CYCLE_NUMBER",
    "SUB_CYCLE_NUMBER",
    "TARGET_NAME",
    "ANTIBODY_NAME",
    "RRID_IDENTIFIER",
    "FLUOROPHORE",
    "CLONE",
    "LOT",
    "CATALOG_NUMBER",
    "EXCITATION_WAVELENGTH",
    "EMISSION_WAVELENGTH",
    "EXCITATION_BANDWIDTH",
    "EMISSION_BANDWIDTH",
    "METAL_ISOTOPE_ELEMENT_ABBREVIATION",
    "METAL_ISOTOPE_ELEMENT_MASS",
    "OLIGO_BARCODE_UPPER_STRAND",
    "OLIGO_BARCODE_LOWER_STRAND",
    "DILUTION",
    "CONCENTRATION",
]

HE_METADATA_FIELDS = [
    "WUSTL Participant",
    "WUSTL Specimen",
    "WUSTL Path",
    "HTAN_DATA_FILE_ID",
    "HTAN_PARENT_ID",
    "CITATION_OR_DOI",
    "DE_IDENTIFICATION_METHOD_DESCRIPTION",
    "DE_IDENTIFICATION_METHOD_TYPE",
    "DE_IDENTIFICATION_SOFTWARE",
    "DE_IDENTIFIED",
    "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES",
    "HAS_SLIDE_LABEL",
    "IMAGE_MODALITY",
    "IMAGING_EQUIPMENT_MANUFACTURER",
    "IMAGING_EQUIPMENT_MODEL",
    "IMAGING_PROTOCOL",
    "IMAGING_SOFTWARE",
    "IMMERSION",
    "LENS_NUMERICAL_APERTURE",
    "LICENSE",
    "NOMINAL_MAGNIFICATION",
    "OBJECTIVE",
    "PASSED_QC",
    "QC_COMMENT",
    "SLIDE_LABEL_REDACTED",
    "SPECIES",
    "STAINING_METHOD",
    "ANNOTATION_TYPE",
    "FILENAME",
    "FILE_FORMAT",
    "HAS_ANNOTATIONS",
]

HE_REQUIRED_FIELDS = [
    "HTAN_DATA_FILE_ID",
    "HTAN_PARENT_ID",
    "CITATION_OR_DOI",
    "DE_IDENTIFICATION_METHOD_TYPE",
    "DE_IDENTIFIED",
    "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES",
    "HAS_SLIDE_LABEL",
    "IMAGE_MODALITY",
    "IMAGING_EQUIPMENT_MANUFACTURER",
    "LICENSE",
    "PASSED_QC",
    "SPECIES",
    "STAINING_METHOD",
    "FILENAME",
    "FILE_FORMAT",
    "HAS_ANNOTATIONS",
]


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


def resolve_path(path_text: str, base_dir: Path) -> Path:
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def resolve_command_path(path_text: str, base_dir: Path) -> str:
    text = as_clean_string(path_text)
    if not text:
        return ""
    if "/" in text or text.startswith(".") or text.startswith("~"):
        return str(resolve_path(text, base_dir))
    return text


def maybe_number_to_string(value: object) -> str:
    if value is None:
        return ""
    text = as_clean_string(value)
    if text == "None":
        return ""
    return text


def is_ome_tiff(filename: str) -> bool:
    return filename.lower().endswith(OME_TIFF_EXTENSIONS)


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def csv_read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = []
        for row in reader:
            rows.append({key: as_clean_string(value) for key, value in row.items()})
    return rows, fieldnames


def csv_write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


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


def normalize_platform(value: str) -> str:
    text = as_clean_string(value)
    lookup = {
        "xenium": "10x Genomics Xenium",
        "visium": "10x Genomics Visium",
        "visium hd": "10x Genomics Visium HD",
    }
    return lookup.get(text.lower(), text)


def normalize_segmented_object_type(value: str) -> str:
    text = as_clean_string(value)
    lookup = {
        "whole cell": "Whole cell",
        "nucleus": "nucleus",
        "cytoplasm": "cytoplasm",
    }
    return lookup.get(text.lower(), text)


def determine_file_format(filename: str) -> str:
    return "tar.gz" if filename.endswith(".tar.gz") else "gz"


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

    platform = normalize_platform(row.get("PLATFORM"))
    if platform in {"10x Genomics Xenium", "10x Genomics Visium", "10x Genomics Visium HD"}:
        required.append("SLIDE_SERIAL_NUMBER")
    if platform in {"10x Genomics Visium", "10x Genomics Visium HD"}:
        required.extend(["CAPTURE_AREA", "CYTASSIST_USED", "GENOMIC_REFERENCE"])
    if as_clean_string(row.get("SPATIAL_ASSAY_TYPE")) == "capture-based":
        required.extend(["SEQUENCING_INSTRUMENT", "SEQUENCING_CONFIGURATION", "SEQUENCING_DEPTH"])

    deduped = []
    seen = set()
    for field in required:
        if field not in seen:
            seen.add(field)
            deduped.append(field)
    return deduped


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


@dataclass
class CosMxTrackingRecord:
    specimen: str
    output_path: str
    rds_path: str
    fov_spec: str
    panel_hint: str
    sample_run_name: str


@dataclass
class CosMxTrackingReference:
    by_specimen: dict[str, CosMxTrackingRecord]
    by_output_path: dict[str, list[CosMxTrackingRecord]]


@dataclass
class ExtractionResult:
    derived: dict[str, str]
    authoritative_fields: set[str]
    panel_export: dict | None
    warnings: list[str]


def empty_panel_reference() -> PanelReference:
    return PanelReference(definitions_by_key={}, specimen_to_panel_key={})


def empty_cosmx_tracking_reference() -> CosMxTrackingReference:
    return CosMxTrackingReference(by_specimen={}, by_output_path={})


def parse_panel_reference(path: Path | None) -> PanelReference:
    if path is None or not path.exists():
        return empty_panel_reference()

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

    return PanelReference(definitions_by_key=definitions, specimen_to_panel_key=specimen_to_panel_key)


def parse_cosmx_tracking_reference(path: Path | None) -> CosMxTrackingReference:
    if path is None or not path.exists():
        return empty_cosmx_tracking_reference()

    rows, _ = csv_read_rows(path)
    by_specimen: dict[str, CosMxTrackingRecord] = {}
    by_output_path: dict[str, list[CosMxTrackingRecord]] = {}

    for row in rows:
        specimen = as_clean_string(row.get("Section_ID_clean"))
        output_path = as_clean_string(row.get("Output file path"))
        record = CosMxTrackingRecord(
            specimen=specimen,
            output_path=output_path,
            rds_path=as_clean_string(row.get("AtoMx RDS (Do not use for analysis)")),
            fov_spec=as_clean_string(row.get("FOVs")),
            panel_hint=as_clean_string(row.get("Panel")),
            sample_run_name=as_clean_string(row.get("sample_run_name")) or as_clean_string(row.get("folder_name")),
        )
        if specimen:
            by_specimen[specimen] = record
        if output_path:
            by_output_path.setdefault(output_path, []).append(record)

    return CosMxTrackingReference(by_specimen=by_specimen, by_output_path=by_output_path)


def build_codex_rows_from_survey(
    survey_csv: Path,
    include_specimens: set[str],
    experiment_filter: str,
    level_filter: str,
) -> list[dict[str, str]]:
    survey_rows, _ = csv_read_rows(survey_csv)
    rows: list[dict[str, str]] = []

    for survey_row in survey_rows:
        experiment = as_clean_string(survey_row.get("experiment")).lower()
        level = as_clean_string(survey_row.get("level")).lower()
        if experiment != experiment_filter.lower():
            continue
        if level_filter and level != level_filter.lower():
            continue

        specimen = as_clean_string(survey_row.get("Specimen name"))
        if include_specimens and specimen not in include_specimens:
            continue

        row = {field: "" for field in CODEX_METADATA_FIELDS}
        row["WUSTL Participant"] = as_clean_string(survey_row.get("participant name"))
        row["WUSTL Specimen"] = specimen
        row["WUSTL Path"] = as_clean_string(survey_row.get("path"))
        rows.append(row)

    return rows


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
            if parts:
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

        html_names = []
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


def determine_archive_stem(source_path_text: str) -> str:
    path = Path(source_path_text)
    name = path.name
    if name.endswith(".tar.gz"):
        return name[:-7]
    if name.endswith(".gz"):
        return name[:-3]
    if name == "outs":
        return path.parent.name
    return name


def determine_filename(source_path_text: str) -> str:
    name = Path(source_path_text).name
    if name.endswith(".tar.gz") or name.endswith(".gz"):
        return name
    return f"{determine_archive_stem(source_path_text)}.tar.gz"


def count_clusters_from_csv(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        clusters = {as_clean_string(row.get("Cluster")) for row in reader if as_clean_string(row.get("Cluster"))}
    return len(clusters) or None


def choose_first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def derive_xenium_panel_name(row: dict[str, str], panel_reference: PanelReference, gene_panel: dict | None) -> str:
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


def derive_xenium_panel_export(
    specimen: str,
    panel_name: str,
    panel_reference: PanelReference,
    gene_panel: dict | None,
) -> dict | None:
    if not gene_panel:
        return None

    targets = ((gene_panel.get("payload") or {}).get("targets") or [])
    unique_valid: OrderedDict[tuple[str, str], dict[str, str]] = OrderedDict()
    excluded_rows = []

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
                unique_valid[key] = {"GENE_SYMBOL": gene_symbol, "GENE_ID": gene_id}
        else:
            excluded_rows.append(
                {
                    "PANEL_NAME": panel_name,
                    "GENE_SYMBOL": gene_symbol,
                    "GENE_ID": gene_id,
                    "REASON": "GENE_ID is not ENSG* or numeric",
                }
            )

    ref_def = None
    panel_key = panel_reference.specimen_to_panel_key.get(specimen)
    if panel_key:
        ref_def = panel_reference.definitions_by_key.get(panel_key)

    return {
        "panel_name": panel_name,
        "panel_id": safe_panel_id(panel_name or specimen),
        "source_specimen": specimen,
        "panel_type": ref_def.panel_type if ref_def else "",
        "predesigned_panel": ref_def.predesigned_panel if ref_def else "",
        "expected_targets_from_reference": ref_def.n_gene_total if ref_def else "",
        "gene_targets_from_source": str(len(unique_valid)),
        "rows": list(unique_valid.values()),
        "excluded_rows": excluded_rows,
    }


def extract_xenium_row(
    row: dict[str, str],
    source_path: Path,
    panel_reference: PanelReference,
) -> ExtractionResult:
    warnings: list[str] = []
    bundle = BundleReader(source_path)
    try:
        experiment = read_json(bundle, "experiment.xenium")
        metrics = read_single_row_csv(bundle, "metrics_summary.csv") if bundle.exists("metrics_summary.csv") else {}
        gene_panel = read_json(bundle, "gene_panel.json") if bundle.exists("gene_panel.json") else None
        cluster_path = "analysis/clustering/gene_expression_graphclust/clusters.csv"
        pca_path = "analysis/pca/gene_expression_10_components/projection.csv"
        umap_path = "analysis/umap/gene_expression_2_components/projection.csv"

        chemistry = as_clean_string(experiment.get("chemistry_version"))
        if not chemistry and gene_panel:
            chemistry = as_clean_string((((gene_panel.get("payload") or {}).get("chemistry") or {}).get("version")))
        if not chemistry and as_clean_string(row.get("ASSAY_CHEMISTRY_VERSION")) == "Xenium Prime":
            chemistry = "v2"

        panel_name = derive_xenium_panel_name(row, panel_reference, gene_panel)
        gene_targets = ""
        if gene_panel:
            targets = ((gene_panel.get("payload") or {}).get("targets") or [])
            gene_targets = str(sum(1 for target in targets if ((target.get("type") or {}).get("descriptor") == "gene")))

        derived = {
            "ASSAY_CHEMISTRY_VERSION": chemistry,
            "BUNDLE_CONTENTS": ",".join([name for name in XENIUM_BUNDLE_CONTENTS if bundle.exists(name)]),
            "CLUSTERING_METHOD": "Graph-Based" if bundle.exists(cluster_path) else "",
            "DIMENSIONALITY_REDUCTION_METHOD": "PCA" if bundle.exists(pca_path) else ("UMAP" if bundle.exists(umap_path) else ""),
            "FILE_FORMAT": determine_file_format(determine_filename(str(source_path))),
            "FILENAME": determine_filename(str(source_path)),
            "HAS_DIMENSIONALITY_REDUCTION": "True" if (bundle.exists(pca_path) or bundle.exists(umap_path)) else "False",
            "NUMBER_OF_CLUSTERS": as_clean_string(
                count_clusters_from_csv(source_path / cluster_path) if source_path.is_dir() else ""
            ),
            "NUMBER_OF_SEGMENTED_CELLS": as_clean_string(
                experiment.get("num_cells") if not is_blank(experiment.get("num_cells")) else metrics.get("num_cells_detected")
            ),
            "PANEL_NAME": panel_name if as_clean_string(row.get("TRANSCRIPTOME_TYPE")) == "Targeted" else "",
            "PANEL_SIZE_TOTAL_TARGETS": gene_targets,
            "PLATFORM": "10x Genomics Xenium",
            "PORTAL_PREVIEW_FILE": bundle.find_root_html(),
            "RUN_ID": Path(as_clean_string(row.get("WUSTL Path"))).parent.name or as_clean_string(experiment.get("run_name")),
            "SLIDE_SERIAL_NUMBER": as_clean_string(experiment.get("slide_id")),
            "SOFTWARE_AND_VERSION": "; ".join(
                [
                    part
                    for part in [
                        f"analysis_sw_version={as_clean_string(experiment.get('analysis_sw_version'))}" if as_clean_string(experiment.get("analysis_sw_version")) else "",
                        f"instrument_sw_version={as_clean_string(experiment.get('instrument_sw_version'))}" if as_clean_string(experiment.get("instrument_sw_version")) else "",
                    ]
                    if part
                ]
            ),
            "SPATIAL_ASSAY_TYPE": "In situ",
        }

        if not bundle.exists("gene_panel.json"):
            warnings.append("Missing gene_panel.json")
        if not derived["NUMBER_OF_CLUSTERS"] and parse_bool(row.get("HAS_CLUSTERING")):
            warnings.append("Clustering declared but graphclust cluster file was not found")

        panel_export = derive_xenium_panel_export(
            specimen=as_clean_string(row.get("WUSTL Specimen")),
            panel_name=panel_name,
            panel_reference=panel_reference,
            gene_panel=gene_panel,
        )

        authoritative_fields = {
            "ASSAY_CHEMISTRY_VERSION",
            "BUNDLE_CONTENTS",
            "CLUSTERING_METHOD",
            "DIMENSIONALITY_REDUCTION_METHOD",
            "FILE_FORMAT",
            "FILENAME",
            "HAS_DIMENSIONALITY_REDUCTION",
            "NUMBER_OF_CLUSTERS",
            "NUMBER_OF_SEGMENTED_CELLS",
            "PLATFORM",
            "PORTAL_PREVIEW_FILE",
            "RUN_ID",
            "SLIDE_SERIAL_NUMBER",
            "SOFTWARE_AND_VERSION",
            "SPATIAL_ASSAY_TYPE",
        }
        if as_clean_string(row.get("TRANSCRIPTOME_TYPE")) == "Targeted":
            authoritative_fields.add("PANEL_NAME")

        return ExtractionResult(
            derived=derived,
            authoritative_fields=authoritative_fields,
            panel_export=panel_export,
            warnings=warnings,
        )
    finally:
        bundle.close()


def resolve_visium_roots(source_path: Path) -> tuple[Path, Path]:
    if source_path.name == "outs":
        return source_path.parent, source_path
    if (source_path / "outs").exists():
        return source_path, source_path / "outs"
    return source_path, source_path


def read_metrics_csv(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader, None)
    return {key: as_clean_string(value) for key, value in (row or {}).items()}


def read_probe_set(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    headers: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("#"):
                if "=" in line:
                    key, value = line[1:].split("=", 1)
                    headers[key.strip()] = value.strip()
                continue
            fieldnames = [field.strip() for field in line.split(",")]
            reader = csv.DictReader(handle, fieldnames=fieldnames)
            rows.extend({key: as_clean_string(value) for key, value in row.items()} for row in reader)
            break
    return headers, rows


def read_molecule_info_json(path: Path) -> dict:
    with h5py.File(path, "r") as handle:
        raw = handle["metrics_json"][()]
        if hasattr(raw, "decode"):
            return json.loads(raw.decode())
        return json.loads(str(raw))


def read_feature_genomes(path: Path) -> list[str]:
    with h5py.File(path, "r") as handle:
        genomes = handle["features"]["genome"][()]
    cleaned = []
    for genome in genomes:
        value = genome.decode() if isinstance(genome, bytes) else str(genome)
        value = value.strip()
        if value:
            cleaned.append(value)
    return sorted(set(cleaned))


def derive_slide_serial_and_capture_area(value: str) -> tuple[str, str]:
    text = as_clean_string(value)
    if not text or "-" not in text:
        return text, ""
    slide_serial, capture_area = text.rsplit("-", 1)
    return slide_serial, capture_area


def derive_visium_genomic_reference(
    run_root: Path,
    outs_root: Path,
    metrics_json: dict,
    probe_set_headers: dict[str, str],
) -> str:
    invocation_candidates = [run_root / "_invocation", outs_root / "_invocation", run_root.parent / "_invocation"]
    for candidate in invocation_candidates:
        if candidate.exists():
            text = candidate.read_text()
            match = re.search(r'reference_path\s*=\s*"([^"]+)"', text)
            if match:
                basename = Path(match.group(1)).name
                if basename.startswith("refdata-gex-"):
                    return basename.replace("refdata-gex-", "", 1)
                return basename

    if probe_set_headers.get("reference_genome") and probe_set_headers.get("reference_version"):
        return f"{probe_set_headers['reference_genome']}-{probe_set_headers['reference_version']}"

    target_panel_path = as_clean_string(metrics_json.get("target_panel_path"))
    match = re.search(r"(GRCh\d+-[0-9A-Za-z.-]+|GRCm\d+-[0-9A-Za-z.-]+)", target_panel_path)
    if match:
        return match.group(1)

    feature_h5 = choose_first_existing(
        [
            outs_root / "filtered_feature_bc_matrix.h5",
            outs_root / "binned_outputs" / "square_008um" / "filtered_feature_bc_matrix.h5",
            outs_root / "binned_outputs" / "square_016um" / "filtered_feature_bc_matrix.h5",
        ]
    )
    if feature_h5 is not None:
        genomes = read_feature_genomes(feature_h5)
        if genomes:
            return ",".join(genomes)

    return ""


def derive_visium_panel_export(specimen: str, probe_set_path: Path) -> dict | None:
    if not probe_set_path.exists():
        return None

    headers, rows = read_probe_set(probe_set_path)
    panel_name = headers.get("panel_name") or probe_set_path.stem
    unique_valid: OrderedDict[tuple[str, str], dict[str, str]] = OrderedDict()
    excluded_rows = []

    for row in rows:
        if row.get("included", "TRUE").upper() == "FALSE":
            continue
        gene_id = as_clean_string(row.get("gene_id"))
        deprecated_match = re.fullmatch(r"DEPRECATED_(ENSG\d+)", gene_id)
        if deprecated_match:
            gene_id = deprecated_match.group(1)
        probe_id = as_clean_string(row.get("probe_id"))
        gene_symbol = ""
        parts = probe_id.split("|")
        if len(parts) >= 2:
            gene_symbol = parts[1]
        key = (gene_symbol, gene_id)
        if VALID_HTAN_GENE_ID.fullmatch(gene_id):
            if key not in unique_valid:
                unique_valid[key] = {"GENE_SYMBOL": gene_symbol, "GENE_ID": gene_id}
        else:
            excluded_rows.append(
                {
                    "PANEL_NAME": panel_name,
                    "GENE_SYMBOL": gene_symbol,
                    "GENE_ID": gene_id,
                    "REASON": "GENE_ID is not ENSG* or numeric",
                }
            )

    return {
        "panel_name": panel_name,
        "panel_id": safe_panel_id(panel_name),
        "source_specimen": specimen,
        "panel_type": headers.get("panel_type", ""),
        "predesigned_panel": "",
        "expected_targets_from_reference": str(len(unique_valid)),
        "gene_targets_from_source": str(len(unique_valid)),
        "rows": list(unique_valid.values()),
        "excluded_rows": excluded_rows,
    }


def extract_visium_row(row: dict[str, str], source_path: Path) -> ExtractionResult:
    warnings: list[str] = []
    run_root, outs_root = resolve_visium_roots(source_path)

    metrics_path = outs_root / "metrics_summary.csv"
    if not metrics_path.exists():
        return ExtractionResult(derived={}, authoritative_fields=set(), panel_export=None, warnings=["Missing metrics_summary.csv"])

    metrics = read_metrics_csv(metrics_path)
    molecule_info_path = outs_root / "molecule_info.h5"
    metrics_json = read_molecule_info_json(molecule_info_path) if molecule_info_path.exists() else {}

    probe_set_path = outs_root / "probe_set.csv"
    probe_set_headers, probe_set_rows = ({}, [])
    if probe_set_path.exists():
        probe_set_headers, probe_set_rows = read_probe_set(probe_set_path)

    is_hd = (outs_root / "binned_outputs").exists() or "visium hd" in as_clean_string(metrics_json.get("chemistry_description")).lower()
    platform = "10x Genomics Visium HD" if is_hd else "10x Genomics Visium"
    analysis_root = choose_first_existing(
        [
            outs_root / "binned_outputs" / "square_008um" / "analysis",
            outs_root / "binned_outputs" / "square_016um" / "analysis",
            outs_root / "analysis",
        ]
    )

    graphclust_path = choose_first_existing(
        [
            (analysis_root / "clustering" / "gene_expression_graphclust" / "clusters.csv") if analysis_root else Path("__missing__"),
            (analysis_root / "clustering" / "graphclust" / "clusters.csv") if analysis_root else Path("__missing__"),
        ]
    )
    pca_path = choose_first_existing(
        [
            (analysis_root / "pca" / "gene_expression_10_components" / "projection.csv") if analysis_root else Path("__missing__"),
            (analysis_root / "pca" / "10_components" / "projection.csv") if analysis_root else Path("__missing__"),
        ]
    )
    umap_path = choose_first_existing(
        [
            (analysis_root / "umap" / "gene_expression_2_components" / "projection.csv") if analysis_root else Path("__missing__"),
            (analysis_root / "umap" / "2_components" / "projection.csv") if analysis_root else Path("__missing__"),
        ]
    )
    tsne_path = choose_first_existing(
        [
            (analysis_root / "tsne" / "gene_expression_2_components" / "projection.csv") if analysis_root else Path("__missing__"),
            (analysis_root / "tsne" / "2_components" / "projection.csv") if analysis_root else Path("__missing__"),
        ]
    )

    slide_serial_capture_area = as_clean_string(metrics_json.get("slide_serial_capture_area"))
    if not slide_serial_capture_area:
        invocation_candidates = [run_root / "_invocation", outs_root / "_invocation", run_root.parent / "_invocation"]
        for candidate in invocation_candidates:
            if candidate.exists():
                match = re.search(r'slide_serial_capture_area\s*=\s*"([^"]+)"', candidate.read_text())
                if match:
                    slide_serial_capture_area = match.group(1)
                    break

    slide_serial_number, capture_area = derive_slide_serial_and_capture_area(slide_serial_capture_area)
    cytassist_used = (outs_root / "spatial" / "cytassist_image.tiff").exists()
    if not cytassist_used:
        invocation_candidates = [run_root / "_invocation", outs_root / "_invocation", run_root.parent / "_invocation"]
        for candidate in invocation_candidates:
            if candidate.exists():
                match = re.search(r"cytassist_image_paths\s*=\s*\[(.*?)\]", candidate.read_text(), flags=re.DOTALL)
                if match and '"' in match.group(1):
                    cytassist_used = True
                    break

    panel_size_total_targets = as_clean_string(row.get("PANEL_SIZE_TOTAL_TARGETS"))
    if is_blank(panel_size_total_targets):
        if not is_blank(metrics_json.get("target_panel_gene_count")):
            panel_size_total_targets = as_clean_string(metrics_json.get("target_panel_gene_count"))
        elif probe_set_rows:
            panel_size_total_targets = str(len({as_clean_string(probe_row.get("gene_id")) for probe_row in probe_set_rows if as_clean_string(probe_row.get("gene_id"))}))

    derived = {
        "BUNDLE_CONTENTS": ",".join([name for name in VISIUM_BUNDLE_CONTENTS if (outs_root / name).exists()]),
        "CAPTURE_AREA": capture_area,
        "CLUSTERING_METHOD": "Graph-Based" if graphclust_path is not None and graphclust_path.exists() else "",
        "CYTASSIST_USED": "True" if cytassist_used else "False",
        "DIMENSIONALITY_REDUCTION_METHOD": (
            "PCA"
            if pca_path is not None and pca_path.exists()
            else ("UMAP" if umap_path is not None and umap_path.exists() else ("t-SNE" if tsne_path is not None and tsne_path.exists() else ""))
        ),
        "FILE_FORMAT": determine_file_format(determine_filename(str(source_path))),
        "FILENAME": determine_filename(str(source_path)),
        "GENOMIC_REFERENCE": derive_visium_genomic_reference(run_root, outs_root, metrics_json, probe_set_headers),
        "HAS_DIMENSIONALITY_REDUCTION": "True" if any(path is not None and path.exists() for path in [pca_path, umap_path, tsne_path]) else "False",
        "NUMBER_OF_CLUSTERS": as_clean_string(count_clusters_from_csv(graphclust_path) if graphclust_path is not None else ""),
        "PANEL_SIZE_TOTAL_TARGETS": panel_size_total_targets,
        "PLATFORM": platform,
        "PORTAL_PREVIEW_FILE": "web_summary.html" if (outs_root / "web_summary.html").exists() else "",
        "QC_SPATIAL_UNIT": "8um bin" if is_hd else "spot",
        "RUN_ID": run_root.name,
        "SLIDE_SERIAL_NUMBER": slide_serial_number,
        "SOFTWARE_AND_VERSION": as_clean_string(metrics_json.get("cellranger_version")),
    }

    if not derived["NUMBER_OF_CLUSTERS"] and parse_bool(row.get("HAS_CLUSTERING")):
        warnings.append("Clustering declared but graphclust cluster file was not found")
    if not metrics_json:
        warnings.append("Missing or unreadable molecule_info.h5 metrics_json")

    panel_export = derive_visium_panel_export(as_clean_string(row.get("WUSTL Specimen")), probe_set_path)

    return ExtractionResult(
        derived=derived,
        authoritative_fields={
            "BUNDLE_CONTENTS",
            "CAPTURE_AREA",
            "CLUSTERING_METHOD",
            "CYTASSIST_USED",
            "DIMENSIONALITY_REDUCTION_METHOD",
            "FILE_FORMAT",
            "FILENAME",
            "GENOMIC_REFERENCE",
            "HAS_DIMENSIONALITY_REDUCTION",
            "NUMBER_OF_CLUSTERS",
            "PLATFORM",
            "PORTAL_PREVIEW_FILE",
            "QC_SPATIAL_UNIT",
            "RUN_ID",
            "SLIDE_SERIAL_NUMBER",
            "SOFTWARE_AND_VERSION",
        },
        panel_export=panel_export,
        warnings=warnings,
    )


def find_cosmx_tracking_record(row: dict[str, str], tracking_reference: CosMxTrackingReference) -> CosMxTrackingRecord | None:
    specimen = as_clean_string(row.get("WUSTL Specimen"))
    if specimen and specimen in tracking_reference.by_specimen:
        return tracking_reference.by_specimen[specimen]

    output_path = as_clean_string(row.get("WUSTL Path"))
    matches = tracking_reference.by_output_path.get(output_path, [])
    if len(matches) == 1:
        return matches[0]
    return None


def resolve_cosmx_rds_path(source_path: Path, tracking_record: CosMxTrackingRecord | None) -> tuple[Path | None, list[str]]:
    warnings: list[str] = []
    if source_path.is_file():
        if source_path.suffix.lower() == ".rds":
            return source_path, warnings
        warnings.append("Source path points to a file, but it is not an .RDS file")
        return None, warnings

    if not source_path.is_dir():
        warnings.append("Source path is neither a directory nor an .RDS file")
        return None, warnings

    tracked_name = Path(tracking_record.rds_path).name if tracking_record and tracking_record.rds_path else ""
    candidates = sorted(path for path in source_path.iterdir() if path.is_file() and path.suffix.lower() == ".rds")

    if tracked_name:
        tracked_matches = [path for path in candidates if path.name == tracked_name]
        if len(tracked_matches) == 1:
            return tracked_matches[0], warnings

    if len(candidates) == 1:
        return candidates[0], warnings

    if tracking_record and tracking_record.rds_path:
        tracked_path = Path(tracking_record.rds_path)
        if tracked_path.exists():
            warnings.append("Using tracking CSV .RDS path because the metadata source directory did not resolve cleanly")
            return tracked_path, warnings
        tracked_in_dir = source_path / tracked_path.name
        if tracked_in_dir.exists():
            return tracked_in_dir, warnings

    if len(candidates) > 1:
        warnings.append("Multiple top-level .RDS files were found; unable to choose one automatically")
    else:
        warnings.append("No top-level .RDS file was found")
    return None, warnings


def cosmx_panel_size_from_hint(panel_hint: str) -> str:
    lookup = {
        "1k": "1000",
        "6k": "6175",
    }
    return lookup.get(as_clean_string(panel_hint).lower(), "")


def clean_cosmx_panel_name(panel_name: str, panel_hint: str) -> str:
    cleaned = re.sub(r"\s+", " ", as_clean_string(panel_name)).strip(" ;,")
    if cleaned:
        return cleaned
    hint = as_clean_string(panel_hint)
    return f"CosMx {hint}" if hint else ""


def run_cosmx_r_extractor(rscript_path: str, helper_script: Path, rds_path: Path, fov_spec: str) -> dict:
    completed = subprocess.run(
        [rscript_path, str(helper_script), str(rds_path), fov_spec],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        error_text = as_clean_string(completed.stderr) or as_clean_string(completed.stdout) or f"Rscript exited with code {completed.returncode}"
        raise RuntimeError(error_text)

    stdout = as_clean_string(completed.stdout)
    if not stdout:
        raise RuntimeError("Rscript returned no JSON output")
    try:
        return json.loads(stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Could not parse Rscript JSON output: {exc}") from exc


def extract_cosmx_row(
    row: dict[str, str],
    source_path: Path,
    tracking_reference: CosMxTrackingReference,
    rscript_path: str,
    helper_script: Path,
    cache: dict[tuple[str, str], dict],
) -> ExtractionResult:
    warnings: list[str] = []
    tracking_record = find_cosmx_tracking_record(row, tracking_reference)
    if tracking_record is None:
        warnings.append("No matching CosMx tracking row was found; cannot determine specimen-specific FOV subset")

    rds_path, path_warnings = resolve_cosmx_rds_path(source_path, tracking_record)
    warnings.extend(path_warnings)
    if rds_path is None:
        return ExtractionResult(derived={}, authoritative_fields=set(), panel_export=None, warnings=warnings)

    run_id = (
        as_clean_string(tracking_record.sample_run_name) if tracking_record else ""
    ) or (source_path.parent.name if source_path.is_file() else source_path.name)
    derived = {
        "BUNDLE_CONTENTS": rds_path.name,
        "FILE_FORMAT": "rds",
        "FILENAME": rds_path.name,
        "PLATFORM": "CosMx",
        "RUN_ID": run_id,
    }
    authoritative_fields = {"BUNDLE_CONTENTS", "FILE_FORMAT", "FILENAME", "PLATFORM", "RUN_ID"}

    if tracking_record is None:
        return ExtractionResult(derived=derived, authoritative_fields=authoritative_fields, panel_export=None, warnings=warnings)

    cache_key = (str(rds_path.resolve()), tracking_record.fov_spec)
    if cache_key not in cache:
        cache[cache_key] = run_cosmx_r_extractor(rscript_path, helper_script, rds_path, tracking_record.fov_spec)
    cosmx_data = cache[cache_key]

    for warning in cosmx_data.get("warnings", []):
        if as_clean_string(warning):
            warnings.append(as_clean_string(warning))

    panel_size = as_clean_string(cosmx_data.get("panel_size_total_targets")) or cosmx_panel_size_from_hint(tracking_record.panel_hint)
    derived.update(
        {
            "ASSAY_CHEMISTRY_VERSION": as_clean_string(cosmx_data.get("assay_chemistry_version")),
            "CELL_SEGMENTATION_METHOD": as_clean_string(cosmx_data.get("cell_segmentation_method")) or "CosMx cell segmentation",
            "CELL_SEGMENTED_OBJECT_TYPE": as_clean_string(cosmx_data.get("cell_segmented_object_type")) or "Whole cell",
            "DIMENSIONALITY_REDUCTION_METHOD": as_clean_string(cosmx_data.get("dimensionality_reduction_method")),
            "HAS_CELL_SEGMENTATION": as_clean_string(cosmx_data.get("has_cell_segmentation")) or "True",
            "HAS_CLUSTERING": "False",
            "HAS_DIMENSIONALITY_REDUCTION": as_clean_string(cosmx_data.get("has_dimensionality_reduction")),
            "NUMBER_OF_SEGMENTED_CELLS": as_clean_string(cosmx_data.get("number_of_segmented_cells")),
            "PANEL_NAME": clean_cosmx_panel_name(as_clean_string(cosmx_data.get("panel_name")), tracking_record.panel_hint),
            "PANEL_SIZE_TOTAL_TARGETS": panel_size,
            "PROTEIN_MEASURED": as_clean_string(cosmx_data.get("protein_measured")) or "False",
            "QC_FEATURE_NUMBER": as_clean_string(cosmx_data.get("qc_feature_number")),
            "QC_MEAN_READS_PER_FEATURE": as_clean_string(cosmx_data.get("qc_mean_reads_per_feature")),
            "QC_SPATIAL_UNIT": as_clean_string(cosmx_data.get("qc_spatial_unit")) or "cell",
            "QC_TOTAL_GENES_DETECTED": as_clean_string(cosmx_data.get("qc_total_genes_detected")) or panel_size,
            "QC_TOTAL_NUMBER_OF_READS": as_clean_string(cosmx_data.get("qc_total_number_of_reads")),
            "RNA_MEASURED": as_clean_string(cosmx_data.get("rna_measured")) or "True",
            "SAME_SECTION_IMAGING_CHANNELS": as_clean_string(cosmx_data.get("same_section_imaging_channels")),
            "SLIDE_SERIAL_NUMBER": as_clean_string(cosmx_data.get("slide_serial_number")),
            "SOFTWARE_AND_VERSION": as_clean_string(cosmx_data.get("software_and_version")),
            "SPATIAL_ASSAY_TYPE": as_clean_string(cosmx_data.get("spatial_assay_type")) or "In situ",
            "TRANSCRIPTOME_TYPE": as_clean_string(cosmx_data.get("transcriptome_type")) or "Targeted",
        }
    )
    authoritative_fields.update(
        {
            "ASSAY_CHEMISTRY_VERSION",
            "CELL_SEGMENTATION_METHOD",
            "CELL_SEGMENTED_OBJECT_TYPE",
            "DIMENSIONALITY_REDUCTION_METHOD",
            "HAS_CELL_SEGMENTATION",
            "HAS_CLUSTERING",
            "HAS_DIMENSIONALITY_REDUCTION",
            "NUMBER_OF_SEGMENTED_CELLS",
            "PANEL_NAME",
            "PANEL_SIZE_TOTAL_TARGETS",
            "PROTEIN_MEASURED",
            "QC_FEATURE_NUMBER",
            "QC_MEAN_READS_PER_FEATURE",
            "QC_SPATIAL_UNIT",
            "QC_TOTAL_GENES_DETECTED",
            "QC_TOTAL_NUMBER_OF_READS",
            "RNA_MEASURED",
            "SAME_SECTION_IMAGING_CHANNELS",
            "SLIDE_SERIAL_NUMBER",
            "SOFTWARE_AND_VERSION",
            "SPATIAL_ASSAY_TYPE",
            "TRANSCRIPTOME_TYPE",
        }
    )

    return ExtractionResult(derived=derived, authoritative_fields=authoritative_fields, panel_export=None, warnings=warnings)


def clean_codex_target_name(channel_name: str) -> str:
    text = as_clean_string(channel_name)
    text = re.sub(r"\s*[-/]?\(D\)\s*$", "", text)
    return text.strip()


def extract_codex_row(source_path: Path) -> tuple[dict[str, str], list[dict[str, str]], list[str]]:
    from ome_types import from_tiff

    warnings: list[str] = []
    ome = from_tiff(source_path)
    if not ome.images:
        raise ValueError("OME metadata does not contain any images")

    image = ome.images[0]
    pixels = image.pixels
    if pixels is None:
        raise ValueError("OME metadata does not contain pixel metadata")

    derived = {
        "FILE_FORMAT": "ome-tiff" if is_ome_tiff(source_path.name) else source_path.suffix.lstrip("."),
        "FILENAME": source_path.name,
        "IMAGING_ASSAY_TYPE": "CODEX",
        "PHYSICAL_SIZE_X": maybe_number_to_string(pixels.physical_size_x),
        "PHYSICAL_SIZE_Y": maybe_number_to_string(pixels.physical_size_y),
        "PHYSICAL_SIZE_Z": maybe_number_to_string(pixels.physical_size_z) or "0",
        "SIZE_C": maybe_number_to_string(pixels.size_c),
        "SIZE_T": maybe_number_to_string(pixels.size_t),
        "SIZE_X": maybe_number_to_string(pixels.size_x),
        "SIZE_Y": maybe_number_to_string(pixels.size_y),
        "SIZE_Z": maybe_number_to_string(pixels.size_z),
        "CHANNEL_METADATA_ID": "",
        "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES": "Pathological",
        "DE_IDENTIFICATION_METHOD_TYPE": "Automatic",
        "LICENSE": "CC BY 4.0",
        "IMAGE_MODALITY": "SM",
        "IMAGING_EQUIPMENT_MANUFACTURER": "Akoya",
        "CITATION_OR_DOI": "https://doi.org/10.1158/2159-8290.CD-26-0012",
        "STAINING_METHOD": "CODEX",
        "OBJECTIVE": "",
        "NOMINAL_MAGNIFICATION": "",
        "PASSED_QC": "TRUE",
        "QC_COMMENT": "",
        "SPECIES": "9606 (Homo sapiens)",
        "HAS_SLIDE_LABEL": "FALSE",
        "DE_IDENTIFIED": "TRUE",
    }

    channel_rows: list[dict[str, str]] = []
    for index, channel in enumerate(pixels.channels):
        channel_id = as_clean_string(channel.id) or f"Channel:{index}"
        channel_name = as_clean_string(channel.name) or channel_id
        channel_rows.append(
            {
                "CHANNEL_ID": channel_id,
                "CHANNEL_NAME": channel_name,
                "CYCLE_NUMBER": "",
                "SUB_CYCLE_NUMBER": "",
                "TARGET_NAME": clean_codex_target_name(channel_name),
                "ANTIBODY_NAME": "",
                "RRID_IDENTIFIER": "",
                "FLUOROPHORE": as_clean_string(channel.fluor),
                "CLONE": "",
                "LOT": "",
                "CATALOG_NUMBER": "",
                "EXCITATION_WAVELENGTH": maybe_number_to_string(channel.excitation_wavelength),
                "EMISSION_WAVELENGTH": maybe_number_to_string(channel.emission_wavelength),
                "EXCITATION_BANDWIDTH": "",
                "EMISSION_BANDWIDTH": "",
                "METAL_ISOTOPE_ELEMENT_ABBREVIATION": "",
                "METAL_ISOTOPE_ELEMENT_MASS": "",
                "OLIGO_BARCODE_UPPER_STRAND": "",
                "OLIGO_BARCODE_LOWER_STRAND": "",
                "DILUTION": "",
                "CONCENTRATION": "",
            }
        )

    return derived, channel_rows, warnings


def required_fields_for_codex_row(row: dict[str, str]) -> list[str]:
    required = list(CODEX_REQUIRED_FIELDS)
    if parse_bool(row.get("HAS_SLIDE_LABEL")):
        required.append("SLIDE_LABEL_REDACTED")
    return required


def extract_he_row(source_path: Path) -> tuple[dict[str, str], list[str]]:
    from ome_types import from_tiff

    ome = from_tiff(source_path)
    warnings: list[str] = []

    derived = {
        "CITATION_OR_DOI": "https://doi.org/10.1158/2159-8290.CD-26-0012",
        "DE_IDENTIFICATION_METHOD_TYPE": "Automatic",
        "DE_IDENTIFIED": "TRUE",
        "EXPERIMENTAL_STRATEGY_AND_DATA_SUBTYPES": "Pathological",
        "HAS_SLIDE_LABEL": "FALSE",
        "IMAGE_MODALITY": "SM",
        "IMAGING_EQUIPMENT_MANUFACTURER": "Akoya",
        "IMAGING_SOFTWARE": as_clean_string(getattr(ome, "creator", "")),
        "LICENSE": "CC BY 4.0",
        "NOMINAL_MAGNIFICATION": "",
        "OBJECTIVE": "",
        "PASSED_QC": "TRUE",
        "QC_COMMENT": "",
        "SPECIES": "9606 (Homo sapiens)",
        "STAINING_METHOD": "H&E",
        "FILENAME": source_path.name,
        "FILE_FORMAT": "ome-tiff" if is_ome_tiff(source_path.name) else source_path.suffix.lstrip("."),
        "HAS_ANNOTATIONS": "FALSE",
    }

    return derived, warnings


def required_fields_for_he_row(row: dict[str, str]) -> list[str]:
    required = list(HE_REQUIRED_FIELDS)
    if parse_bool(row.get("HAS_SLIDE_LABEL")):
        required.append("SLIDE_LABEL_REDACTED")
    if parse_bool(row.get("HAS_ANNOTATIONS")):
        required.append("ANNOTATION_TYPE")

    deduped = []
    seen = set()
    for field in required:
        if field not in seen:
            seen.add(field)
            deduped.append(field)
    return deduped


def process_codex_job(job: dict, config_dir: Path) -> dict[str, object]:
    survey_csv = resolve_path(job["survey_csv"], config_dir)
    output_dir = resolve_path(job.get("output_dir", "outputs_codex"), config_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    include_specimens = set(job.get("include_specimens") or [])
    path_overrides = {key: as_clean_string(value) for key, value in (job.get("path_overrides_by_specimen") or {}).items()}
    experiment_filter = as_clean_string(job.get("experiment_filter", "codex")) or "codex"
    level_filter = as_clean_string(job.get("level_filter", "level2")) or "level2"

    metadata_rows = build_codex_rows_from_survey(survey_csv, include_specimens, experiment_filter, level_filter)
    fieldnames = list(CODEX_METADATA_FIELDS)

    output_filename = job.get("metadata_output_filename") or "HTAN2 3D Prostate Breast_ Multiplex Microscopy Level 2 - CODEX.generated.csv"
    metadata_output_path = output_dir / output_filename
    fill_log_path = output_dir / (Path(output_filename).stem + ".fill_log.csv")
    unresolved_path = output_dir / (Path(output_filename).stem + ".unresolved_required_fields.csv")
    channel_manifest_path = output_dir / (Path(output_filename).stem + ".channel_manifest.csv")
    channels_dir = output_dir / "channels"
    channels_dir.mkdir(parents=True, exist_ok=True)

    job_name = as_clean_string(job.get("name") or "codex")
    log(f"=== Job: {job_name} (codex) ===")
    log(f"Loaded survey CSV: {survey_csv}")
    log(f"Output directory: {output_dir}")
    log(f"Specimens selected for processing: {len(metadata_rows)}")

    fill_log_rows: list[dict[str, str]] = []
    unresolved_rows: list[dict[str, str]] = []
    channel_manifest_rows: list[dict[str, str]] = []

    for index, row in enumerate(metadata_rows, start=1):
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        source_path_text = path_overrides.get(specimen, as_clean_string(row.get("WUSTL Path")))
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
                    "DETAIL": "Source CODEX image path does not exist",
                }
            )
            continue

        try:
            derived, channel_rows, warnings = extract_codex_row(source_path)
        except Exception as exc:
            log(f"  Failed to read OME metadata: {exc}")
            unresolved_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "__file__",
                    "VALUE": source_path.name,
                    "DETAIL": f"Could not read OME metadata: {exc}",
                }
            )
            continue

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
            old_value = as_clean_string(row.get(field, ""))
            if old_value != value:
                fill_log_rows.append(
                    {
                        "WUSTL Specimen": specimen,
                        "FIELD": field,
                        "OLD_VALUE": old_value,
                        "NEW_VALUE": value,
                        "ACTION": "filled",
                        "SOURCE": "codex_ome",
                    }
                )
            row[field] = value

        channel_filename = f"{safe_panel_id(specimen)}.channels.csv"
        channel_output_path = channels_dir / channel_filename
        csv_write_rows(channel_output_path, channel_rows, CODEX_CHANNEL_FIELDS)
        channel_manifest_rows.append(
            {
                "WUSTL Participant": as_clean_string(row.get("WUSTL Participant")),
                "WUSTL Specimen": specimen,
                "SOURCE_IMAGE_FILE": source_path.name,
                "CHANNEL_COUNT": str(len(channel_rows)),
                "CHANNEL_OUTPUT_FILE": str(channel_output_path),
            }
        )

        for field in required_fields_for_codex_row(row):
            if is_blank(row.get(field)):
                unresolved_rows.append(
                    {
                        "WUSTL Specimen": specimen,
                        "FIELD": field,
                        "VALUE": "",
                        "DETAIL": "Required field still blank after automated template generation",
                    }
                )

        log(
            "  Completed"
            f" | assay={as_clean_string(row.get('IMAGING_ASSAY_TYPE')) or '<blank>'}"
            f" | size_c={as_clean_string(row.get('SIZE_C')) or '<blank>'}"
            f" | channels_csv={channel_filename}"
        )

    csv_write_rows(metadata_output_path, metadata_rows, fieldnames)
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
        channel_manifest_path,
        channel_manifest_rows,
        ["WUSTL Participant", "WUSTL Specimen", "SOURCE_IMAGE_FILE", "CHANNEL_COUNT", "CHANNEL_OUTPUT_FILE"],
    )
    log(f"Wrote channel manifest: {channel_manifest_path}")
    log(f"Generated channel CSVs: {len(channel_manifest_rows)}")

    return {
        "job_name": job_name,
        "modality": "codex",
        "metadata_output": str(metadata_output_path),
        "fill_log_output": str(fill_log_path),
        "unresolved_output": str(unresolved_path),
        "channel_manifest_output": str(channel_manifest_path),
        "processed_specimen_count": len(metadata_rows),
        "generated_channel_count": len(channel_manifest_rows),
    }


def process_he_job(job: dict, config_dir: Path) -> dict[str, object]:
    survey_csv = resolve_path(job["survey_csv"], config_dir)
    output_dir = resolve_path(job.get("output_dir", "outputs_he"), config_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    include_specimens = set(job.get("include_specimens") or [])
    path_overrides = {key: as_clean_string(value) for key, value in (job.get("path_overrides_by_specimen") or {}).items()}
    experiment_filter = as_clean_string(job.get("experiment_filter", "he")) or "he"
    level_filter = as_clean_string(job.get("level_filter", "level2")) or "level2"

    metadata_rows = build_codex_rows_from_survey(survey_csv, include_specimens, experiment_filter, level_filter)
    for row in metadata_rows:
        for field in list(row.keys()):
            if field not in HE_METADATA_FIELDS:
                row.pop(field, None)
        for field in HE_METADATA_FIELDS:
            row.setdefault(field, "")

    fieldnames = list(HE_METADATA_FIELDS)
    output_filename = job.get("metadata_output_filename") or "HTAN2 3D Prostate Breast_ Digital Pathology Level 2 - HE.generated.csv"
    metadata_output_path = output_dir / output_filename
    fill_log_path = output_dir / (Path(output_filename).stem + ".fill_log.csv")
    unresolved_path = output_dir / (Path(output_filename).stem + ".unresolved_required_fields.csv")

    job_name = as_clean_string(job.get("name") or "he")
    log(f"=== Job: {job_name} (he) ===")
    log(f"Loaded survey CSV: {survey_csv}")
    log(f"Output directory: {output_dir}")
    log(f"Specimens selected for processing: {len(metadata_rows)}")

    fill_log_rows: list[dict[str, str]] = []
    unresolved_rows: list[dict[str, str]] = []

    for index, row in enumerate(metadata_rows, start=1):
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        source_path_text = path_overrides.get(specimen, as_clean_string(row.get("WUSTL Path")))
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
                    "DETAIL": "Source H&E image path does not exist",
                }
            )
            continue

        try:
            derived, warnings = extract_he_row(source_path)
        except Exception as exc:
            log(f"  Failed to read OME metadata: {exc}")
            unresolved_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "__file__",
                    "VALUE": source_path.name,
                    "DETAIL": f"Could not read OME metadata: {exc}",
                }
            )
            continue

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
            old_value = as_clean_string(row.get(field, ""))
            if old_value != value:
                fill_log_rows.append(
                    {
                        "WUSTL Specimen": specimen,
                        "FIELD": field,
                        "OLD_VALUE": old_value,
                        "NEW_VALUE": value,
                        "ACTION": "filled",
                        "SOURCE": "he_ome",
                    }
                )
            row[field] = value

        for field in required_fields_for_he_row(row):
            if is_blank(row.get(field)):
                unresolved_rows.append(
                    {
                        "WUSTL Specimen": specimen,
                        "FIELD": field,
                        "VALUE": "",
                        "DETAIL": "Required field still blank after automated template generation",
                    }
                )

        log(
            "  Completed"
            f" | staining={as_clean_string(row.get('STAINING_METHOD')) or '<blank>'}"
            f" | file_format={as_clean_string(row.get('FILE_FORMAT')) or '<blank>'}"
            f" | software={as_clean_string(row.get('IMAGING_SOFTWARE')) or '<blank>'}"
        )

    csv_write_rows(metadata_output_path, metadata_rows, fieldnames)
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

    return {
        "job_name": job_name,
        "modality": "he",
        "metadata_output": str(metadata_output_path),
        "fill_log_output": str(fill_log_path),
        "unresolved_output": str(unresolved_path),
        "processed_specimen_count": len(metadata_rows),
    }


def build_panel_outputs(output_dir: Path, panel_exports: list[dict], hgnc_version: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    excluded_target_rows = []
    seen_panel_ids = set()

    for panel_export in panel_exports:
        if not panel_export:
            continue
        panel_id = panel_export["panel_id"]
        if panel_id in seen_panel_ids:
            continue
        seen_panel_ids.add(panel_id)

        panel_rows = []
        for row in panel_export["rows"]:
            panel_rows.append(
                {
                    "HTAN_PANEL_ID": panel_id,
                    "GENE_SYMBOL": row["GENE_SYMBOL"],
                    "HGNC_VERSION": hgnc_version,
                    "GENE_ID": row["GENE_ID"],
                    "USER_GENE_NAME": "",
                }
            )

        panel_output_path = panels_dir / f"{panel_id}.csv"
        csv_write_rows(
            panel_output_path,
            panel_rows,
            ["HTAN_PANEL_ID", "GENE_SYMBOL", "HGNC_VERSION", "GENE_ID", "USER_GENE_NAME"],
        )

        excluded_output_path = ""
        excluded_rows = []
        for row in panel_export["excluded_rows"]:
            excluded_rows.append(
                {
                    "HTAN_PANEL_ID": panel_id,
                    "PANEL_NAME": panel_export["panel_name"],
                    "GENE_SYMBOL": row["GENE_SYMBOL"],
                    "GENE_ID": row["GENE_ID"],
                    "REASON": row["REASON"],
                }
            )
        if excluded_rows:
            excluded_output = panels_dir / f"{panel_id}.excluded_targets.csv"
            csv_write_rows(
                excluded_output,
                excluded_rows,
                ["HTAN_PANEL_ID", "PANEL_NAME", "GENE_SYMBOL", "GENE_ID", "REASON"],
            )
            excluded_output_path = str(excluded_output)
            excluded_target_rows.extend(excluded_rows)

        manifest_rows.append(
            {
                "HTAN_PANEL_ID": panel_id,
                "PANEL_NAME": panel_export["panel_name"],
                "SOURCE_SPECIMEN": panel_export["source_specimen"],
                "PANEL_TYPE": panel_export["panel_type"],
                "PREDESIGNED_PANEL": panel_export["predesigned_panel"],
                "EXPECTED_TARGETS_FROM_REFERENCE": panel_export["expected_targets_from_reference"],
                "GENE_TARGETS_FROM_SOURCE": panel_export["gene_targets_from_source"],
                "VALID_HTAN_PANEL_ROWS": str(len(panel_rows)),
                "EXCLUDED_NON_STANDARD_TARGETS": str(len(excluded_rows)),
                "PANEL_OUTPUT_FILE": str(panel_output_path),
                "EXCLUDED_TARGETS_FILE": excluded_output_path,
            }
        )

    return manifest_rows, excluded_target_rows


def extract_row(row: dict[str, str], source_path: Path, modality: str, context: dict[str, object]) -> ExtractionResult:
    if modality == "xenium":
        return extract_xenium_row(row, source_path, context["panel_reference"])
    if modality == "visium":
        return extract_visium_row(row, source_path)
    if modality == "cosmx":
        return extract_cosmx_row(
            row,
            source_path,
            context["cosmx_tracking_reference"],
            context["cosmx_rscript_path"],
            context["cosmx_helper_script"],
            context["cosmx_cache"],
        )
    raise ValueError(f"Unsupported modality: {modality}")


def process_job(job: dict, config_dir: Path) -> dict[str, object]:
    modality = as_clean_string(job.get("modality")).lower()
    if not modality:
        raise ValueError("Each job must define a modality")

    if modality == "codex":
        return process_codex_job(job, config_dir)
    if modality == "he":
        return process_he_job(job, config_dir)

    metadata_csv = resolve_path(job["metadata_csv"], config_dir)
    output_dir = resolve_path(job.get("output_dir", f"outputs_{modality}"), config_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    panel_reference_csv = None
    if job.get("panel_reference_csv"):
        panel_reference_csv = resolve_path(job["panel_reference_csv"], config_dir)
    panel_reference = parse_panel_reference(panel_reference_csv)
    context: dict[str, object] = {"panel_reference": panel_reference}

    if modality == "cosmx":
        tracking_csv = resolve_path(job["tracking_csv"], config_dir)
        rscript_path = resolve_command_path(job.get("rscript_path", "Rscript"), config_dir)
        helper_script = Path(__file__).resolve().with_name("extract_cosmx_subset_metadata.R")
        if not helper_script.exists():
            raise FileNotFoundError(f"Missing helper script: {helper_script}")
        context.update(
            {
                "cosmx_tracking_reference": parse_cosmx_tracking_reference(tracking_csv),
                "cosmx_rscript_path": rscript_path,
                "cosmx_helper_script": helper_script,
                "cosmx_cache": {},
            }
        )

    metadata_rows, fieldnames = csv_read_rows(metadata_csv)
    include_specimens = set(job.get("include_specimens") or [])
    path_overrides = {key: as_clean_string(value) for key, value in (job.get("path_overrides_by_specimen") or {}).items()}
    hgnc_version = as_clean_string(job.get("hgnc_version", ""))

    output_filename = job.get("metadata_output_filename") or f"{metadata_csv.stem}.filled.csv"
    metadata_output_path = output_dir / output_filename
    fill_log_path = output_dir / (Path(output_filename).stem + ".fill_log.csv")
    unresolved_path = output_dir / (Path(output_filename).stem + ".unresolved_required_fields.csv")
    panel_manifest_path = output_dir / (Path(output_filename).stem + ".panel_manifest.csv")

    selected_indices = [index for index, row in enumerate(metadata_rows) if not include_specimens or as_clean_string(row.get("WUSTL Specimen")) in include_specimens]

    job_name = as_clean_string(job.get("name") or modality)
    log(f"=== Job: {job_name} ({modality}) ===")
    log(f"Loaded metadata CSV: {metadata_csv}")
    if panel_reference_csv is not None:
        log(f"Loaded panel reference CSV: {panel_reference_csv}")
    if modality == "cosmx":
        log(f"Loaded CosMx tracking CSV: {tracking_csv}")
        log(f"Using Rscript: {rscript_path}")
    log(f"Output directory: {output_dir}")
    log(f"Specimens selected for processing: {len(selected_indices)}")

    fill_log_rows: list[dict[str, str]] = []
    unresolved_rows: list[dict[str, str]] = []
    processed_by_index: dict[int, dict[str, str]] = {}
    panel_exports: list[dict] = []

    processed_counter = 0
    total_selected = len(selected_indices)
    for index, row in enumerate(metadata_rows):
        specimen = as_clean_string(row.get("WUSTL Specimen"))
        if include_specimens and specimen not in include_specimens:
            continue

        processed_counter += 1
        working_row = dict(row)
        source_path_text = path_overrides.get(specimen, as_clean_string(row.get("WUSTL Path")))
        source_path = resolve_path(source_path_text, config_dir)

        log(f"[{processed_counter}/{total_selected}] Processing specimen: {specimen or '<blank>'}")
        log(f"  Source path: {source_path}")

        if not source_path.exists():
            log("  Source path not found; leaving row unresolved")
            unresolved_rows.append(
                {
                    "WUSTL Specimen": specimen,
                    "FIELD": "WUSTL Path",
                    "VALUE": source_path_text,
                    "DETAIL": f"Source {modality} path does not exist",
                }
            )
            processed_by_index[index] = working_row
            continue

        result = extract_row(working_row, source_path, modality, context)
        if result.panel_export:
            panel_exports.append(result.panel_export)

        for warning in result.warnings:
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

        for field, value in result.derived.items():
            if is_blank(value):
                continue
            old_value = as_clean_string(working_row.get(field, ""))
            if field in result.authoritative_fields or is_blank(old_value):
                if old_value != value:
                    fill_log_rows.append(
                        {
                            "WUSTL Specimen": specimen,
                            "FIELD": field,
                            "OLD_VALUE": old_value,
                            "NEW_VALUE": value,
                            "ACTION": "filled" if is_blank(old_value) else "overrode",
                            "SOURCE": f"{modality}_bundle",
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
            if field == "PANEL_SYNAPSE_ID":
                continue
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
            f" | platform={as_clean_string(working_row.get('PLATFORM')) or '<blank>'}"
            f" | clusters={as_clean_string(working_row.get('NUMBER_OF_CLUSTERS')) or '<blank>'}"
            f" | segmented_cells={as_clean_string(working_row.get('NUMBER_OF_SEGMENTED_CELLS')) or '<blank>'}"
        )
        processed_by_index[index] = working_row

    rows_to_write = []
    for index, row in enumerate(metadata_rows):
        rows_to_write.append(processed_by_index.get(index, row))

    panel_manifest_rows, excluded_target_rows = build_panel_outputs(output_dir, panel_exports, hgnc_version)

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
            "GENE_TARGETS_FROM_SOURCE",
            "VALID_HTAN_PANEL_ROWS",
            "EXCLUDED_NON_STANDARD_TARGETS",
            "PANEL_OUTPUT_FILE",
            "EXCLUDED_TARGETS_FILE",
        ],
    )
    log(f"Wrote panel manifest: {panel_manifest_path}")
    log(f"Generated panel CSVs: {len(panel_manifest_rows)}")

    return {
        "job_name": job_name,
        "modality": modality,
        "metadata_output": str(metadata_output_path),
        "fill_log_output": str(fill_log_path),
        "unresolved_output": str(unresolved_path),
        "panel_manifest_output": str(panel_manifest_path),
        "processed_specimen_count": len(selected_indices),
        "generated_panel_count": len(panel_manifest_rows),
        "excluded_nonstandard_panel_targets": len(excluded_target_rows),
    }


def resolve_jobs(config: dict) -> list[dict]:
    if "jobs" in config:
        return list(config["jobs"])

    single_job = dict(config)
    if "modality" not in single_job:
        single_job["modality"] = "xenium"
    if "name" not in single_job:
        single_job["name"] = single_job["modality"]
    return [single_job]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to the JSON configuration file.")
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    config_dir = config_path.parent
    config = json.loads(config_path.read_text())

    summaries = []
    for job in resolve_jobs(config):
        summaries.append(process_job(job, config_dir))

    print(json.dumps(summaries if len(summaries) > 1 else summaries[0], indent=2))


if __name__ == "__main__":
    main()

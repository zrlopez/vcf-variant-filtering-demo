"""annotate_clinvar.py — ClinVar clinical significance annotation for VCF variants.

Enriches filtered VCF records with ClinVar clinical significance labels
by performing a lookup against a local ClinVar summary TSV (downloaded
from NCBI FTP; not included in this repository due to file size).

Annotation strategy:
  - Primary key: (CHROM, POS, REF, ALT) — exact positional match
  - Fallback key: rsID (when present in both VCF and ClinVar)
  - Unknown variants receive significance = 'not_found'

Output: VCF records with an additional 'clinvar_significance' key.

Usage::

    from scripts.annotate_clinvar import annotate_records, load_clinvar_index

    index = load_clinvar_index(Path("data/clinvar_summary.tsv"))
    annotated = annotate_records(records, index)

CLI usage::

    python scripts/annotate_clinvar.py \\
        --vcf     results/filtered.vcf \\
        --clinvar data/clinvar_summary.tsv \\
        --output  results/annotated.vcf
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

import click

from scripts.utils import parse_vcf_records, write_vcf_records

# ClinVar significance values we surface in reports (canonical form)
CLINSIG_CANONICAL: dict[str, str] = {
    "pathogenic": "Pathogenic",
    "likely pathogenic": "Likely_pathogenic",
    "uncertain significance": "VUS",
    "likely benign": "Likely_benign",
    "benign": "Benign",
    "conflicting interpretations of pathogenicity": "Conflicting",
    "not provided": "Not_provided",
}

# Positional index key type alias
_PosKey = tuple[str, int, str, str]  # (CHROM, POS, REF, ALT)
_RsKey = str  # rsID e.g. "rs12345"


# ---------------------------------------------------------------------------
# ClinVar index loader
# ---------------------------------------------------------------------------

def load_clinvar_index(
    clinvar_tsv: Path,
) -> dict[_PosKey | _RsKey, str]:
    """Load a ClinVar summary TSV into an in-memory lookup dict.

    Expects a tab-delimited file with at minimum these columns:
        Chromosome, Start, ReferenceAllele, AlternateAllele,
        ClinicalSignificance, RS# (db SNP)

    The index is keyed by both positional tuple AND rsID for maximum
    recall when either coordinate system is available.

    Args:
        clinvar_tsv: Path to the ClinVar summary file.

    Returns:
        Dict mapping (chrom, pos, ref, alt) → significance string
        and rsID → significance string.

    Raises:
        FileNotFoundError: If *clinvar_tsv* does not exist.
        ValueError: If the TSV is missing required column headers.
    """
    if not clinvar_tsv.exists():
        raise FileNotFoundError(f"ClinVar TSV not found: {clinvar_tsv}")

    required_cols = {
        "Chromosome", "Start", "ReferenceAllele",
        "AlternateAllele", "ClinicalSignificance",
    }
    index: dict[Any, str] = {}

    with clinvar_tsv.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("ClinVar TSV appears empty or has no header row.")
        missing = required_cols - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"ClinVar TSV missing required columns: {sorted(missing)}"
            )

        for row in reader:
            raw_sig = row["ClinicalSignificance"].strip().lower()
            sig = CLINSIG_CANONICAL.get(raw_sig, "Other")

            try:
                pos_key: _PosKey = (
                    row["Chromosome"].strip(),
                    int(row["Start"]),
                    row["ReferenceAllele"].strip(),
                    row["AlternateAllele"].strip(),
                )
                index[pos_key] = sig
            except (ValueError, KeyError):
                pass  # Skip malformed rows; do not fail entire load

            rs_raw = row.get("RS# (db SNP)", "").strip()
            if rs_raw and rs_raw not in ("-1", ".", ""):
                rs_key: _RsKey = f"rs{rs_raw}" if not rs_raw.startswith("rs") else rs_raw
                index[rs_key] = sig

    return index


# ---------------------------------------------------------------------------
# Annotation engine
# ---------------------------------------------------------------------------

def annotate_records(
    records: list[dict[str, Any]],
    clinvar_index: dict[Any, str],
) -> list[dict[str, Any]]:
    """Attach ClinVar clinical significance to each variant record.

    Attempts lookup in order:
      1. Positional key (CHROM, POS, REF, ALT) — highest specificity
      2. rsID key (if record['id'] is a valid rs identifier)

    Records with no match receive 'not_found'.

    Args:
        records:       Variant record dicts from parse_vcf_records() or
                       filter_variants().
        clinvar_index: Lookup dict from load_clinvar_index().

    Returns:
        A new list with 'clinvar_significance' added to every record.
        Input records are NOT mutated.
    """
    annotated: list[dict[str, Any]] = []
    for record in records:
        # Try positional lookup first
        pos_key: _PosKey = (
            record["chrom"],
            record["pos"],
            record["ref"],
            record["alt"],
        )
        sig = clinvar_index.get(pos_key)

        # Fall back to rsID lookup
        if sig is None:
            rs_id = record.get("id", "") or ""
            if rs_id.startswith("rs"):
                sig = clinvar_index.get(rs_id)

        enriched = {**record, "clinvar_significance": sig or "not_found"}
        annotated.append(enriched)

    return annotated


def get_annotation_summary(
    annotated_records: list[dict[str, Any]],
) -> dict[str, int]:
    """Return a frequency count of clinvar_significance values.

    Useful for generating per-run QC summaries and structured log output.
    """
    summary: dict[str, int] = {}
    for r in annotated_records:
        sig = r.get("clinvar_significance", "unknown")
        summary[sig] = summary.get(sig, 0) + 1
    return dict(sorted(summary.items()))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

@click.command()
@click.option(
    "--vcf",
    "vcf_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Filtered VCF input file.",
)
@click.option(
    "--clinvar",
    "clinvar_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="ClinVar summary TSV (downloaded from NCBI FTP).",
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Annotated output VCF file path.",
)
def main(vcf_path: Path, clinvar_path: Path, output_path: Path) -> None:
    """Annotate filtered VCF variants with ClinVar clinical significance."""
    try:
        records = parse_vcf_records(vcf_path)
        click.echo(f"Loaded {len(records):,} variants from {vcf_path}")

        clinvar_index = load_clinvar_index(clinvar_path)
        click.echo(f"Loaded ClinVar index: {len(clinvar_index):,} entries")

        annotated = annotate_records(records, clinvar_index)
        write_vcf_records(annotated, output_path)

        summary = get_annotation_summary(annotated)
        click.echo("\nAnnotation summary:")
        for sig, count in summary.items():
            click.echo(f"  {sig:<30} {count:>6,}")
        click.echo(f"\nWritten to: {output_path}")

    except (FileNotFoundError, ValueError) as exc:
        click.echo(f"ERROR: {exc}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

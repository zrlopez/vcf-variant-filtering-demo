"""filter_variants.py — VCF variant filtering pipeline.

Applies three independent filters to VCF variants:
  1. Minor Allele Frequency (MAF) — removes rare variants below threshold
  2. Quality score (QUAL) — removes low-confidence calls
  3. FILTER field — optionally restricts to PASS / '.' variants only

All filtering logic is encapsulated in pure functions for testability.
The CLI entry point is thin: it delegates to run_filter_pipeline() and
handles user-facing error reporting.

Usage::

    python scripts/filter_variants.py \\
        --input  data/sample.vcf \\
        --output results/filtered.vcf \\
        --maf 0.01 \\
        --min-qual 30

    # Skip FILTER-field check:
    python scripts/filter_variants.py --input ... --output ... --no-require-pass
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from scripts.utils import (
    parse_vcf_records,
    passes_maf_threshold,
    passes_qual_threshold,
    write_vcf_records,
)


# ---------------------------------------------------------------------------
# Core filtering logic — pure functions, no I/O, fully unit-testable
# ---------------------------------------------------------------------------

def filter_variants(
    records: list[dict[str, Any]],
    maf_threshold: float = 0.01,
    min_qual: float = 30.0,
    require_pass: bool = True,
) -> list[dict[str, Any]]:
    """Apply MAF, QUAL, and FILTER criteria to a list of variant records.

    Each filter is applied in order; a record is removed on the first
    failure, without evaluating remaining criteria (short-circuit).

    Args:
        records:        List of variant dicts from parse_vcf_records().
        maf_threshold:  Minimum MAF, inclusive.  Must be in (0, 0.5].
        min_qual:       Minimum QUAL score, inclusive.  Must be >= 0.
        require_pass:   If True, only variants with FILTER == 'PASS' or '.'
                        are retained.

    Returns:
        A new list containing only records that passed all active filters.

    Raises:
        ValueError: If maf_threshold is outside (0, 0.5] or min_qual < 0.
    """
    if not (0 < maf_threshold <= 0.5):
        raise ValueError(
            f"maf_threshold must be in (0, 0.5]; got {maf_threshold!r}"
        )
    if min_qual < 0:
        raise ValueError(f"min_qual must be >= 0; got {min_qual!r}")

    kept: list[dict[str, Any]] = []
    for record in records:
        # ── FILTER field gate ─────────────────────────────────
        if require_pass:
            filt = record.get("filter", "PASS")
            if filt not in ("PASS", ".", ""):
                continue
        # ── MAF gate ────────────────────────────────────────
        if not passes_maf_threshold(record, maf_threshold):
            continue
        # ── QUAL gate ───────────────────────────────────────
        if not passes_qual_threshold(record, min_qual):
            continue
        kept.append(record)

    return kept


def run_filter_pipeline(
    input_path: Path,
    output_path: Path,
    maf_threshold: float = 0.01,
    min_qual: float = 30.0,
    require_pass: bool = True,
) -> dict[str, Any]:
    """End-to-end filtering pipeline with a structured result dict.

    Orchestrates parse → filter → write, and returns provenance metadata
    suitable for logging or downstream report generation.

    Returns:
        dict with keys:
            total_input (int), total_output (int), removed (int),
            removal_pct (float), output_path (str)
    """
    records = parse_vcf_records(input_path)
    total_input = len(records)

    filtered = filter_variants(
        records,
        maf_threshold=maf_threshold,
        min_qual=min_qual,
        require_pass=require_pass,
    )
    total_output = len(filtered)
    write_vcf_records(filtered, output_path)

    removal_pct = (
        round((total_input - total_output) / total_input * 100, 2)
        if total_input > 0
        else 0.0
    )
    return {
        "total_input": total_input,
        "total_output": total_output,
        "removed": total_input - total_output,
        "removal_pct": removal_pct,
        "output_path": str(output_path),
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

@click.command()
@click.option(
    "--input", "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Input VCF file path.",
)
@click.option(
    "--output", "output_path",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output VCF file path (created if absent).",
)
@click.option(
    "--maf",
    default=0.01,
    show_default=True,
    type=click.FloatRange(min=1e-9, max=0.5),
    help="Minimum minor allele frequency (inclusive).",
)
@click.option(
    "--min-qual",
    default=30.0,
    show_default=True,
    type=click.FloatRange(min=0.0),
    help="Minimum QUAL score (inclusive).",
)
@click.option(
    "--require-pass/--no-require-pass",
    default=True,
    show_default=True,
    help="Retain only PASS/. variants (default: on).",
)
def main(
    input_path: Path,
    output_path: Path,
    maf: float,
    min_qual: float,
    require_pass: bool,
) -> None:
    """Filter VCF variants by MAF, QUAL, and FILTER field criteria."""
    try:
        result = run_filter_pipeline(
            input_path=input_path,
            output_path=output_path,
            maf_threshold=maf,
            min_qual=min_qual,
            require_pass=require_pass,
        )
    except (FileNotFoundError, ValueError) as exc:
        click.echo(f"ERROR: {exc}", err=True)
        sys.exit(1)

    click.echo(f"Input variants : {result['total_input']:>8,}")
    click.echo(f"Output variants: {result['total_output']:>8,}")
    click.echo(f"Removed        : {result['removed']:>8,}  ({result['removal_pct']}%)")
    click.echo(f"Written to     : {result['output_path']}")


if __name__ == "__main__":
    main()

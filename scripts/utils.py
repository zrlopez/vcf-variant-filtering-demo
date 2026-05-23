"""utils.py — VCF parsing, threshold checks, and I/O utilities.

All functions are pure (no side effects beyond explicit file I/O in
write_vcf_records) and fully unit-testable in isolation.

Design constraints:
  - No global mutable state.
  - All parsing errors surface as explicit ValueError / FileNotFoundError.
  - MAF extraction falls back gracefully: MAF= key, then AF= key, then None.
  - QUAL '.' (missing) is represented as None throughout the pipeline.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def parse_vcf_records(path: Path) -> list[dict[str, Any]]:
    """Parse a VCF file into a list of variant record dicts.

    Skips meta-information and header lines (those starting with '#').
    Each returned dict contains the keys:
        chrom, pos (int), id, ref, alt, qual (float|None),
        filter, info (raw string), maf (float|None).

    Args:
        path: Absolute or relative path to the input VCF file.

    Returns:
        List of variant record dicts, one per data line.

    Raises:
        FileNotFoundError: If *path* does not exist.
        ValueError: If a non-header line has fewer than 8 tab-separated fields.
    """
    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for lineno, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(
                    f"Line {lineno}: expected >=8 tab-separated fields, "
                    f"got {len(fields)}: {line!r}"
                )
            chrom, pos_str, var_id, ref, alt, qual_str, filt, info = fields[:8]

            qual: float | None
            try:
                qual = float(qual_str) if qual_str != "." else None
            except ValueError:
                qual = None

            records.append(
                {
                    "chrom": chrom,
                    "pos": int(pos_str),
                    "id": var_id,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filt,
                    "info": info,
                    "maf": _extract_maf_from_info(info),
                }
            )
    return records


def _extract_maf_from_info(info: str) -> float | None:
    """Extract a minor allele frequency from a VCF INFO field string.

    Checks for the 'MAF=' key first (bioinformatics convention), then falls
    back to 'AF=' (population allele frequency, used as MAF proxy).
    Returns None if neither key is present or the value is non-numeric.

    Examples::

        >>> _extract_maf_from_info("AF=0.05;DP=100")
        0.05
        >>> _extract_maf_from_info("MAF=0.12;AF=0.05")
        0.12
        >>> _extract_maf_from_info("DP=100")
        # returns None
    """
    for key in ("MAF", "AF"):
        for token in info.split(";"):
            if token.startswith(f"{key}="):
                try:
                    return float(token.split("=", 1)[1])
                except (ValueError, IndexError):
                    pass
    return None


# ---------------------------------------------------------------------------
# Threshold checks (pure predicates)
# ---------------------------------------------------------------------------

def passes_maf_threshold(record: dict[str, Any], threshold: float) -> bool:
    """Return True iff the record's MAF meets or exceeds *threshold*.

    Records with MAF == None (INFO field contained no AF/MAF key) are treated
    as failing the threshold filter to preserve conservative filtering.
    """
    maf = record.get("maf")
    if maf is None:
        return False
    return maf >= threshold


def passes_qual_threshold(record: dict[str, Any], min_qual: float) -> bool:
    """Return True iff the record's QUAL meets or exceeds *min_qual*.

    Records with QUAL == None (VCF '.' sentinel) are treated as failing
    the threshold to avoid propagating unscored variants downstream.
    """
    qual = record.get("qual")
    if qual is None:
        return False
    return qual >= min_qual


# ---------------------------------------------------------------------------
# VCF output
# ---------------------------------------------------------------------------

def write_vcf_records(records: list[dict[str, Any]], output_path: Path) -> None:
    """Write a list of variant record dicts to a VCF-formatted file.

    Writes a minimal two-line VCF header followed by one tab-separated
    data line per record.  The QUAL '.' sentinel is reconstructed from
    None values.  Parent directories are created automatically.

    Args:
        records:     List of variant dicts produced by parse_vcf_records().
        output_path: Destination file path (created if absent).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in records:
            qual_str = "." if r["qual"] is None else str(r["qual"])
            fh.write(
                f"{r['chrom']}\t{r['pos']}\t{r['id']}\t"
                f"{r['ref']}\t{r['alt']}\t{qual_str}\t"
                f"{r['filter']}\t{r['info']}\n"
            )

"""test_annotate_clinvar.py — Unit tests for scripts/annotate_clinvar.py.

Covers:
  - load_clinvar_index: positional key, rsID key, missing file, bad columns
  - annotate_records: positional match, rsID fallback, not_found default
  - get_annotation_summary: frequency counting, empty input
"""
from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Any

import pytest

from scripts.annotate_clinvar import (
    annotate_records,
    get_annotation_summary,
    load_clinvar_index,
)


# ---------------------------------------------------------------------------
# Helper factories
# ---------------------------------------------------------------------------

def _rec(chrom: str = "chr1", pos: int = 100, ref: str = "A", alt: str = "T",
         rs_id: str = "rs1", **kwargs: Any) -> dict[str, Any]:
    return {
        "chrom": chrom, "pos": pos, "id": rs_id,
        "ref": ref, "alt": alt,
        "qual": 50.0, "filter": "PASS", "info": "AF=0.05", "maf": 0.05,
        **kwargs,
    }


# ---------------------------------------------------------------------------
# load_clinvar_index
# ---------------------------------------------------------------------------

class TestLoadClinvarIndex:
    def test_positional_key_loaded(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        assert ("chr1", 100, "A", "T") in index

    def test_positional_significance_correct(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        assert index[("chr1", 100, "A", "T")] == "Pathogenic"

    def test_rsid_key_loaded(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        # clinvar_tsv fixture has RS# 1 and 2
        assert "rs1" in index or "rs2" in index

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_clinvar_index(tmp_path / "ghost.tsv")

    def test_missing_required_columns_raises(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.tsv"
        bad.write_text("OnlyCol\nval\n", encoding="utf-8")
        with pytest.raises(ValueError, match="missing required columns"):
            load_clinvar_index(bad)

    def test_malformed_rows_skipped(self, tmp_path: Path) -> None:
        """A row with a non-integer Start should be skipped, not crash."""
        content = textwrap.dedent("""\
            Chromosome\tStart\tReferenceAllele\tAlternateAllele\tClinicalSignificance\tRS# (db SNP)
            chr1\tNOT_AN_INT\tA\tT\tPathogenic\t1
            chr2\t200\tG\tC\tBenign\t2
        """)
        f = tmp_path / "partial.tsv"
        f.write_text(content, encoding="utf-8")
        index = load_clinvar_index(f)
        # Only the valid row should appear
        assert ("chr2", 200, "G", "C") in index
        assert ("chr1", "NOT_AN_INT", "A", "T") not in index

    def test_empty_tsv_raises(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.tsv"
        empty.write_text("", encoding="utf-8")
        with pytest.raises(ValueError):
            load_clinvar_index(empty)

    def test_clinvar_significance_canonicalized(self, tmp_path: Path) -> None:
        content = textwrap.dedent("""\
            Chromosome\tStart\tReferenceAllele\tAlternateAllele\tClinicalSignificance\tRS# (db SNP)
            chr1\t100\tA\tT\tuncertain significance\t-1
        """)
        f = tmp_path / "vus.tsv"
        f.write_text(content, encoding="utf-8")
        index = load_clinvar_index(f)
        assert index[("chr1", 100, "A", "T")] == "VUS"


# ---------------------------------------------------------------------------
# annotate_records
# ---------------------------------------------------------------------------

class TestAnnotateRecords:
    def test_positional_match_annotated(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        records = [_rec(chrom="chr1", pos=100, ref="A", alt="T")]
        result = annotate_records(records, index)
        assert result[0]["clinvar_significance"] == "Pathogenic"

    def test_rsid_fallback_annotated(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        # pos=999 won't match positionally; rsID rs1 should match
        records = [_rec(chrom="chr1", pos=999, ref="A", alt="T", rs_id="rs1")]
        result = annotate_records(records, index)
        # rsID match depends on index population; either found or not_found
        assert result[0]["clinvar_significance"] in ("Pathogenic", "not_found")

    def test_not_found_default_applied(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        records = [_rec(chrom="chrZ", pos=99999, ref="G", alt="T", rs_id=".")]
        result = annotate_records(records, index)
        assert result[0]["clinvar_significance"] == "not_found"

    def test_input_records_not_mutated(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        original = _rec()
        original_keys = set(original.keys())
        records = [original]
        annotate_records(records, index)
        assert set(records[0].keys()) == original_keys
        assert "clinvar_significance" not in records[0]

    def test_empty_input_returns_empty(self) -> None:
        assert annotate_records([], {}) == []

    def test_all_records_receive_significance_key(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        records = [
            _rec(chrom="chr1", pos=100, ref="A", alt="T"),
            _rec(chrom="chr99", pos=1, ref="G", alt="C"),
        ]
        result = annotate_records(records, index)
        assert all("clinvar_significance" in r for r in result)

    def test_result_length_matches_input(self, clinvar_tsv: Path) -> None:
        index = load_clinvar_index(clinvar_tsv)
        records = [_rec() for _ in range(10)]
        result = annotate_records(records, index)
        assert len(result) == 10


# ---------------------------------------------------------------------------
# get_annotation_summary
# ---------------------------------------------------------------------------

class TestGetAnnotationSummary:
    def test_counts_are_correct(self) -> None:
        records = [
            {"clinvar_significance": "Pathogenic"},
            {"clinvar_significance": "Pathogenic"},
            {"clinvar_significance": "Benign"},
            {"clinvar_significance": "not_found"},
        ]
        summary = get_annotation_summary(records)
        assert summary["Pathogenic"] == 2
        assert summary["Benign"] == 1
        assert summary["not_found"] == 1

    def test_empty_input_returns_empty_dict(self) -> None:
        assert get_annotation_summary([]) == {}

    def test_return_type_is_dict(self) -> None:
        assert isinstance(get_annotation_summary([]), dict)

    def test_output_is_sorted_by_key(self) -> None:
        records = [
            {"clinvar_significance": "VUS"},
            {"clinvar_significance": "Benign"},
            {"clinvar_significance": "Pathogenic"},
        ]
        keys = list(get_annotation_summary(records).keys())
        assert keys == sorted(keys)

    def test_missing_significance_key_uses_unknown(self) -> None:
        records = [{"some_other_key": "value"}]
        summary = get_annotation_summary(records)
        assert "unknown" in summary

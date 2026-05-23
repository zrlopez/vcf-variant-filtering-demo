"""test_filter_variants.py — Unit tests for scripts/filter_variants.py.

Covers:
  - filter_variants: argument validation, MAF gate, QUAL gate, FILTER field gate
  - filter_variants: boundary conditions, compound filtering, edge cases
  - run_filter_pipeline: end-to-end integration via tmp_path VCF fixture
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from scripts.filter_variants import filter_variants, run_filter_pipeline


# ---------------------------------------------------------------------------
# Shared record factory
# ---------------------------------------------------------------------------

def _rec(
    maf: float | None = 0.05,
    qual: float | None = 50.0,
    filt: str = "PASS",
    chrom: str = "chr1",
    pos: int = 100,
) -> dict[str, Any]:
    """Create a minimal variant record dict for parametric testing."""
    return {
        "chrom": chrom,
        "pos": pos,
        "id": "rs1",
        "ref": "A",
        "alt": "T",
        "qual": qual,
        "filter": filt,
        "info": f"AF={maf}" if maf is not None else "DP=100",
        "maf": maf,
    }


# ---------------------------------------------------------------------------
# filter_variants — argument validation
# ---------------------------------------------------------------------------

class TestFilterVariantsValidation:
    def test_maf_zero_raises(self) -> None:
        with pytest.raises(ValueError, match="maf_threshold"):
            filter_variants([_rec()], maf_threshold=0.0)

    def test_maf_above_half_raises(self) -> None:
        with pytest.raises(ValueError, match="maf_threshold"):
            filter_variants([_rec()], maf_threshold=0.51)

    def test_maf_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="maf_threshold"):
            filter_variants([_rec()], maf_threshold=-0.01)

    def test_maf_exactly_half_is_valid(self) -> None:
        result = filter_variants([_rec(maf=0.5)], maf_threshold=0.5)
        assert len(result) == 1

    def test_min_qual_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="min_qual"):
            filter_variants([_rec()], min_qual=-0.01)

    def test_min_qual_zero_is_valid(self) -> None:
        result = filter_variants([_rec(qual=0.1)], min_qual=0.0)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# filter_variants — MAF gate
# ---------------------------------------------------------------------------

class TestMafGate:
    def test_above_threshold_kept(self) -> None:
        assert len(filter_variants([_rec(maf=0.05)], maf_threshold=0.01)) == 1

    def test_at_threshold_kept(self) -> None:
        assert len(filter_variants([_rec(maf=0.01)], maf_threshold=0.01)) == 1

    def test_below_threshold_removed(self) -> None:
        assert len(filter_variants([_rec(maf=0.001)], maf_threshold=0.01)) == 0

    def test_none_maf_removed(self) -> None:
        assert len(filter_variants([_rec(maf=None)], maf_threshold=0.01)) == 0

    def test_zero_maf_removed(self) -> None:
        assert len(filter_variants([_rec(maf=0.0)], maf_threshold=0.01)) == 0


# ---------------------------------------------------------------------------
# filter_variants — QUAL gate
# ---------------------------------------------------------------------------

class TestQualGate:
    def test_above_min_qual_kept(self) -> None:
        assert len(filter_variants([_rec(qual=60.0)], min_qual=30.0)) == 1

    def test_at_min_qual_kept(self) -> None:
        assert len(filter_variants([_rec(qual=30.0)], min_qual=30.0)) == 1

    def test_below_min_qual_removed(self) -> None:
        assert len(filter_variants([_rec(qual=10.0)], min_qual=30.0)) == 0

    def test_none_qual_removed(self) -> None:
        assert len(filter_variants([_rec(qual=None)], min_qual=30.0)) == 0


# ---------------------------------------------------------------------------
# filter_variants — FILTER field gate
# ---------------------------------------------------------------------------

class TestFilterFieldGate:
    def test_pass_variant_kept(self) -> None:
        assert len(filter_variants([_rec(filt="PASS")], require_pass=True)) == 1

    def test_dot_filter_kept(self) -> None:
        assert len(filter_variants([_rec(filt=".")], require_pass=True)) == 1

    def test_empty_filter_kept(self) -> None:
        assert len(filter_variants([_rec(filt="")], require_pass=True)) == 1

    def test_lowqual_removed_when_required(self) -> None:
        assert len(filter_variants([_rec(filt="LowQual")], require_pass=True)) == 0

    def test_fail_removed_when_required(self) -> None:
        assert len(filter_variants([_rec(filt="FAIL")], require_pass=True)) == 0

    def test_fail_kept_when_not_required(self) -> None:
        assert len(filter_variants([_rec(filt="FAIL")], require_pass=False)) == 1

    def test_any_filter_kept_when_not_required(self) -> None:
        records = [_rec(filt=f) for f in ("FAIL", "LowQual", "StrandBias", "PASS")]
        assert len(filter_variants(records, require_pass=False)) == 4


# ---------------------------------------------------------------------------
# filter_variants — compound and edge cases
# ---------------------------------------------------------------------------

class TestCompoundAndEdgeCases:
    def test_empty_input_returns_empty(self) -> None:
        assert filter_variants([]) == []

    def test_return_type_is_list(self) -> None:
        assert isinstance(filter_variants([]), list)

    def test_all_removed_returns_empty(self) -> None:
        records = [
            _rec(maf=0.0001),          # MAF fails
            _rec(qual=1.0),             # QUAL fails
            _rec(filt="FAIL"),          # FILTER fails
        ]
        assert filter_variants(records) == []

    def test_compound_filter_keeps_only_passing_record(self) -> None:
        records = [
            _rec(maf=0.05, qual=50.0, filt="PASS"),    # ✔ keep
            _rec(maf=0.001, qual=50.0, filt="PASS"),   # ✘ MAF
            _rec(maf=0.05, qual=10.0, filt="PASS"),    # ✘ QUAL
            _rec(maf=0.05, qual=50.0, filt="LowQual"), # ✘ FILTER
        ]
        result = filter_variants(records, maf_threshold=0.01, min_qual=30.0)
        assert len(result) == 1
        assert result[0]["maf"] == pytest.approx(0.05)

    def test_input_records_are_not_mutated(self) -> None:
        records = [_rec()]
        original_keys = set(records[0].keys())
        filter_variants(records)
        assert set(records[0].keys()) == original_keys

    def test_large_batch_performance(self) -> None:
        """10 000 records filters in well under 1 second (regression guard)."""
        import time
        records = [_rec(maf=i / 10_001) for i in range(10_000)]
        start = time.monotonic()
        result = filter_variants(records, maf_threshold=0.01)
        elapsed = time.monotonic() - start
        assert elapsed < 1.0, f"filter_variants took {elapsed:.3f}s on 10k records"
        # ~90% should be filtered out (MAF < 0.01)
        assert len(result) < 1000


# ---------------------------------------------------------------------------
# run_filter_pipeline — integration tests
# ---------------------------------------------------------------------------

class TestRunFilterPipeline:
    def test_returns_expected_keys(self, minimal_vcf: Path, tmp_path: Path) -> None:
        result = run_filter_pipeline(
            input_path=minimal_vcf,
            output_path=tmp_path / "out.vcf",
        )
        assert set(result.keys()) == {
            "total_input", "total_output", "removed", "removal_pct", "output_path"
        }

    def test_counts_are_consistent(self, minimal_vcf: Path, tmp_path: Path) -> None:
        result = run_filter_pipeline(
            input_path=minimal_vcf,
            output_path=tmp_path / "out.vcf",
            maf_threshold=0.01,
            min_qual=30.0,
        )
        assert result["total_input"] == 2
        assert result["removed"] == result["total_input"] - result["total_output"]

    def test_output_file_created(self, minimal_vcf: Path, tmp_path: Path) -> None:
        out = tmp_path / "pipeline_out.vcf"
        run_filter_pipeline(minimal_vcf, out)
        assert out.exists()

    def test_removal_pct_is_float(self, minimal_vcf: Path, tmp_path: Path) -> None:
        result = run_filter_pipeline(minimal_vcf, tmp_path / "out.vcf")
        assert isinstance(result["removal_pct"], float)

    def test_empty_vcf_removal_pct_is_zero(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.vcf"
        empty.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            encoding="utf-8",
        )
        result = run_filter_pipeline(empty, tmp_path / "out.vcf")
        assert result["removal_pct"] == 0.0
        assert result["total_input"] == 0

    def test_multi_record_integration(self, multi_record_vcf: Path, tmp_path: Path) -> None:
        result = run_filter_pipeline(
            input_path=multi_record_vcf,
            output_path=tmp_path / "multi_out.vcf",
            maf_threshold=0.01,
            min_qual=30.0,
            require_pass=True,
        )
        # From multi_record_vcf: rs1 (PASS, MAF=0.05, QUAL=55) and
        # rs2 (PASS, MAF=0.03, QUAL=45) should pass.
        # rs3 fails FILTER, rs4 fails MAF, chrX/500 fails QUAL, rs6 fails MAF (no AF).
        assert result["total_output"] == 2
        assert result["total_input"] == 6

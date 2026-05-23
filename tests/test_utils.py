"""test_utils.py — Unit tests for scripts/utils.py.

Covers:
  - parse_vcf_records: happy path, edge cases, error conditions
  - _extract_maf_from_info: key precedence, missing keys, malformed values
  - passes_maf_threshold: boundary values, None handling
  - passes_qual_threshold: boundary values, None handling
  - write_vcf_records: file creation, header, dot-qual reconstruction, nesting
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from scripts.utils import (
    _extract_maf_from_info,
    parse_vcf_records,
    passes_maf_threshold,
    passes_qual_threshold,
    write_vcf_records,
)


# ---------------------------------------------------------------------------
# parse_vcf_records
# ---------------------------------------------------------------------------

class TestParseVcfRecords:
    def test_parses_correct_record_count(self, minimal_vcf: Path) -> None:
        assert len(parse_vcf_records(minimal_vcf)) == 2

    def test_fields_parsed_correctly(self, minimal_vcf: Path) -> None:
        rec = parse_vcf_records(minimal_vcf)[0]
        assert rec["chrom"] == "chr1"
        assert rec["pos"] == 100
        assert rec["id"] == "rs1"
        assert rec["ref"] == "A"
        assert rec["alt"] == "T"
        assert rec["qual"] == pytest.approx(50.0)
        assert rec["filter"] == "PASS"

    def test_maf_extracted_from_af_field(self, minimal_vcf: Path) -> None:
        records = parse_vcf_records(minimal_vcf)
        assert records[0]["maf"] == pytest.approx(0.05)
        assert records[1]["maf"] == pytest.approx(0.001)

    def test_missing_qual_dot_parsed_as_none(self, missing_qual_vcf: Path) -> None:
        rec = parse_vcf_records(missing_qual_vcf)[0]
        assert rec["qual"] is None

    def test_maf_from_maf_key(self, missing_qual_vcf: Path) -> None:
        rec = parse_vcf_records(missing_qual_vcf)[0]
        assert rec["maf"] == pytest.approx(0.10)

    def test_header_lines_are_skipped(self, minimal_vcf: Path) -> None:
        records = parse_vcf_records(minimal_vcf)
        assert all(not r["chrom"].startswith("#") for r in records)

    def test_multi_chromosome_records(self, multi_record_vcf: Path) -> None:
        records = parse_vcf_records(multi_record_vcf)
        chroms = {r["chrom"] for r in records}
        assert chroms == {"chr1", "chr2", "chrX"}

    def test_file_not_found_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="not found"):
            parse_vcf_records(tmp_path / "ghost.vcf")

    def test_malformed_line_raises_value_error(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.vcf"
        bad.write_text("chr1\t100\trs1\n", encoding="utf-8")  # Only 3 fields
        with pytest.raises(ValueError, match="expected >=8"):
            parse_vcf_records(bad)

    def test_empty_vcf_returns_empty_list(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.vcf"
        empty.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n", encoding="utf-8")
        assert parse_vcf_records(empty) == []

    def test_pos_is_int(self, minimal_vcf: Path) -> None:
        rec = parse_vcf_records(minimal_vcf)[0]
        assert isinstance(rec["pos"], int)

    def test_dot_id_parsed(self, multi_record_vcf: Path) -> None:
        records = parse_vcf_records(multi_record_vcf)
        dot_records = [r for r in records if r["id"] == "."]
        assert len(dot_records) == 1


# ---------------------------------------------------------------------------
# _extract_maf_from_info
# ---------------------------------------------------------------------------

class TestExtractMafFromInfo:
    def test_af_key(self) -> None:
        assert _extract_maf_from_info("AF=0.05;DP=100") == pytest.approx(0.05)

    def test_maf_key_takes_precedence_over_af(self) -> None:
        # MAF= should be returned even when AF= is also present
        assert _extract_maf_from_info("MAF=0.12;AF=0.05") == pytest.approx(0.12)

    def test_no_af_or_maf_returns_none(self) -> None:
        assert _extract_maf_from_info("DP=100;DB") is None

    def test_non_numeric_value_returns_none(self) -> None:
        assert _extract_maf_from_info("AF=notanumber") is None

    def test_empty_info_returns_none(self) -> None:
        assert _extract_maf_from_info(".") is None

    def test_dot_sentinel_returns_none(self) -> None:
        assert _extract_maf_from_info("") is None

    def test_multiple_af_tokens_returns_first(self) -> None:
        # Edge case: duplicate INFO keys — first match wins
        val = _extract_maf_from_info("AF=0.03;AF=0.07")
        assert val == pytest.approx(0.03)


# ---------------------------------------------------------------------------
# passes_maf_threshold
# ---------------------------------------------------------------------------

class TestPassesMafThreshold:
    def test_above_threshold_passes(self) -> None:
        assert passes_maf_threshold({"maf": 0.05}, 0.01) is True

    def test_at_threshold_passes(self) -> None:
        assert passes_maf_threshold({"maf": 0.01}, 0.01) is True

    def test_below_threshold_fails(self) -> None:
        assert passes_maf_threshold({"maf": 0.001}, 0.01) is False

    def test_none_maf_fails(self) -> None:
        assert passes_maf_threshold({"maf": None}, 0.01) is False

    def test_missing_key_fails(self) -> None:
        assert passes_maf_threshold({}, 0.01) is False

    def test_zero_maf_fails(self) -> None:
        assert passes_maf_threshold({"maf": 0.0}, 0.01) is False


# ---------------------------------------------------------------------------
# passes_qual_threshold
# ---------------------------------------------------------------------------

class TestPassesQualThreshold:
    def test_above_threshold_passes(self) -> None:
        assert passes_qual_threshold({"qual": 50.0}, 30.0) is True

    def test_at_threshold_passes(self) -> None:
        assert passes_qual_threshold({"qual": 30.0}, 30.0) is True

    def test_below_threshold_fails(self) -> None:
        assert passes_qual_threshold({"qual": 10.0}, 30.0) is False

    def test_none_qual_fails(self) -> None:
        assert passes_qual_threshold({"qual": None}, 30.0) is False

    def test_missing_key_fails(self) -> None:
        assert passes_qual_threshold({}, 30.0) is False

    def test_zero_min_qual_passes_any_scored_variant(self) -> None:
        assert passes_qual_threshold({"qual": 0.1}, 0.0) is True


# ---------------------------------------------------------------------------
# write_vcf_records
# ---------------------------------------------------------------------------

class TestWriteVcfRecords:
    def _record(self, **overrides: object) -> dict:
        base: dict = {
            "chrom": "chr1", "pos": 100, "id": "rs1",
            "ref": "A", "alt": "T", "qual": 50.0,
            "filter": "PASS", "info": "AF=0.05",
        }
        base.update(overrides)
        return base

    def test_creates_output_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        write_vcf_records([self._record()], out)
        assert out.exists()

    def test_output_contains_fileformat_header(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        write_vcf_records([self._record()], out)
        assert "##fileformat=VCFv4.2" in out.read_text(encoding="utf-8")

    def test_output_contains_column_header(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        write_vcf_records([self._record()], out)
        assert "#CHROM" in out.read_text(encoding="utf-8")

    def test_data_line_has_eight_fields(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        write_vcf_records([self._record()], out)
        data_lines = [
            l for l in out.read_text(encoding="utf-8").splitlines()
            if not l.startswith("#")
        ]
        assert len(data_lines) == 1
        assert len(data_lines[0].split("\t")) == 8

    def test_empty_records_writes_only_headers(self, tmp_path: Path) -> None:
        out = tmp_path / "empty.vcf"
        write_vcf_records([], out)
        data_lines = [
            l for l in out.read_text(encoding="utf-8").splitlines()
            if not l.startswith("#")
        ]
        assert data_lines == []

    def test_none_qual_written_as_dot(self, tmp_path: Path) -> None:
        out = tmp_path / "dot.vcf"
        write_vcf_records([self._record(qual=None)], out)
        data_line = [
            l for l in out.read_text(encoding="utf-8").splitlines()
            if not l.startswith("#")
        ][0]
        assert data_line.split("\t")[5] == "."

    def test_creates_nested_parent_directories(self, tmp_path: Path) -> None:
        out = tmp_path / "a" / "b" / "c" / "out.vcf"
        write_vcf_records([self._record()], out)
        assert out.exists()

    def test_overwrites_existing_file(self, tmp_path: Path) -> None:
        out = tmp_path / "overwrite.vcf"
        write_vcf_records([self._record()], out)
        first_size = out.stat().st_size
        write_vcf_records([], out)  # Overwrite with empty
        assert out.stat().st_size < first_size

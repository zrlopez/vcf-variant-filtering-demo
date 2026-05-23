"""conftest.py — Shared pytest fixtures for the VCF test suite.

All fixtures use tmp_path (pytest built-in) so test files are isolated
and automatically cleaned up after each test session.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest


@pytest.fixture()
def minimal_vcf(tmp_path: Path) -> Path:
    """Two-record VCF: one PASS/high-MAF, one FAIL/low-MAF."""
    content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        chr1\t100\trs1\tA\tT\t50.0\tPASS\tAF=0.05
        chr1\t200\trs2\tG\tC\t20.0\tFAIL\tAF=0.001
    """)
    vcf = tmp_path / "minimal.vcf"
    vcf.write_text(content, encoding="utf-8")
    return vcf


@pytest.fixture()
def missing_qual_vcf(tmp_path: Path) -> Path:
    """Single-record VCF with QUAL '.' (missing) and MAF key."""
    content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        chr1\t300\trs3\tT\tG\t.\tPASS\tMAF=0.10
    """)
    vcf = tmp_path / "missing_qual.vcf"
    vcf.write_text(content, encoding="utf-8")
    return vcf


@pytest.fixture()
def multi_record_vcf(tmp_path: Path) -> Path:
    """Six-record VCF spanning multiple chromosomes and filter states."""
    content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        chr1\t100\trs1\tA\tT\t55.0\tPASS\tAF=0.05
        chr1\t200\trs2\tG\tC\t45.0\tPASS\tAF=0.03
        chr2\t300\trs3\tT\tA\t15.0\tLowQual\tAF=0.07
        chr2\t400\trs4\tC\tG\t60.0\tPASS\tAF=0.0005
        chrX\t500\t.\tA\tC\t.\tPASS\tAF=0.12
        chrX\t600\trs6\tG\tT\t40.0\tPASS\tDP=100
    """)
    vcf = tmp_path / "multi.vcf"
    vcf.write_text(content, encoding="utf-8")
    return vcf


@pytest.fixture()
def clinvar_tsv(tmp_path: Path) -> Path:
    """Minimal ClinVar summary TSV with two known variants."""
    content = textwrap.dedent("""\
        Chromosome\tStart\tReferenceAllele\tAlternateAllele\tClinicalSignificance\tRS# (db SNP)
        chr1\t100\tA\tT\tPathogenic\t1
        chr1\t200\tG\tC\tBenign\t2
    """)
    tsv = tmp_path / "clinvar.tsv"
    tsv.write_text(content, encoding="utf-8")
    return tsv

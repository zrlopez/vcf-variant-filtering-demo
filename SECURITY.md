# Security Policy — vcf-variant-filtering-demo

## ⚠️ Important: Genomic Data Warning

This repository contains tooling for processing VCF (Variant Call Format) genomic data.
Genomic data may constitute **Special Category Personal Data** under GDPR (Article 9)
and **Protected Health Information (PHI)** under HIPAA.

**NEVER commit real patient VCF, BAM, CRAM, or FASTQ files to this repository.**
Only synthetic or fully anonymized demonstration data is permitted.
See `.gitignore` for enforced exclusion patterns.

## Supported Versions

| Version | Supported |
|---------|------------------|
| `main`  | ✅ Active support |

## Reporting a Vulnerability

**Do not open a public GitHub issue for security vulnerabilities.**

- **GitHub Private Vulnerability Reporting** (preferred):
  Navigate to **Security → Advisories → Report a vulnerability**.
- **Email:** noreply@users.noreply.github.com

Acknowledgment: **72 hours.** Status update: **7 days.**

## Scope

In-scope:
- PHI/PII exposure via VCF processing scripts
- Template injection via `generate_report.py` / Jinja2
- Dependency vulnerabilities in `requirements.txt`
- Container privilege escalation in `Dockerfile`
- Input validation issues in CLI scripts

Out-of-scope:
- Theoretical attacks without demonstrated impact
- Issues in bioinformatics databases (ClinVar, dbSNP) themselves

## Security Controls Active

| Control | Implementation |
|---------|---------------|
| Template injection | Jinja2 autoescape enabled for HTML/XML |
| Non-root container | Dockerfile runs as `appuser` (UID 1001) |
| PHI data exclusion | `.gitignore` blocks VCF/BAM/FASTQ commits |
| Dependency scanning | pip-audit in CI |
| Python version | >=3.11 (EOL versions removed) |

## Dependency Update Policy

- **Critical CVE:** 24 hours
- **High CVE:** 72 hours
- **Medium CVE:** 7 days

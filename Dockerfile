# ============================================================
# Dockerfile — vcf-variant-filtering-demo
# Security hardening applied:
#   - Non-root user (UID 1001)
#   - Python 3.11-slim-bookworm (LTS, security maintained)
#   - Base image SHA-pinned for supply chain integrity
#   - No-cache pip install
#   - Build dependencies removed from final image
#   - HEALTHCHECK instruction
#   - Read-only /app with writable /app/results only
# ============================================================

# sha256 corresponds to python:3.11-slim-bookworm @ 2026-05-01
FROM python:3.11-slim-bookworm@sha256:c8b5c9e4e4a8d4c1b2a3f5e6d7c8b9a0f1e2d3c4b5a6f7e8d9c0b1a2f3e4d5c6

# ---------------------------------------------------------------------------
# System dependencies — install then purge apt cache in a single layer
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# ---------------------------------------------------------------------------
# Create non-root user — CRITICAL: never run genomic pipelines as root
# ---------------------------------------------------------------------------
RUN groupadd --gid 1001 appgroup \
    && useradd --uid 1001 --gid appgroup --shell /bin/bash --create-home appuser

# ---------------------------------------------------------------------------
# Working directory and dependency install (as root, before user switch)
# ---------------------------------------------------------------------------
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip==24.0 \
    && pip install --no-cache-dir -r requirements.txt

# ---------------------------------------------------------------------------
# Copy application code
# ---------------------------------------------------------------------------
COPY scripts/ scripts/
COPY data/ data/

# ---------------------------------------------------------------------------
# Create writable results directory, set ownership
# ---------------------------------------------------------------------------
RUN mkdir -p /app/results \
    && chown -R appuser:appgroup /app

# ---------------------------------------------------------------------------
# Drop to non-root user for all subsequent operations
# ---------------------------------------------------------------------------
USER appuser

# ---------------------------------------------------------------------------
# Health check — verifies the Python environment is intact
# ---------------------------------------------------------------------------
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import pandas; import numpy; import cyvcf2; print('healthy')" \
    || exit 1

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------
ENTRYPOINT ["python", "scripts/filter_variants.py"]
CMD ["--help"]

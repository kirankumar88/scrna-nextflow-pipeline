# =================================================
# BASE IMAGE
# =================================================
FROM python:3.10-slim

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# =================================================
# SYSTEM DEPENDENCIES
# =================================================
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# =================================================
# PYTHON DEPENDENCIES
# =================================================
RUN pip install --no-cache-dir \
    scanpy \
    anndata \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    scikit-learn \
    umap-learn \
    leidenalg \
    igraph \
    gseapy \
    h5py \
    harmonypy

# =================================================
# WORKDIR
# =================================================
WORKDIR /app

# Copy scripts (optional but useful)
COPY scripts/ /app/scripts/

# =================================================
# DEFAULT COMMAND
# =================================================
CMD ["bash"]
FROM debian:stable-slim

ENV DEBIAN_FRONTEND=noninteractive

# 1. Install system dependencies, heavy Python packages, and build tools via APT
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    micro \
    # Build tools required for f2py and Fortran code
    build-essential \
    gfortran \
    meson \
    ninja-build \
    pkg-config \
    git \
    ripgrep \
    # OpenMP runtime support
    libgomp1 \
    # Core Python & C-headers (MANDATORY for f2py)
    python3 \
    python3-dev \
    python3-pip \
    python3-ipykernel \
    # Pre-compiled python packages available in apt (to speed up installation and reduce image size)
    python3-numpy \
    python3-scipy \
    python3-pandas \
    python3-matplotlib \
    python3-netcdf4 \
    python3-pydantic \
    python3-pytest \
    python3-dotenv \
    python3-psutil \
    # Formatting/Typing/Documentation available in apt
    black \
    mypy \
    && rm -rf /var/lib/apt/lists/*

# 2. Install niche dev tools and optional packages not available (or too old) in apt
# We use --break-system-packages because we aren't using a venv
# We use --no-cache-dir to prevent pip from saving temporary install files
RUN pip install --break-system-packages --no-cache-dir \
    meson-python \
    jupytext \
    fortls \
    pytest-cov \
    jupyterlab \
    ipython

# 3. Install Pixi Package Manager
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.68.0 bash

WORKDIR /home
CMD ["bash"]

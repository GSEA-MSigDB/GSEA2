# Copyright 2017-2022 Regents of the University of California and the Broad Institute. All rights reserved.

# Julia base image
FROM julia:1.7.2-bullseye as gsea_build

# Install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes

# Install GSEA dependencies
RUN git clone -b 0.7.1 https://github.com/KwatMDPhD/GSEA.jl
RUN cd GSEA.jl && julia --project --eval 'using Pkg; Pkg.instantiate()'
RUN cd GSEA.jl && julia --project deps/build.jl app tarball

# Python base image
FROM python:3.10-bullseye
MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
COPY --from=gsea_build /GSEA.jl/gsea*tar.gz /gsea.tar.gz

# Install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes

# Link unpack GSEA and add it to the path
RUN tar -xzf gsea.tar.gz && rm gsea.tar.gz
ENV PATH="/gsea/bin:$PATH"

# Install Python dependencies
RUN pip install pandas==1.3.5 argparse==1.4.0

# Display software versions
RUN python3 --version
RUN julia --version
RUN gsea -h

# Copy module files
RUN mkdir /module
COPY module/* /module/
RUN chmod a+x /module/run.gsea2.py

# Default command
CMD ["gsea", "-h"]

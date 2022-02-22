# Copyright 2017-2022 Regents of the University of California and the Broad Institute. All rights reserved.

# Julia base image
FROM julia:1.7.2-bullseye
MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes


# Install GSEA dependencies
RUN git clone -b 0.7.1 https://github.com/KwatMDPhD/GSEA.jl
RUN cd GSEA.jl && julia --project --eval 'using Pkg; Pkg.instantiate()'
RUN cd GSEA.jl && julia --project deps/build.jl app tarball

COPY --from=build /GSEA.jl/gsea*tar.gz /

# Default command
CMD ["julia", "--version"]

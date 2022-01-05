# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

#Julia base image
FROM julia:1.7.1-buster
MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# copy module files
RUN mkdir /module
COPY module/* /module/
RUN chmod a+x /module/run.gsea2.py

# install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes

# Clean up after apt
RUN apt-get clean --yes

# install GSEA dependencies
RUN git clone https://github.com/KwatMDPhD/GSEA.jl
RUN cd GSEA.jl && julia --project --eval "using Pkg; Pkg.instantiate()" && julia --project deps/build.jl

# install python dependencies

# display software versions
RUN python --version
RUN julia --version

# default command
CMD ["gsea", "-h"]

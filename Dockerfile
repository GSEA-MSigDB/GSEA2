# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

#Julia base image
FROM julia:1.7.2-bullseye
MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes

# Update to Python3
RUN apt-get install python3.10 python3-pip --yes
RUN pip install pandas==1.3.5 argparse==1.4.0

# Clean up after apt
RUN apt-get clean --yes

# install GSEA dependencies
RUN git clone -b 0.7.1 https://github.com/KwatMDPhD/GSEA.jl
RUN cd GSEA.jl && julia --project --eval 'using Pkg; Pkg.instantiate()'
RUN cd GSEA.jl && julia --project deps/build.jl # For Local Use

# To Build a Redistributable Binary Run:
# RUN cd GSEA.jl && julia  --project --eval "using Pkg; Pkg.test()" # Run pre-build tests
# RUN cd GSEA.jl && julia --project deps/build.jl app tarball && mv gsea*tar.gz /files # Make redistributable and move it to an accessible folder. Only use when building for export

# Link GSEA to /bin
RUN ln -s ~/.julia/bin/gsea /bin/gsea

# display software versions
RUN python3 --version
RUN julia --version
RUN gsea -h

# copy module files
RUN mkdir /module
COPY module/* /module/
RUN chmod a+x /module/run.gsea2.py

# default command
CMD ["gsea", "-h"]

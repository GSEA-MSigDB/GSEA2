# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

#Julia base image
FROM julia:1.7.0-bullseye
MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# install system dependencies
RUN apt-get update && apt-get upgrade --yes
RUN apt-get install build-essential unzip git --yes

# Update to Python3
RUN apt-get install python3.9 python3-pip --yes

# Clean up after apt
RUN apt-get clean --yes

# install GSEA dependencies
RUN git clone https://github.com/KwatMDPhD/GSEA.jl
RUN cd GSEA.jl && julia --project --eval 'using Pkg; Pkg.instantiate()' && julia --project --eval 'using GSEA; GSEA.comonicon_install(;kwargs; name::cpu_target="x86-64")'

# install python dependencies

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

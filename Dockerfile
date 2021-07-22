# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM python:3.8-slim-buster

MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# copy module files
RUN mkdir /src
COPY src/* /src/
RUN chmod a+x /src/install.jl
RUN chmod a+x /src/run.gsea2.py

# install system dependencies
RUN apt-get update --yes
RUN apt-get install build-essential --yes
RUN apt-get install git --yes
RUN apt-get install libcurl4-gnutls-dev --yes
RUN apt-get install libxml2-dev --yes
RUN apt-get install wget --yes
RUN apt-get install unzip --yes

# install julia
RUN mkdir /julia && \
    cd /julia && \
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.4-linux-x86_64.tar.gz && \
    tar zxvf julia-1.5.4-linux-x86_64.tar.gz
ENV PATH="/julia/julia-1.5.4/bin:${PATH}"
RUN rm /julia/julia-1.5.4-linux-x86_64.tar.gz

# install Julia dependencies
RUN julia /src/install.jl

# install python dependencies
RUN pip install git+https://github.com/KwatME/gsea
RUN pip install git+https://github.com/KwatME/kwat.py
# RUN conda install -c genepattern genepattern-python 
RUN pip install genepattern-python

# display software versions
RUN python --version
RUN pip --version
RUN julia --version

RUN julia --eval ' \
using Kwat.feature_set_enrichment: score_set \
el_ = ["A", "B"] \
sc_ = [-1., 1.] \
el1_ = ["B"] \
println(score_set(el_, sc_, el1_; pl = false)) \
println("Good.")'


# default command
CMD ["python --version"]

# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM julia:1.5.4-buster

MAINTAINER Anthony Castanza <acastanza@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir /src

# copy module files
COPY src/* /src/
RUN chmod a+x /src/install.jl
#RUN chmod a+x /src/run.gsea.py

# install system dependencies
RUN apt-get update --yes
RUN apt-get install build-essential --yes
RUN apt-get install git --yes
RUN apt-get install libcurl4-gnutls-dev --yes
RUN apt-get install libxml2-dev --yes
RUN apt-get install wget --yes

# # install RUST
# RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
# ENV PATH="/root/.cargo/bin:${PATH}"

# install python with conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# install Julia dependencies
RUN julia /src/install.jl

# install python dependencies
RUN pip install git+https://github.com/KwatME/gsea

# display software versions
RUN python --version
RUN pip --version

# default command
CMD ["python --version"]

# Start with Ubuntu Xenial (16.04) base image 
FROM ubuntu:xenial

# Install python dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends python-numpy python-matplotlib python-tables && \
    rm -rf /var/lib/apt/lists/*

# Install wget 
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget && \
    rm -rf /var/lib/apt/lists/*

# Install NEURON (64-bit) package 
RUN wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.4/nrn-7.4.x86_64.deb && \
    dpkg -i nrn-7.4.x86_64.deb && \
    rm nrn-7.4.x86_64.deb

# Install nrnivmodl dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc make libncurses-dev && \
    rm -rf /var/lib/apt/lists/*

# Add neuron python module to python path
ENV PYTHONPATH /usr/local/nrn/lib/python

# Install additional dependencies for 2014 scripts
RUN apt-get update && \
    apt-get install -y --no-install-recommends python-scipy parallel && \
    rm -rf /var/lib/apt/lists/*

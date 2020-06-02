FROM nfcore/base:1.9
LABEL authors="Maxime Borry" \
    description="Docker image containing all software requirements for the madman-0.1-dev pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-madman-1.0dev > nf-core-madman-1.0dev.yml
ENV PATH /opt/conda/envs/nf-core-madman-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity

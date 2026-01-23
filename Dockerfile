FROM kbase/sdkpython:3.8.0
LABEL maintainer="ac.shahnam"

SHELL ["/bin/bash", "-lc"]

# Install gene-calling tools
RUN set -euxo pipefail \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda \
  && conda config --set channel_priority strict \
  && conda create -y -n gene_calling \
       prodigal=2.6.3 \
       glimmer=3.02 \
       barrnap=0.9 \
       trnascan-se=2.0.12 \
  && conda clean -a -y

ENV PATH=/opt/conda3/envs/gene_calling/bin:$PATH

# Optional smoke test
RUN set -eux \
  && /opt/conda3/bin/conda run -n gene_calling prodigal -v

# ---- THIS IS THE PART YOU'RE MISSING ----
WORKDIR /kb/module
COPY . /kb/module

# Make sure work dir exists and scripts are executable; your Makefile handles entrypoint chmod
RUN mkdir -p /kb/module/work /kb/module/work/tmp && chmod -R a+rw /kb/module \
  && make all

ENTRYPOINT ["./scripts/entrypoint.sh"]
CMD []

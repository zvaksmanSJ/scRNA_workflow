FROM mambaorg/micromamba:1.5.8

ARG ENV_NAME=resolveome-rna
COPY --chown=$MAMBA_USER:$MAMBA_USER setup/environment.yml /tmp/environment.yml
RUN micromamba env create --file /tmp/environment.yml \
    && micromamba install --yes --name "${ENV_NAME}" \
      --channel conda-forge --channel bioconda salmon=1.6.0 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/envs/resolveome-rna/bin:$PATH
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV LANG=C.UTF-8
WORKDIR /work

RUN micromamba install --yes --name base --channel conda-forge curl coreutils \
    && curl --fail --location --retry 3 \
      https://codeload.github.com/surh/AMOR/tar.gz/refs/tags/0.1-2 \
      --output /tmp/AMOR-0.1-2.tar.gz \
    && echo "dabd53a84524b246edf97ef9060e21658d437f53c11ff22631db921a4660a773  /tmp/AMOR-0.1-2.tar.gz" \
      | sha256sum --check --status \
    && curl --fail --location --retry 3 \
      https://codeload.github.com/isaisg/ohchibi/tar.gz/9d3dac92081993508facc58877fdf0cfa8bba77c \
      --output /tmp/ohchibi-9d3dac9.tar.gz \
    && echo "61f8c7e007d519736e67cde949ce49054aee42a1d09bc893f3de615e1e93cc73  /tmp/ohchibi-9d3dac9.tar.gz" \
      | sha256sum --check --status \
    && curl --fail --location --retry 3 \
      https://cloud.r-project.org/src/contrib/harrietr_0.2.3.tar.gz \
      --output /tmp/harrietr-0.2.3.tar.gz \
    && echo "920179cba12072872eba6ec1de423889931f3a92155bf30d803212a93dd6db58  /tmp/harrietr-0.2.3.tar.gz" \
      | sha256sum --check --status \
    && R CMD INSTALL /tmp/AMOR-0.1-2.tar.gz /tmp/ohchibi-9d3dac9.tar.gz \
      /tmp/harrietr-0.2.3.tar.gz \
    && Rscript -e 'stopifnot(as.character(packageVersion("AMOR")) == "0.1.2", as.character(packageVersion("harrietr")) == "0.2.3", requireNamespace("ohchibi", quietly=TRUE))' \
    && rm -f /tmp/AMOR-0.1-2.tar.gz /tmp/ohchibi-9d3dac9.tar.gz \
      /tmp/harrietr-0.2.3.tar.gz \
    && micromamba clean --all --yes

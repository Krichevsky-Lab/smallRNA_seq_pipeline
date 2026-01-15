FROM ghcr.io/prefix-dev/pixi:0.45.0-plucky

WORKDIR /opt/smallRNA_seq_pipeline
COPY pixi.toml pixi.lock ./

# Install dependencies and remove cache
# post-link scripts necessary for geneinfodb package
RUN pixi config set --local run-post-link-scripts insecure && \
    pixi install --locked && \
    rm -rf ~/.cache/rattler

COPY scripts/ scripts/

ENTRYPOINT ["pixi", "run"]
CMD ["bash"]

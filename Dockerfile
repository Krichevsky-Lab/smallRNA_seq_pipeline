FROM ghcr.io/prefix-dev/pixi:0.41.4

WORKDIR /opt/smallRNA_seq_pipeline
COPY pixi.toml pixi.lock ./
RUN pixi install --locked

COPY . .

ENTRYPOINT ["pixi", "run"]
CMD ["bash"]

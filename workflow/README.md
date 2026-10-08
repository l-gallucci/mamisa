# MaMISA Snakemake workflow

Runs the MAG quality-control tail of MaMISA with **one conda env per rule**, so
CheckM2, GTDB-Tk and GUNC (which pin conflicting dependencies) never share an
env. Snakemake builds each env automatically on first run.

## Install Snakemake

```bash
conda install -c conda-forge -c bioconda snakemake-minimal
# or: pip install snakemake
```

## Configure

The quick way — let MaMISA do the edits:

```bash
mamisa setup-workflow --genomes-dir bins/ --genome-ext fa --threads 40 --install
mamisa fetch-databases --gtdbtk-data /data/gtdbtk_r220 \
                       --gunc-db /data/gunc/gunc_db_progenomes2.1.dmnd   # or --download
```

`setup-workflow` points each `envs/*.yaml` at this checkout (`pip -> -e <repo>`)
and writes the paths/threads into `config.yaml`; `--install` builds every
per-rule conda env. `fetch-databases` registers existing databases into the
config, or downloads missing ones with each tool's own downloader.

Or edit by hand: set `genomes_dir`, `genome_ext`, the database paths
(`gtdbtk_data`, `gunc_db`), `threads`, `sample_regex` / `tax_level` in
`config.yaml`, and in each `envs/*.yaml` replace `- mamisa` with
`-e /absolute/path/to/this/repo` until MaMISA is on PyPI.

## Run

```bash
# from the repo root
snakemake --use-conda --cores 20 -s workflow/Snakefile --configfile workflow/config.yaml

# preview the plan without running
snakemake -n -s workflow/Snakefile --configfile workflow/config.yaml
```

## DAG

```
genomes_dir ─┬─> run-checkm2 ─┐
             ├─> run-gtdbtk ──┼─> organize-mags (HQ/MQ/LQ + rename)
             └─> run-gunc ────┘
                      └──────────> check-chimeras (GC + contamination + GUNC)
```

Binning and anvi'o misassembly detection are assumed done upstream; this
workflow starts from a directory of MAG FASTAs.

## Outputs (under `results/`)

| Path | Produced by |
|---|---|
| `checkm2/quality_report.tsv` | run-checkm2 |
| `gtdbtk/` | run-gtdbtk |
| `gunc/` | run-gunc |
| `organized/Selected/{HQ,MQ,LQ}/` | organize-mags (renamed genomes) |
| `organized/merged_quality.tsv` | organize-mags |
| `chimera_report.tsv` | check-chimeras |

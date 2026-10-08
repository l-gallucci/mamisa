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

Edit `workflow/config.yaml`: set `genomes_dir`, `genome_ext`, the database paths
(`gtdbtk_data`, `gunc_db`), `threads`, and the `sample_regex` / `tax_level` used
when `organize-mags` renames outputs.

> Until MaMISA is on PyPI, edit each `workflow/envs/*.yaml` and replace the
> `pip: - mamisa` line with `-e /absolute/path/to/this/repo`.

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

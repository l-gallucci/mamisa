# Running on an HPC cluster (SLURM)

Create the conda envs once on a shared filesystem (they are reused by every job).
Two ways to submit.

## 1. Snakemake submits each rule as its own SLURM job (recommended)

Each tool gets its own CPUs/RAM/time:

```bash
pip install snakemake-executor-plugin-slurm
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
    --executor slurm --jobs 20 \
    --default-resources slurm_partition=standard runtime=1440 \
    --set-resources run_gtdbtk:mem_mb=320000 run_gunc:mem_mb=64000 \
                    run_checkm2:mem_mb=64000
```

SLURM settings (partition, time, memory, account, mail) come from
`--default-resources` / `--set-resources`, or from a profile YAML. Example
profile (`~/.config/snakemake/slurm/config.yaml`):

```yaml
executor: slurm
jobs: 20
default-resources:
  slurm_partition: standard
  runtime: 1440            # minutes
  slurm_account: my_account
set-resources:
  run_gtdbtk:
    mem_mb: 320000
    cpus_per_task: 40
  run_gunc:
    mem_mb: 64000
  run_checkm2:
    mem_mb: 64000
# mail via slurm_extra, e.g.:
# default-resources:
#   slurm_extra: "'--mail-type=ALL --mail-user=you@example.org'"
```

Then just `snakemake --profile slurm -s workflow/Snakefile --configfile workflow/config.yaml`.

## 2. A single sbatch script (fine for one sample / modest data)

Runs Option A's `conda run` commands in order inside one job:

```bash
#!/bin/bash
#SBATCH -J mamisa -c 40 --mem=320G -t 24:00:00
source ~/miniconda3/etc/profile.d/conda.sh
conda run -n checkm2 mamisa run-checkm2 --genome-dir bins/ -o results/checkm2 --threads 40
conda run -n gtdbtk  env GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH mamisa run-gtdbtk --genome-dir bins/ -o results/gtdbtk --cpus 40
conda run -n gunc    mamisa run-gunc    --genome-dir bins/ -o results/gunc --file-suffix .fa --threads 40 --db-file $GUNC_DB
# ... organize-mags, check-chimeras in the mamisa env
```

## Notes

- GTDB-Tk needs a lot of RAM (pplacer: ~150-320 GB depending on release; r232
  requires at least 140 GB). GUNC and CheckM2 are lighter (~16-64 GB).
- `run-gtdbtk --skip-existing` reuses an existing `gtdbtk.*.summary.tsv` instead
  of re-classifying. Useful for the sbatch form; inside Snakemake it has limited
  effect because Snakemake deletes a rule's output dir before re-running it.
- If your cluster provides tools as `module load` instead of conda, drop the
  `conda run -n <env>` prefix and `module load <tool>` before each step (Option A
  / sbatch). Snakemake named-env mode assumes conda, so prefer the sbatch form
  with modules.

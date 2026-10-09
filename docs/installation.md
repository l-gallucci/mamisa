# Installation

MaMISA itself is pure Python with no external dependencies (Python >= 3.9), so it
installs anywhere in seconds:

```bash
git clone https://github.com/l-gallucci/mamisa.git
cd mamisa
pip install -e .

# Verify
mamisa --version
mamisa --help
```

That is enough to run the dependency-free commands (`organize-mags`,
`check-chimeras` without GUNC, and the helpers). The heavy commands additionally
call external tools (see below).

## How the environments fit together

The key fact: MaMISA is tiny and has no dependencies, but the tools it drives
(CheckM2, GTDB-Tk, GUNC, anvi'o) pin conflicting versions and cannot live in one
conda env. So the layout is:

- one conda env per heavy tool (unavoidable, standard in metagenomics), and
- MaMISA installed into each env you run from (a 1-second `pip install`, since it
  has no dependencies). The wrapper `mamisa run-gunc` must run inside the env that
  has `gunc`, so `mamisa` has to be importable there too.

| Conda env | Contains | MaMISA commands that run in it |
|---|---|---|
| `mamisa` (light) | MaMISA + samtools, BLAST+, bedtools, Kraken2, Kaiju, Meryl | `organize-mags`, `check-chimeras`, `classify-clipping`, `check-read-chimeras`, `check-zero-coverage`, `process-large-contigs`, and the `setup-workflow` / `fetch-databases` / `check-envs` helpers |
| `checkm2` | CheckM2 (+ MaMISA) | `run-checkm2` |
| `gtdbtk` | GTDB-Tk (+ MaMISA) | `run-gtdbtk` |
| `gunc` | GUNC + DIAMOND + Prodigal (+ MaMISA) | `run-gunc` |
| `mimag` | barrnap + tRNAscan-SE (+ MaMISA) | `run-mimag-rna` |
| `anvio-9` | anvi'o (upstream, no MaMISA needed) | produces the `*-clipping.txt` inputs |

Do I need all of them? Only the ones whose commands you use. Just want
`organize-mags` / `check-chimeras`? The `mamisa` env alone is enough. Want the full
assembly-to-taxonomy pipeline? Create all of them.

Snakemake is optional: it is only an orchestrator on top of these envs (see
[pipeline.md](pipeline.md)). You never put MaMISA in a Snakemake env; Snakemake
just calls the envs above.

## External tools (required per command)

MaMISA's wrappers shell out to standard bioinformatics tools. A tool only needs
to be on `PATH` for the command that uses it.

| Command | External tool(s) |
|---|---|
| `process-large-contigs` | CheckM2 |
| `run-checkm2` | CheckM2 |
| `run-gtdbtk` | GTDB-Tk (>=2.5 uses skani; Mash/`--mash-db` no longer used) |
| `run-gunc` | GUNC + DIAMOND + Prodigal |
| `run-mimag-rna` | barrnap + tRNAscan-SE |
| `check-read-chimeras` | samtools, Kraken2 or Kaiju |
| `check-chimeras` | (none; `--gunc-dir` consumes GUNC output) |
| `classify-clipping` | samtools (+ BLAST+ for `--self-blast`) |
| `check-zero-coverage` | BLAST+ (`blastn`/`makeblastdb`), optional Meryl |
| `filter-misassemblies` | anvi'o (upstream misassembly detection) |
| `organize-mags` (was `filter-checkm2`) | (none; consumes CheckM2 + GTDB-Tk output) |

## Create the environments

Create only the ones you need (see the table above). `mamba` or `conda` both work.

```bash
# 1) Light env: MaMISA + the small CLI tools. Needed for most commands.
mamba create -n mamisa -c conda-forge -c bioconda \
    python=3.10 samtools blast bedtools kraken2 kaiju meryl
conda activate mamisa
pip install -e .                     # installs MaMISA (clone this repo first)

# 2) CheckM2
mamba create -n checkm2 -c conda-forge -c bioconda checkm2
conda run -n checkm2 pip install -e .    # add MaMISA (zero-dep, instant)

# 3) GTDB-Tk
mamba create -n gtdbtk -c conda-forge -c bioconda gtdbtk
conda run -n gtdbtk pip install -e .

# 4) GUNC - OLD tool: PIN the versions or its recipe pulls pandas 3.x +
#    a too-new Python that crash it at runtime.
mamba create -n gunc -c conda-forge -c bioconda \
    "gunc=1.1.1" "python=3.10" "pandas>=2,<3" "numpy<2" diamond prodigal
conda run -n gunc pip install -e .

# 5) MIMAG rRNA/tRNA - small env for run-mimag-rna (optional, full-MIMAG HQ)
mamba create -n mimag -c conda-forge -c bioconda barrnap trnascan-se
conda run -n mimag pip install -e .

# 6) anvi'o - upstream only, produces *-clipping.txt (no MaMISA needed)
mamba create -n anvio-9 -c conda-forge -c bioconda anvio=9
```

Check everything is compatible at any time:

```bash
mamisa check-envs --mamisa-env mamisa --checkm2-env checkm2 \
                  --gtdbtk-env gtdbtk --gunc-env gunc --mimag-env mimag
```

Verified (2026-10): GUNC 1.1.1's Bioconda recipe declares `pandas>=2.0.0` with no
upper bound, so an unpinned install grabs pandas 3.x and Python 3.14+, which crash
GUNC. The pinned recipe above resolves to python 3.10, pandas 2.3, numpy 1.26,
diamond 2.1.24, prodigal 2.6.3 and runs cleanly. `mamisa run-gunc` also runs an
env-health check and refuses to launch a broken GUNC.

## Databases

The path KIND is checked per tool: GTDB-Tk = a directory, GUNC and CheckM2 = a
`.dmnd` file.

```bash
# already have them? just point MaMISA at them (validated, written to config)
mamisa fetch-databases \
    --gtdbtk-data /data/gtdbtk_r220 \                      # DIRECTORY
    --gunc-db     /data/gunc/gunc_db_progenomes2.1.dmnd \  # FILE
    --checkm2-db  /data/checkm2/uniref100.KO.1.dmnd        # FILE

# or download the missing ones (each runs inside its tool's env)
mamisa fetch-databases --download --db-dir /data/mamisa_dbs
```

## Running the tests

```bash
pip install -e ".[dev]"   # installs pytest
pytest                    # 96 unit tests, no external tools needed
```

The suite covers the dependency-free logic: contig-id parsing, CheckM2 tiering,
GC/GUNC chimera scoring, the SA-tag join detector, clipping classification,
organize-mags renaming, MIMAG rRNA/tRNA parsing, per-tool database-path
validation, and the config/env helpers.

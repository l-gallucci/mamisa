# MaMISA — Manage Misassemblies

A comprehensive toolkit for metagenomic assembly quality control and filtering.

## Overview

MaMISA provides a set of commands that cover the full workflow from raw assembly to
taxonomically classified, quality-filtered genomes.

### Commands at a glance

**Assembly QC & misassembly handling**

| Command | What it does | Example |
|---|---|---|
| `process-large-contigs` | Extract, QC and filter large contigs before binning | `mamisa process-large-contigs -a assembly.fa -o large/` |
| `classify-clipping` | Label each clipping position (repeat/deletion/**chimeric_join**/…) from BAM evidence | `mamisa classify-clipping --bam map.bam -m mis/ -o clip.tsv` |
| `filter-misassemblies` | Split or remove contigs at misassembly positions | `mamisa filter-misassemblies -a assembly.fa -m mis/ -o clean.fa` |
| `remove-hq-contigs` | Remove known HQ-genome contigs from an assembly | `mamisa remove-hq-contigs -a assembly.fa --hq-dir hq/ -o out.fa` |
| `check-zero-coverage` | Validate no-coverage regions via BLAST + k-mers | `mamisa check-zero-coverage --assembly assembly.fa ...` |

**Chimera detection**

| Command | What it does | Example |
|---|---|---|
| `check-read-chimeras` | Chimeras from read-level taxonomy (Kraken2 + BAM) | `mamisa check-read-chimeras --bam map.bam --kraken2-output k2.txt -o chim/` |
| `check-chimeras` | Chimeric MAGs from GC + contamination + **GUNC** + GTDB signals | `mamisa check-chimeras --bins-dir bins/ --gunc-dir gunc/ -o report.tsv` |
| `run-gunc` | Gene-level chimerism/contamination (GUNC wrapper) | `mamisa run-gunc --genome-dir bins/ -o gunc/ --db-file gunc.dmnd` |

**MAG quality & taxonomy**

| Command | What it does | Example |
|---|---|---|
| `run-checkm2` | Completeness/contamination (CheckM2 wrapper) | `mamisa run-checkm2 --genome-dir bins/ -o checkm2/ --threads 40` |
| `run-gtdbtk` | Taxonomy classification (GTDB-Tk wrapper) | `mamisa run-gtdbtk --genome-dir bins/ -o gtdbtk/ --cpus 40` |
| `run-mimag-rna` | MIMAG rRNA (5S/16S/23S) + tRNA criterion (barrnap + tRNAscan-SE) | `mamisa run-mimag-rna --genome-dir bins/ -o mimag_rna/ --threads 8` |
| `organize-mags` | Split genomes into HQ/MQ/LQ, rename `<sample>__<taxon>__<orig>`, optional full-MIMAG HQ gate | `mamisa organize-mags --checkm2-root checkm2/ --genomes-dir bins/ --gtdbtk-dir gtdbtk/ -o filtered/ --rename --copy` |

**Setup & environment helpers**

| Command | What it does | Example |
|---|---|---|
| `setup-workflow` | Configure the Snakemake workflow (env yamls + config), optionally build envs | `mamisa setup-workflow --genomes-dir bins/ --gunc-env gunc --checkm2-env checkm2` |
| `fetch-databases` | Register existing databases, or download missing ones | `mamisa fetch-databases --gtdbtk-data /db/gtdbtk --gunc-db /db/gunc.dmnd` |
| `check-envs` | Report tool versions per conda env + MaMISA compatibility | `mamisa check-envs --gunc-env gunc --checkm2-env checkm2` |

See [Step-by-Step Guide](#step-by-step-guide) for the full pipeline and
[Command Reference](#command-reference) for every flag.

---

## Installation

MaMISA itself is **pure Python with no external dependencies** (Python ≥ 3.9), so
it installs anywhere in seconds:

```bash
git clone https://github.com/l-gallucci/mamisa.git
cd mamisa
pip install -e .

# Verify
mamisa --version
mamisa --help
```

That is enough to run the **dependency-free commands** (`organize-mags`,
`check-chimeras` without GUNC, and the helpers). The heavy commands additionally
call external tools — read on.

### How the environments fit together (read this first)

The key fact: **MaMISA is tiny and has no dependencies**, but the tools it drives
(CheckM2, GTDB-Tk, GUNC, anvi'o) pin **conflicting** versions and **cannot live in
one conda env**. So the layout is:

- **one conda env per heavy tool** (unavoidable, standard in metagenomics), and
- **MaMISA installed into each env you run from** — a 1-second `pip install`, since
  it has no dependencies. The wrapper `mamisa run-gunc` must run *inside* the env
  that has `gunc`, so `mamisa` has to be importable there too.

| Conda env | Contains | MaMISA commands that run in it |
|---|---|---|
| `mamisa` (light) | MaMISA + samtools, BLAST+, bedtools, Kraken2, Meryl | `organize-mags`, `check-chimeras`, `classify-clipping`, `check-read-chimeras`, `check-zero-coverage`, `process-large-contigs`, and the `setup-workflow` / `fetch-databases` / `check-envs` helpers |
| `checkm2` | CheckM2 (+ MaMISA) | `run-checkm2` |
| `gtdbtk` | GTDB-Tk (+ MaMISA) | `run-gtdbtk` |
| `gunc` | GUNC + DIAMOND + Prodigal (+ MaMISA) | `run-gunc` |
| `mimag` | barrnap + tRNAscan-SE (+ MaMISA) | `run-mimag-rna` |
| `anvio-9` | anvi'o (upstream, no MaMISA needed) | produces the `*-clipping.txt` inputs |

**Do I need all of them?** Only the ones whose commands you use. Just want
`organize-mags` / `check-chimeras`? The `mamisa` env alone is enough. Want the full
assembly→taxonomy pipeline? Create all of them.

**Snakemake is optional** — it is only an orchestrator *on top* of these envs
(see [Run the pipeline](#run-the-pipeline)). You never put MaMISA "in a Snakemake
env"; Snakemake just calls the envs above.

### Running the tests

```bash
pip install -e ".[dev]"   # installs pytest
pytest                    # 96 unit tests, no external tools needed
```

The suite covers the dependency-free logic: contig-id parsing, CheckM2 tiering,
GC/GUNC chimera scoring, the SA-tag join detector, clipping classification,
organize-mags renaming, per-tool database-path validation, and the config/env
helpers.

### External tools (required per command)

MaMISA's wrappers shell out to standard bioinformatics tools. A tool only needs
to be on `PATH` for the command that uses it.

| Command | External tool(s) |
|---|---|
| `process-large-contigs` | CheckM2 |
| `run-checkm2` | CheckM2 |
| `run-gtdbtk` | GTDB-Tk (>=2.5 uses skani; Mash/`--mash-db` no longer used) |
| `run-gunc` | GUNC + DIAMOND + Prodigal |
| `run-mimag-rna` | barrnap + tRNAscan-SE |
| `check-read-chimeras` | samtools, Kraken2 |
| `check-chimeras` | (none; `--gunc-dir` consumes GUNC output) |
| `classify-clipping` | samtools (+ BLAST+ for `--self-blast`) |
| `check-zero-coverage` | BLAST+ (`blastn`/`makeblastdb`), optional Meryl |
| `filter-misassemblies` | anvi'o (upstream misassembly detection) |
| `organize-mags` (was `filter-checkm2`) | (none; consumes CheckM2 + GTDB-Tk output) |

### Create the environments

Create only the ones you need (see the table above). `mamba` or `conda` both work.

```bash
# 1) Light env: MaMISA + the small CLI tools. Needed for most commands.
mamba create -n mamisa -c conda-forge -c bioconda \
    python=3.10 samtools blast bedtools kraken2 meryl
conda activate mamisa
pip install -e .                     # installs MaMISA (clone this repo first)

# 2) CheckM2
mamba create -n checkm2 -c conda-forge -c bioconda checkm2
conda run -n checkm2 pip install -e .    # add MaMISA (zero-dep, instant)

# 3) GTDB-Tk
mamba create -n gtdbtk -c conda-forge -c bioconda gtdbtk
conda run -n gtdbtk pip install -e .

# 4) GUNC — OLD tool: PIN the versions or its recipe pulls pandas 3.x +
#    a too-new Python that crash it at runtime.
mamba create -n gunc -c conda-forge -c bioconda \
    "gunc=1.1.1" "python=3.10" "pandas>=2,<3" "numpy<2" diamond prodigal
conda run -n gunc pip install -e .

# 5) MIMAG rRNA/tRNA — small env for run-mimag-rna (optional, full-MIMAG HQ)
mamba create -n mimag -c conda-forge -c bioconda barrnap trnascan-se
conda run -n mimag pip install -e .

# 6) anvi'o — upstream only, produces *-clipping.txt (no MaMISA needed)
mamba create -n anvio-9 -c conda-forge -c bioconda anvio=9
```

Check everything is compatible at any time:

```bash
mamisa check-envs --mamisa-env mamisa --checkm2-env checkm2 \
                  --gtdbtk-env gtdbtk --gunc-env gunc --mimag-env mimag
```

> **Verified (2026-10):** GUNC 1.1.1's Bioconda recipe declares `pandas>=2.0.0`
> with no upper bound, so an unpinned install grabs pandas 3.x and Python 3.14+,
> which crash GUNC. The pinned recipe above resolves to python 3.10, pandas 2.3,
> numpy 1.26, diamond 2.1.24, prodigal 2.6.3 and runs cleanly. `mamisa run-gunc`
> also runs an env-health check and refuses to launch a broken GUNC.

### Databases

The path KIND is checked per tool: **GTDB-Tk = a directory**, **GUNC and CheckM2 =
a `.dmnd` file**.

```bash
# already have them? just point MaMISA at them (validated, written to config)
mamisa fetch-databases \
    --gtdbtk-data /data/gtdbtk_r220 \                      # DIRECTORY
    --gunc-db     /data/gunc/gunc_db_progenomes2.1.dmnd \  # FILE
    --checkm2-db  /data/checkm2/uniref100.KO.1.dmnd        # FILE

# or download the missing ones (each runs inside its tool's env)
mamisa fetch-databases --download --db-dir /data/mamisa_dbs
```

### Run the pipeline

Two equivalent ways. Both use the same envs created above.

**Option A — run commands yourself (simplest).** Call each command in its env:

```bash
conda run -n checkm2 mamisa run-checkm2 --genome-dir bins/ -o results/checkm2 --threads 40
conda run -n gtdbtk  env GTDBTK_DATA_PATH=/data/gtdbtk_r220 \
                     mamisa run-gtdbtk  --genome-dir bins/ -o results/gtdbtk --cpus 40
conda run -n gunc    mamisa run-gunc    --genome-dir bins/ -o results/gunc \
                     --file-suffix .fa --threads 40 --db-file /data/gunc/gunc_db_progenomes2.1.dmnd
conda run -n mimag   mamisa run-mimag-rna --genome-dir bins/ -o results/mimag_rna --threads 40  # optional
conda run -n mamisa  mamisa organize-mags --checkm2-root results/checkm2 --genomes-dir bins/ \
                     --gtdbtk-dir results/gtdbtk -o results/filtered --rename --copy \
                     --mimag-rna-dir results/mimag_rna          # drop this flag to skip the rRNA/tRNA HQ gate
conda run -n mamisa  mamisa check-chimeras --bins-dir bins/ --gunc-dir results/gunc \
                     --checkm2-report results/checkm2/quality_report.tsv -o results/chimera_report.tsv
```

> MIMAG HQ also requires 5S/16S/23S rRNA + tRNAs for ≥18 amino acids, which
> CheckM2 does not check. `run-mimag-rna` adds that evidence and
> `organize-mags --mimag-rna-dir` demotes HQ genomes that fail it to MQ. Skip it
> and HQ is completeness/contamination only.

**Option B — Snakemake orchestrates it.** One command runs the whole DAG. First
configure once:

```bash
conda install -c conda-forge -c bioconda snakemake-minimal   # in base or its own env

# tell the workflow which envs + databases to use (sanity-checks them)
mamisa setup-workflow --genomes-dir bins/ --genome-ext fa --threads 40 \
    --mamisa-env mamisa --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc

# then run (named-env mode: Snakemake wraps each rule in `conda run -n`)
snakemake --cores 40 -s workflow/Snakefile --configfile workflow/config.yaml
```

Don't have the envs and want Snakemake to **build them for you** instead? Use
`mamisa setup-workflow --genomes-dir bins/ --install` and run with
`snakemake --use-conda ...`. Details in [`workflow/README.md`](workflow/README.md).

### Running on an HPC cluster (SLURM)

Create the conda envs once on a shared filesystem (they are reused by every job).
Two ways to submit:

**1. One Snakemake run that submits each rule as its own SLURM job** (recommended —
each tool gets its own CPUs/RAM/time):

```bash
pip install snakemake-executor-plugin-slurm
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
    --executor slurm --jobs 20 \
    --default-resources slurm_partition=standard runtime=1440 \
    --set-resources run_gtdbtk:mem_mb=320000 run_gunc:mem_mb=64000 \
                    run_checkm2:mem_mb=64000
```

**2. A single `sbatch` script** that runs Option A's `conda run` commands in order
(fine for one sample / modest data):

```bash
#!/bin/bash
#SBATCH -J mamisa -c 40 --mem=320G -t 24:00:00
source ~/miniconda3/etc/profile.d/conda.sh
conda run -n checkm2 mamisa run-checkm2 --genome-dir bins/ -o results/checkm2 --threads 40
conda run -n gtdbtk  env GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH mamisa run-gtdbtk --genome-dir bins/ -o results/gtdbtk --cpus 40
conda run -n gunc    mamisa run-gunc    --genome-dir bins/ -o results/gunc --file-suffix .fa --threads 40 --db-file $GUNC_DB
# ... organize-mags, check-chimeras in the mamisa env
```

Notes: **GTDB-Tk needs a lot of RAM** (pplacer: ~150–320 GB depending on release);
GUNC and CheckM2 are lighter (~16–64 GB). If your cluster provides the tools as
`module load` instead of conda, drop the `conda run -n <env>` prefix and `module
load <tool>` before each step (Option A / sbatch); Snakemake named-env mode assumes
conda, so prefer the sbatch form with modules.

---

## Complete Workflow

```
assembly.fa
    │
    ▼
[1] process-large-contigs       Separate very large contigs, run CheckM2 on each,
    │                           extract clean HQ genomes, return updated assembly.
    │
    ▼
[2] anvi-script-find-misassemblies   (external — anvi'o)
    │                           Detect soft-clipping positions in the BAM.
    │                           Produces  *-clipping.txt  files.
    │
    ▼
[3] check-read-chimeras         (optional, recommended)
    │                           BAM + Kraken2 per-read taxonomy.
    │                           Flags contigs whose reads come from ≥2 organisms.
    │                           Outputs  chimera_read_report.tsv
    │                                    chimera_read_windows.tsv
    │
    ├──────────────────────────►[3b] check-chimeras   (optional)
    │                                GC-based chimera detection on bins
    │                                (run after binning if preferred).
    │
    ▼
[4] classify-clipping           (optional, recommended)
    │                           Uses BAM evidence to label each clipping position:
    │                             end_artefact / repeat_collapse / deletion_artefact /
    │                             chimera_candidate / sv_candidate / low_confidence
    │                           Cross-references chimera_read_windows.tsv when provided.
    │
    ▼
[5] filter-misassemblies        Split or remove misassembled contigs.
    │                           Preserves HQ genomes. Integrates chimera report.
    │                           Handles HQ circular contigs (--split-hq-circular).
    │
    ▼
[6] binning                     (external — MetaBAT2, MaxBin2, etc.)
    │
    ▼
[7] checkm2 predict             (external — CheckM2)
    │
    ▼
[8] filter-checkm2              Organise bins into HQ / MQ / LQ quality tiers.
    │
    ▼
[9] run-gtdbtk                  Taxonomic classification of selected genomes.
```

---

## Step-by-Step Guide

### Step 1 — Process large contigs

Very long contigs (> 300 kbp by default) often represent single complete genomes and
should be handled separately before binning. This command:

1. Separates large and regular contigs
2. Runs CheckM2 on each large contig individually
3. Classifies each one as: *extract HQ*, *keep for splitting*, or *low quality*
4. Returns an updated assembly ready for misassembly detection

```bash
mamisa process-large-contigs \
    --assembly assembly.fa \
    --misassemblies misasm_dir/ \
    --output-dir 01_large_contigs/ \
    --max-length 300000 \
    --min-completeness 50 \
    --max-contamination 10 \
    --threads 40
```

Output:
```
01_large_contigs/
├── 01_extracted/
│   ├── large_contigs.fa          # contigs > max-length
│   └── assembly_regular.fa       # contigs ≤ max-length
├── 02_individual/                 # one .fa per large contig (CheckM2 input)
├── 03_checkm2/                    # CheckM2 results
│   └── quality_report.tsv
├── HQ_extracted/                  # clean HQ genomes extracted here
├── filtering_decisions.tsv
└── assembly_for_filtering.fa      # use this in Step 5
```

---

### Step 2 — Detect misassemblies with anvi'o (external)

```bash
anvi-script-find-misassemblies \
    -b mapping.bam \
    -o misassemblies/MisAsm \
    -T 40
```

Produces `MisAsm-clipping.txt` in `misassemblies/`.

---

### Step 3 — Detect chimeric contigs (read-level taxonomy)

This command streams the BAM once, cross-references every read against Kraken2
per-read taxonomy output, and builds a per-contig taxonomic profile. Contigs whose
reads originate from more than one organism are flagged as chimeric.

```bash
mamisa check-read-chimeras \
    --bam mapping.bam \
    --kraken2-output kraken2_reads.txt \
    --kraken2-report kraken2_report.txt \
    -o 02_chimera/ \
    --min-mapq 20 \
    --window 10000 \
    --window-step 5000
```

Outputs:
```
02_chimera/
├── chimera_read_report.tsv    # per-contig: risk level, dominant taxon, diversity
└── chimera_read_windows.tsv   # sliding-window detail (for classify-clipping)
```

Chimera risk levels: **High** / **Medium** / **Low** / **Clean** / **Insufficient**

Key scoring signals:
- Dominant-taxon fraction < 80 %: elevated score
- > 1 distinct taxon detected: elevated score
- > 0 windows with taxon shift: highest score

---

### Step 3b — Detect chimeric MAGs (GC-based, optional)

Run after binning when you want a GC-composition check on complete bins
rather than on individual contigs.

```bash
mamisa check-chimeras \
    --bins-dir 06_bins/ \
    --output 02_chimera/gc_chimera_report.tsv \
    --gtdbtk-dir 09_taxonomy/ \
    --checkm2-report 07_checkm2/quality_report.tsv \
    --gc-window 5000 \
    --gc-step 2500
```

---

### Step 4 — Classify clipping positions

Before splitting, characterise each anvi'o-reported clipping position using
BAM evidence so you can prioritise which splits are biologically meaningful.

```bash
mamisa classify-clipping \
    --bam mapping.bam \
    --misassemblies misassemblies/ \
    --taxonomy-windows 02_chimera/chimera_read_windows.tsv \
    --output 03_classified/clipping_classified.tsv
```

Output columns:
```
contig  clip_pos  contig_length  local_depth  primary_reads_in_window
contig_mean_depth  depth_ratio  discordant_fraction  large_insert_fraction
strand_fwd_fraction  clipped_base_entropy  near_contig_end  taxonomy_shift
classification  confidence  evidence
```

Classification labels:

| Label | Biological meaning |
|---|---|
| `end_artefact` | Near contig terminus; assembly edge noise |
| `repeat_collapse` | Coverage spike + low-entropy clipped bases; collapsed repeat |
| `deletion_artefact` | Local depth drop; internal deletion or coverage collapse |
| `chimera_candidate` | Discordant pairs + large inserts ± taxonomy shift |
| `sv_candidate` | SV signal (discordant) without depth anomaly or taxon shift |
| `low_confidence` | Ambiguous or insufficient evidence |

---

### Step 5 — Filter misassemblies

Processes the assembly from Step 1 using the clipping data from Step 2.
Optionally integrates the chimera report from Step 3.

```bash
# Conservative: use all anvi'o-reported positions (--safe, the default)
mamisa filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --hq-genomes 01_large_contigs/HQ_extracted/ \
    --output 04_clean/assembly_clean.fa \
    --mode split \
    --preserve-hq-with-issues \
    --chimera-report 02_chimera/chimera_read_report.tsv \
    --split-hq-circular \
    --safe

# Selective: only split at positions with ≥50 clipped reads
mamisa filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --output 04_clean/assembly_clean.fa \
    --min-clip-coverage 50
```

**Clipping position selection:**

| Flag | Behaviour |
|---|---|
| `--safe` (default) | Use ALL positions reported by anvi'o (trust the tool's threshold) |
| `--min-clip-coverage N` | Only split where ≥ N reads are clipped (selective, less aggressive) |

**Contig classification and actions:**

| Contig type | Default action | With `--preserve-hq-with-issues` |
|---|---|---|
| HQ + no clipping | Removed (already extracted) | Removed |
| HQ + clipping zone | Removed | **Split** |
| HQ + circular + clipping | Removed | **Split** (requires `--split-hq-circular`) |
| HQ + chimera flag only | Removed | **Split** |
| Non-HQ + clipping | Split (mode=split) or removed | Same |
| Non-HQ + clean | Kept intact | Kept intact |

---

### Step 6 — Binning (external)

```bash
jgi_summarize_bam_contig_depths --outputDepth depth.txt mapping.bam

metabat2 \
    -i 04_clean/assembly_clean.fa \
    -a depth.txt \
    -o 05_bins/bin \
    -t 40
```

---

### Step 7 — Quality assessment with CheckM2 (external)

```bash
checkm2 predict \
    --threads 40 \
    --input 05_bins/ \
    --output-directory 06_checkm2/ \
    -x fa
```

---

### Step 8 — Filter genomes by quality

```bash
mamisa filter-checkm2 \
    --checkm2-root 06_checkm2/ \
    --genomes-dir 05_bins/ \
    --output 07_filtered/ \
    --tiers HQ,MQ \
    --symlink
```

**MiMAG quality thresholds (defaults):**

| Tier | Completeness | Contamination |
|---|---|---|
| HQ | ≥ 90% | ≤ 5% |
| MQ | ≥ 70% | ≤ 10% |
| LQ | ≥ 50% | ≤ 10% |

---

### Step 9 — Taxonomic classification with GTDB-Tk

```bash
mamisa run-gtdbtk \
    --selected-dir 07_filtered/Selected/ \
    --output 08_taxonomy/ \
    --cpus 40 \
    --tiers HQ,MQ
```

---

## Command Reference

### process-large-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output-dir PATH       Output directory

Thresholds:
  --max-length INT            Large contig threshold (default: 300000)
  --min-completeness FLOAT    Min completeness for HQ (default: 50)
  --max-contamination FLOAT   Max contamination for HQ (default: 10)

CheckM2:
  --threads INT               Threads (default: 1)
  --skip-checkm2              Skip CheckM2 (requires --checkm2-results)
  --checkm2-results PATH      Existing CheckM2 output directory
```

### check-read-chimeras

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  --kraken2-output PATH       Kraken2 per-read classification output
  -o, --output-dir PATH       Output directory

Optional:
  --kraken2-report PATH       Kraken2 report (adds taxon names to output)
  --assembly PATH             Assembly FASTA (for contig lengths if not in BAM)
  --min-mapq INT              Min mapping quality (default: 20)
  --min-reads INT             Min reads per contig to report (default: 10)
  --window INT                Sliding window size in bp (default: 10000)
  --window-step INT           Window step in bp (default: 5000)
  --window-threshold INT      Min contig length for windowed analysis (default: 50000)
  --exclude-unclassified      Exclude taxid=0 reads from diversity calculations
  --dry-run
```

### check-chimeras

```
Required:
  --bins-dir PATH             Directory containing bin FASTA files

Optional:
  -o, --output PATH           Output TSV (default: chimera_report.tsv)
  --gtdbtk-dir PATH           GTDB-Tk output directory (adds taxonomy signals)
  --checkm2-report PATH       CheckM2 report (adds contamination signal)
  --gunc-dir PATH             GUNC output directory (adds gene-level clade-consistency signal)
  --gc-window INT             GC window size in bp (default: 5000)
  --gc-step INT               GC step in bp (default: 2500)
  --taxonomy-level STR        Taxonomy level for comparison (default: phylum)
  --extensions LIST           Genome extensions (default: fa,fasta,fna)
  --dry-run
```

### classify-clipping

```
Required:
  --bam PATH                  Sorted, indexed BAM file
  -m, --misassemblies PATH    Directory with *-clipping.txt files
  -o, --output PATH           Output classification TSV

Optional:
  --min-mapq INT              Min mapping quality (default: 20)
  --window INT                Read window around each position in bp (default: 500)
  --min-clip-coverage N       Only classify positions with ≥N clipped reads
  --taxonomy-windows TSV      chimera_read_windows.tsv from check-read-chimeras
  --taxonomy-flank BP         Extend shift windows by this many bp (default: 5000)
  --assembly FASTA            Assembly FASTA — required for --self-blast
  --self-blast                Self-BLAST repeat detection (needs --assembly + BLAST+)
  --dry-run
```

Classification labels include `chimeric_join`: when the soft-clipped parts of
reads at a position map (via their `SA:Z:` supplementary-alignment tag)
consistently to **one other contig**, the two contigs are joined there. This is
the read-level join signal from Trigodet et al. 2025 and needs no extra flag —
it is computed automatically from the BAM. See columns `sa_partner_contig` /
`sa_partner_fraction` in the output.

### filter-misassemblies

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  -m, --misassemblies PATH    Directory with *-clipping.txt files

Optional:
  -g, --hq-genomes PATH       Directory with HQ genome files
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum fragment length (default: 2500)
  --mode {remove,split}       Misassembly handling (default: split)
  --preserve-hq-with-issues   Split HQ contigs with clipping instead of removing

Clipping selection (mutually exclusive):
  --safe                      Use ALL anvi'o-reported positions (default)
  --min-clip-coverage N       Only split where ≥N reads are clipped

Chimera awareness:
  --chimera-report PATH       TSV from check-chimeras or check-read-chimeras
  --chimera-risk-threshold    Min risk level to act on (default: Medium)

HQ circular handling:
  --split-hq-circular         Force-split HQ circular contigs with clipping
  --hq-circular-min-length    Length threshold for circular detection (default: 200000)

Other:
  --dry-run
  --stats PATH                Save statistics to TSV
```

### remove-hq-contigs

```
Required:
  -a, --assembly PATH         Input assembly FASTA
  --hq-dir PATH               Directory with HQ genome files
    OR
  --hq-list PATH              Text file with HQ contig IDs (one per line)

Optional:
  -o, --output PATH           Output filtered assembly
  -l, --min-length INT        Minimum contig length (default: 0)
  --dry-run
  --stats PATH                Save statistics to TSV
```

### organize-mags  (formerly `filter-checkm2`)

Splits genomes into `Selected/HQ`, `Selected/MQ`, `Selected/LQ` by MiMAG tier,
and can rename each output `<sample>__<taxon>__<original>` using the sample
parsed from the filename and the GTDB-Tk taxonomy. `filter-checkm2` still works
as a deprecated alias.

```
Required:
  --checkm2-root PATH         Root directory with CheckM2 results
  --genomes-dir PATH          Directory with genome files
  -o, --output PATH           Output directory

Quality thresholds:
  --hq-comp-min FLOAT         HQ min completeness (default: 90)
  --hq-cont-max FLOAT         HQ max contamination (default: 5)
  --mq-comp-min FLOAT         MQ min completeness (default: 70)
  --mq-cont-max FLOAT         MQ max contamination (default: 10)
  --lq-comp-min FLOAT         LQ min completeness (default: 50)
  --lq-cont-max FLOAT         LQ max contamination (default: 10)

Renaming (sample + taxonomy):
  --rename                    Rename outputs as <sample>__<taxon>__<original>
  --gtdbtk-dir PATH           GTDB-Tk summaries, for the <taxon> part of the name
  --sample-regex RE           Capture group for the sample id
                              (default: ^([^._]+) = token before first . or _)
  --tax-level LEVEL           Rank used in name: domain..species (default: genus;
                              falls back up ranks, NoTax if no GTDB hit)

MIMAG rRNA/tRNA HQ gate:
  --mimag-rna-dir PATH        run-mimag-rna output (mimag_rna_summary.tsv);
                              demotes HQ genomes lacking 5S/16S/23S rRNA or
                              >=18 tRNA amino acids to MQ

Other:
  --tiers LIST                Tiers to select (default: HQ,MQ,LQ)
  --extensions LIST           Genome extensions (default: fa,fasta,fna,...)
  --symlink | --copy          Link or copy files (default: symlink)
  --name-map / --strip-prefix / --strip-suffix / --name-prefix / --name-suffix
  --dry-run
```

> **MiMAG note:** the defaults above split the MIMAG *medium* band into MQ/LQ.
> For strict MIMAG tiers (HQ comp>90/cont<5, MQ comp≥50/cont<10, LQ comp<50) run
> with `--mq-comp-min 50 --mq-cont-max 10 --lq-comp-min 0 --lq-cont-max 10`.
> Full MIMAG HQ additionally needs rRNA/tRNA — add `--mimag-rna-dir` (above).

### run-gtdbtk

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --cpus INT                  CPUs for GTDB-Tk (default: 1)
  --place-species             Pass --place_species (GTDB-Tk >=2.7): place in the
                              pplacer tree even when skani classifies
  --mash-db PATH              DEPRECATED/ignored — GTDB-Tk removed Mash in v2.5.0
                              (now uses skani)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR           Extra arguments passed to gtdbtk
```

> **GTDB-Tk version note (latest 2.7.2):** since v2.5.0 GTDB-Tk screens with
> **skani** and no longer uses Mash, so `--mash_db` is gone. The old
> `--skip_ani_screen` was replaced by `--place_species` in v2.7.0. MaMISA parses
> both the old (`fastani_*`) and new (`closest_genome_*`) summary column names,
> so reports from any recent GTDB-Tk version work.

### run-mimag-rna

Checks the MIMAG rRNA (5S/16S/23S) + tRNA (≥18 amino acids) high-quality
criterion that CheckM2 cannot evaluate, via barrnap + tRNAscan-SE. Writes
`mimag_rna_summary.tsv`; feed it to `organize-mags --mimag-rna-dir` to gate HQ.

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --kingdom {bac,arc,euk}     Kingdom for barrnap/tRNAscan-SE (default: bac)
  --threads INT               Threads for barrnap (default: 1)
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --keep-intermediate         Keep per-genome barrnap/tRNAscan output files
```

### run-checkm2

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --extension STR             Genome file extension (default: fa)
  --threads INT               Threads for CheckM2 (default: 1)
  --database PATH             CheckM2 diamond database
  --force                     Overwrite existing CheckM2 output
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --checkm2-args STR          Extra arguments passed to checkm2
```

### run-gunc

Runs GUNC gene-level chimerism/contamination detection, complementing
`check-chimeras`. Feed the output back with `check-chimeras --gunc-dir`.
The command runs an env-health check and refuses to launch a GUNC broken by an
incompatible pandas/numpy/python (GUNC is old — see Installation).

```
Required (one of):
  --selected-dir PATH         Directory with HQ/, MQ/, LQ/ subdirectories
  --genome-dir PATH           Single directory with genome files

Required:
  -o, --output PATH           Output directory

Optional:
  --file-suffix STR           Genome file suffix (default: .fa)
  --threads INT               Threads (default: 1)
  --db-file PATH              GUNC diamond database (.dmnd); else uses $GUNC_DB
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gunc-args STR             Extra arguments passed to gunc run
```

### setup-workflow

Configures the Snakemake workflow: points `workflow/envs/*.yaml` at this checkout,
writes paths into `workflow/config.yaml`, validates database kinds, and (when env
names are given) sanity-checks them. `--install` builds every per-rule conda env.

```
Optional:
  --repo PATH                 MaMISA checkout (default: auto-detected)
  --workflow-dir PATH         Workflow dir (default: <repo>/workflow)
  --genomes-dir PATH          Input genomes directory (config)
  --genome-ext STR            Genome extension (config)
  --outdir PATH / --threads INT
  --gtdbtk-data DIR           GTDB-Tk data DIRECTORY
  --gunc-db FILE              GUNC .dmnd FILE
  --checkm2-db FILE           CheckM2 .dmnd FILE
  --sample-regex RE / --tax-level LEVEL
  --mamisa-env / --checkm2-env / --gtdbtk-env / --gunc-env NAME
                              Use these existing envs (sets use_named_envs)
  --install                   Build per-rule conda envs (snakemake)
  --no-check                  Skip the env sanity check
  --dry-run
```

### fetch-databases

Registers databases you already have (validated, written to config) or downloads
missing ones with each tool's own downloader. Path KIND is enforced: GTDB-Tk =
directory, GUNC and CheckM2 = `.dmnd` file.

```
Optional:
  --config PATH               config.yaml to update (default: workflow/config.yaml)
  --db-dir PATH               Where to download (default: ./databases)
  --tools LIST                Which DBs (default: gtdbtk,gunc,checkm2)
  --download                  Download any DB not supplied as an existing path
  --gtdbtk-data DIR           Existing GTDB-Tk data DIRECTORY
  --gunc-db FILE              Existing GUNC .dmnd FILE
  --checkm2-db FILE           Existing CheckM2 .dmnd FILE
  --gtdbtk-env / --gunc-env / --checkm2-env NAME   Env to run each downloader in
  --dry-run
```

### check-envs

Sanity-checks existing conda envs: reports each tool's version, whether MaMISA is
importable there, and (for GUNC) pandas/numpy/python compatibility. Prints the
exact install/fix command for anything that fails.

```
Optional (at least one):
  --mamisa-env NAME           Env with the light CLI tools + MaMISA
  --checkm2-env NAME          Env with CheckM2
  --gtdbtk-env NAME           Env with GTDB-Tk
  --gunc-env NAME             Env with GUNC
  --repo PATH                 MaMISA checkout used in fix hints
```

---

## Understanding Misassembly Detection

MaMISA uses soft-clipping information from read mapping to detect misassemblies. A
**clipping position** is a genomic location where reads are predominantly soft-clipped —
strong evidence that two unrelated sequences were joined during assembly.

The `classify-clipping` command adds mechanistic insight to each clipping position using
BAM-derived signals:

| Signal | What it detects |
|---|---|
| **depth_ratio** | Coverage spike (repeat collapse) or drop (deletion) |
| **discordant_fraction** | Reads whose mates map to a different contig or in wrong orientation |
| **large_insert_fraction** | Pairs with abnormally large insert sizes (> mean + 3σ) |
| **clipped_base_entropy** | Repetitive vs. diverse sequence at the break point |
| **sa_partner_contig** | Clipped reads whose supplementary alignment maps to one other contig → contig join (`chimeric_join`) |

When `check-read-chimeras` output is provided, taxonomy shifts near a clipping position
provide an additional, strong signal for chimera classification.

BAM alignment categories used throughout:
- **Depth counting**: primary + secondary (supplementary excluded, `-F 2048`)
- **Pair statistics**: primary only (FLAG `0x100 == 0`)
- **Supplementary**: excluded entirely

---

## Citation

If you use MaMISA in your research, please cite:

```
[citation here]
```

## License

MIT License — see LICENSE for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

Issues and feature requests: https://github.com/l-gallucci/mamisa/issues

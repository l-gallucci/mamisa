# MaMISA — Manage Misassemblies

A comprehensive toolkit for metagenomic assembly quality control and filtering.

## Overview

MaMISA provides a set of commands that cover the full workflow from raw assembly to
taxonomically classified, quality-filtered genomes:

| Command | Purpose |
|---|---|
| `process-large-contigs` | Extract, QC, and filter large contigs before binning |
| `check-read-chimeras` | Detect chimeric contigs via read-level taxonomy (Kraken2 + BAM) |
| `check-chimeras` | Detect chimeric MAGs via GC composition and GTDB-Tk signals |
| `classify-clipping` | Classify each clipping position with BAM evidence |
| `filter-misassemblies` | Split or remove contigs at misassembly positions |
| `remove-hq-contigs` | Remove HQ genome contigs from an assembly |
| `filter-checkm2` | Organise genomes into quality tiers from CheckM2 results |
| `run-gtdbtk` | Run GTDB-Tk taxonomy classification |

---

## Installation

```bash
git clone https://github.com/yourusername/mamisa.git
cd mamisa
pip install -e .

# Verify
mamisa --version
mamisa --help
```

### External tools (required per command)

| Command | Tool |
|---|---|
| `process-large-contigs` | CheckM2 |
| `check-read-chimeras` | samtools, Kraken2 |
| `check-chimeras` | (none beyond Python deps) |
| `classify-clipping` | samtools |
| `filter-misassemblies` | anvi'o (for upstream detection) |
| `filter-checkm2` | CheckM2 |
| `run-gtdbtk` | GTDB-Tk |

```bash
# anvi'o
conda create -n anvio-9 -c conda-forge -c bioconda anvio=9

# CheckM2
conda create -n checkm2 -c conda-forge -c bioconda checkm2

# GTDB-Tk
conda create -n gtdbtk-2 -c conda-forge -c bioconda gtdbtk
export GTDBTK_DATA_PATH=/path/to/gtdbtk_data

# samtools (for BAM-based commands)
conda install -c bioconda samtools
```

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
  --dry-run
```

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

### filter-checkm2

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

Other:
  --tiers LIST                Tiers to select (default: HQ,MQ,LQ)
  --extensions LIST           Genome extensions (default: fa,fasta,fna,...)
  --symlink | --copy          Link or copy files (default: symlink)
  --dry-run
```

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
  --mash-db PATH              GTDB-Tk mash database
  --tiers LIST                Tiers to process (default: HQ,MQ,LQ)
  --gtdbtk-args STR           Extra arguments passed to gtdbtk
```

---

## Understanding Misassembly Detection

MaMISA uses soft-clipping information from read mapping to detect misassemblies. A
**clipping position** is a genomic location where reads are predominantly soft-clipped —
strong evidence that two unrelated sequences were joined during assembly.

The `classify-clipping` command adds mechanistic insight to each clipping position using
four BAM-derived signals:

| Signal | What it detects |
|---|---|
| **depth_ratio** | Coverage spike (repeat collapse) or drop (deletion) |
| **discordant_fraction** | Reads whose mates map to a different contig or in wrong orientation |
| **large_insert_fraction** | Pairs with abnormally large insert sizes (> mean + 3σ) |
| **clipped_base_entropy** | Repetitive vs. diverse sequence at the break point |

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

Issues and feature requests: https://github.com/yourusername/mamisa/issues

# Pipeline

Binning and anvi'o misassembly detection are assumed done upstream for the
MAG-quality tail; the misassembly branch starts from an assembly plus a BAM.

## Run options

Both ways use the same conda envs (see [installation.md](installation.md)).

### Option A - run commands yourself (simplest)

Call each command in its env:

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

MIMAG HQ also requires 5S/16S/23S rRNA and tRNAs for at least 18 amino acids,
which CheckM2 does not check. `run-mimag-rna` adds that evidence and
`organize-mags --mimag-rna-dir` demotes HQ genomes that fail it to MQ. Skip it and
HQ is completeness/contamination only.

### Option B - Snakemake orchestrates it

One command runs the whole DAG. Configure once, then run:

```bash
conda install -c conda-forge -c bioconda snakemake-minimal   # in base or its own env

# tell the workflow which envs + databases to use (sanity-checks them)
mamisa setup-workflow --genomes-dir bins/ --genome-ext fa --threads 40 \
    --mamisa-env mamisa --checkm2-env checkm2 --gtdbtk-env gtdbtk --gunc-env gunc

# then run (named-env mode: Snakemake wraps each rule in `conda run -n`)
snakemake --cores 40 -s workflow/Snakefile --configfile workflow/config.yaml
```

Don't have the envs and want Snakemake to build them for you instead? Use
`mamisa setup-workflow --genomes-dir bins/ --install` and run with
`snakemake --use-conda ...`. Details in
[`../workflow/README.md`](../workflow/README.md). For clusters see
[hpc.md](hpc.md).

## Complete workflow (misassembly branch to taxonomy)

```
assembly.fa
    |
    v
[1] process-large-contigs       Separate very large contigs, run CheckM2 on each,
    |                           extract clean HQ genomes, return updated assembly.
    |
    v
[2] anvi-script-find-misassemblies   (external - anvi'o)
    |                           Detect soft-clipping positions in the BAM.
    |                           Produces  *-clipping.txt  files.
    |
    v
[3] check-read-chimeras         (optional, recommended)
    |                           BAM + Kraken2 per-read taxonomy.
    |                           Flags contigs whose reads come from >=2 organisms.
    |
    +-------------------------->[3b] check-chimeras   (optional)
    |                                GC + contamination + GUNC check on bins.
    |
    v
[4] classify-clipping           (optional, recommended)
    |                           Label each clipping position from BAM evidence:
    |                           end_artefact / repeat_collapse / deletion_artefact /
    |                           chimera_candidate / chimeric_join / sv_candidate.
    |
    v
[5] filter-misassemblies        Split or remove misassembled contigs.
    |                           Preserves HQ genomes. Integrates chimera report.
    |
    v
[6] binning                     (external - MetaBAT2, MaxBin2, etc.)
    |
    v
[7] run-checkm2                 (CheckM2) completeness/contamination
    |
    v
[8] organize-mags               Split bins into HQ / MQ / LQ quality tiers.
    |
    v
[9] run-gtdbtk                  Taxonomic classification of selected genomes.
```

## Step-by-step guide

### Step 1 - Process large contigs

Very long contigs (> 300 kbp by default) often represent single complete genomes
and should be handled separately before binning. This command separates large and
regular contigs, runs CheckM2 on each large contig, classifies each as extract HQ,
keep for splitting, or low quality, and returns an updated assembly.

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

Key output: `01_large_contigs/assembly_for_filtering.fa` (use in Step 5) and
`HQ_extracted/` (clean HQ genomes).

### Step 2 - Detect misassemblies with anvi'o (external)

```bash
anvi-script-find-misassemblies \
    -b mapping.bam \
    -o misassemblies/MisAsm \
    -T 40
```

Produces `MisAsm-clipping.txt` in `misassemblies/`. See
[references.md](references.md) (Trigodet et al. 2025) for the method.

### Step 3 - Detect chimeric contigs (read-level taxonomy)

Streams the BAM once, cross-references every read against Kraken2 per-read
taxonomy, builds a per-contig taxonomic profile, and flags contigs whose reads
originate from more than one organism.

```bash
mamisa check-read-chimeras \
    --bam mapping.bam \
    --kraken2-output kraken2_reads.txt \
    --kraken2-report kraken2_report.txt \
    -o 02_chimera/ \
    --min-mapq 20 --window 10000 --window-step 5000
```

Outputs `chimera_read_report.tsv` (per-contig) and `chimera_read_windows.tsv`
(sliding-window detail for classify-clipping).

Classifier choice: `--kraken2-output` (k-mer, fastest) or `--kaiju-output`
(protein-level, more sensitive on divergent/novel organisms). Give exactly one.
For names use `--kraken2-report` with Kraken2, or `--kaiju-names names.dmp` with
Kaiju. BLAST is not used here on purpose: per-read BLAST against nt is orders of
magnitude slower and short reads give low-specificity hits.

### Step 3b - Detect chimeric MAGs (GC + GUNC, optional)

Run after binning for a composition + gene-consistency check on complete bins.

```bash
mamisa check-chimeras \
    --bins-dir 06_bins/ \
    --output 02_chimera/gc_chimera_report.tsv \
    --gtdbtk-dir 09_taxonomy/ \
    --checkm2-report 07_checkm2/quality_report.tsv \
    --gunc-dir 07_gunc/ \
    --gc-window 5000 --gc-step 2500
```

### Step 4 - Classify clipping positions

Characterise each anvi'o-reported clipping position using BAM evidence so you can
prioritise which splits are biologically meaningful.

```bash
mamisa classify-clipping \
    --bam mapping.bam \
    --misassemblies misassemblies/ \
    --taxonomy-windows 02_chimera/chimera_read_windows.tsv \
    --output 03_classified/clipping_classified.tsv
```

Classification labels:

| Label | Biological meaning |
|---|---|
| `end_artefact` | Near contig terminus; assembly edge noise |
| `repeat_collapse` | Coverage spike + low-entropy clipped bases; collapsed repeat |
| `deletion_artefact` | Local depth drop; internal deletion or coverage collapse |
| `chimera_candidate` | Discordant pairs + large inserts, possibly a taxonomy shift |
| `chimeric_join` | Clipped reads' supplementary alignment (`SA:Z:`) maps to one other contig: the two contigs are joined here (Trigodet et al. 2025 read-level signal) |
| `sv_candidate` | SV signal (discordant) without depth anomaly or taxon shift |
| `low_confidence` | Ambiguous or insufficient evidence |

### Step 5 - Filter misassemblies

```bash
# Conservative: use all anvi'o-reported positions (--safe, the default)
mamisa filter-misassemblies \
    --assembly 01_large_contigs/assembly_for_filtering.fa \
    --misassemblies misassemblies/ \
    --hq-genomes 01_large_contigs/HQ_extracted/ \
    --output 04_clean/assembly_clean.fa \
    --mode split --preserve-hq-with-issues \
    --chimera-report 02_chimera/chimera_read_report.tsv \
    --split-hq-circular --safe
```

`--safe` uses every anvi'o-reported position; `--min-clip-coverage N` is the
selective alternative (only split where at least N reads are clipped).

### Step 6 - Binning (external)

```bash
jgi_summarize_bam_contig_depths --outputDepth depth.txt mapping.bam
metabat2 -i 04_clean/assembly_clean.fa -a depth.txt -o 05_bins/bin -t 40
```

### Step 7 - Quality assessment with CheckM2

```bash
conda run -n checkm2 mamisa run-checkm2 --genome-dir 05_bins/ -o 06_checkm2/ --threads 40
```

### Step 8 - Organise genomes by quality

```bash
mamisa organize-mags \
    --checkm2-root 06_checkm2/ --genomes-dir 05_bins/ \
    --output 07_filtered/ --tiers HQ,MQ --symlink
```

MiMAG quality thresholds (defaults):

| Tier | Completeness | Contamination |
|---|---|---|
| HQ | >= 90% | <= 5% |
| MQ | >= 70% | <= 10% |
| LQ | >= 50% | <= 10% |

The defaults split the MIMAG medium band into MQ/LQ. For strict MIMAG tiers run
with `--mq-comp-min 50 --mq-cont-max 10 --lq-comp-min 0 --lq-cont-max 10`. Full
MIMAG HQ also needs rRNA/tRNA (`--mimag-rna-dir`, see Option A).

### Step 9 - Taxonomic classification with GTDB-Tk

```bash
conda run -n gtdbtk mamisa run-gtdbtk \
    --selected-dir 07_filtered/Selected/ \
    --output 08_taxonomy/ --cpus 40 --tiers HQ,MQ
```

## Understanding misassembly detection

MaMISA uses soft-clipping information from read mapping to detect misassemblies. A
clipping position is a genomic location where reads are predominantly soft-clipped:
strong evidence that two unrelated sequences were joined during assembly.

`classify-clipping` adds mechanistic insight to each position using BAM signals:

| Signal | What it detects |
|---|---|
| depth_ratio | Coverage spike (repeat collapse) or drop (deletion) |
| discordant_fraction | Reads whose mates map to a different contig or wrong orientation |
| large_insert_fraction | Pairs with abnormally large insert sizes (> mean + 3 sd) |
| clipped_base_entropy | Repetitive vs diverse sequence at the break point |
| sa_partner_contig | Clipped reads whose supplementary alignment maps to one other contig: contig join (`chimeric_join`) |

When `check-read-chimeras` output is provided, taxonomy shifts near a clipping
position provide an additional, strong signal for chimera classification.

BAM alignment categories used throughout:
- Depth counting: primary + secondary (supplementary excluded, `-F 2048`)
- Pair statistics: primary only (FLAG `0x100 == 0`)
- Supplementary: excluded entirely

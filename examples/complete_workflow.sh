#!/usr/bin/env bash
#
# MaMISA Complete Workflow Example
# This script demonstrates a full pipeline from assembly QC to classified genomes
#

set -euo pipefail

# Configuration
ASSEMBLY="data/assembly.fa"
MAPPING_BAM="data/mapping.bam"
WORKDIR="mamisa_workflow"
THREADS=40

echo "======================================================================"
echo "  MaMISA Complete Workflow Example"
echo "======================================================================"

# Create working directory
mkdir -p "$WORKDIR"
cd "$WORKDIR"

#------------------------------------------------------------------------------
# STEP 1: Detect misassemblies with anvi'o
#------------------------------------------------------------------------------
echo ""
echo "STEP 1: Detecting misassemblies with anvi'o..."
echo "----------------------------------------------------------------------"

mkdir -p 01_misassemblies

# Run anvi-script-find-misassemblies (requires anvi'o)
if command -v anvi-script-find-misassemblies &> /dev/null; then
    anvi-script-find-misassemblies \
        -b "../$MAPPING_BAM" \
        -o 01_misassemblies/MisAsm \
        -T "$THREADS"
    
    echo "✓ Misassembly detection completed"
    echo "  Output: 01_misassemblies/"
    ls -lh 01_misassemblies/
else
    echo "⚠ anvi-script-find-misassemblies not found, skipping..."
    echo "  Please install anvi'o to run this step"
fi

#------------------------------------------------------------------------------
# STEP 2: Initial binning (placeholder - use your favorite binner)
#------------------------------------------------------------------------------
echo ""
echo "STEP 2: Initial binning..."
echo "----------------------------------------------------------------------"
echo "ℹ Run your favorite binner here (MetaBAT2, MaxBin2, CONCOCT, etc.)"
echo "  For example:"
echo "  metabat2 -i assembly.fa -a depth.txt -o bins/bin -t $THREADS"
echo ""
echo "Assuming bins are in: 02_initial_bins/"

mkdir -p 02_initial_bins
# ... your binning commands here ...

#------------------------------------------------------------------------------
# STEP 3: Quality assessment with CheckM2
#------------------------------------------------------------------------------
echo ""
echo "STEP 3: Assessing genome quality with CheckM2..."
echo "----------------------------------------------------------------------"

mkdir -p 03_checkm2

# Run CheckM2 (requires checkm2)
if command -v checkm2 &> /dev/null; then
    checkm2 predict \
        --threads "$THREADS" \
        --input 02_initial_bins/ \
        --output-directory 03_checkm2/ \
        -x fa
    
    echo "✓ CheckM2 completed"
    echo "  Output: 03_checkm2/quality_report.tsv"
else
    echo "⚠ CheckM2 not found, skipping..."
    echo "  Please install CheckM2 to run this step"
fi

#------------------------------------------------------------------------------
# STEP 4: Filter genomes by quality with MaMISA
#------------------------------------------------------------------------------
echo ""
echo "STEP 4: Filtering genomes by quality..."
echo "----------------------------------------------------------------------"

mkdir -p 04_filtered_genomes

mamisa filter-checkm2 \
    --checkm2-root 03_checkm2/ \
    --genomes-dir 02_initial_bins/ \
    --output 04_filtered_genomes/ \
    --tiers HQ,MQ \
    --hq-comp-min 90 \
    --hq-cont-max 5 \
    --mq-comp-min 70 \
    --mq-cont-max 10 \
    --symlink

echo "✓ Quality filtering completed"
echo "  HQ genomes: 04_filtered_genomes/Selected/HQ/"
echo "  MQ genomes: 04_filtered_genomes/Selected/MQ/"

#------------------------------------------------------------------------------
# STEP 5: Remove HQ contigs from assembly
#------------------------------------------------------------------------------
echo ""
echo "STEP 5: Removing HQ contigs from assembly..."
echo "----------------------------------------------------------------------"

mkdir -p 05_assembly_filtered

mamisa remove-hq-contigs \
    --assembly "../$ASSEMBLY" \
    --hq-dir 04_filtered_genomes/Selected/HQ/ \
    --output 05_assembly_filtered/assembly_no_hq.fa \
    --min-length 1000 \
    --stats 05_assembly_filtered/removal_stats.txt

echo "✓ HQ contigs removed"
echo "  Filtered assembly: 05_assembly_filtered/assembly_no_hq.fa"

#------------------------------------------------------------------------------
# STEP 6: Filter assembly for misassemblies
#------------------------------------------------------------------------------
echo ""
echo "STEP 6: Filtering assembly for misassemblies..."
echo "----------------------------------------------------------------------"

mkdir -p 06_assembly_clean

# First, do a dry-run to see statistics
echo "Running dry-run to preview changes..."
mamisa filter-misassemblies \
    --assembly 05_assembly_filtered/assembly_no_hq.fa \
    --misassemblies 01_misassemblies/ \
    --dry-run \
    --min-length 2500

echo ""
echo "Proceeding with actual filtering..."

mamisa filter-misassemblies \
    --assembly 05_assembly_filtered/assembly_no_hq.fa \
    --misassemblies 01_misassemblies/ \
    --output 06_assembly_clean/assembly_clean.fa \
    --mode split \
    --min-length 2500 \
    --stats 06_assembly_clean/filtering_stats.txt

echo "✓ Assembly filtering completed"
echo "  Clean assembly: 06_assembly_clean/assembly_clean.fa"

#------------------------------------------------------------------------------
# STEP 7: Re-bin with cleaned assembly (optional)
#------------------------------------------------------------------------------
echo ""
echo "STEP 7: Re-binning with cleaned assembly (optional)..."
echo "----------------------------------------------------------------------"
echo "ℹ You can now re-run binning on the cleaned assembly"
echo "  This may yield better quality bins"
echo ""
echo "Example:"
echo "  metabat2 -i 06_assembly_clean/assembly_clean.fa -a depth.txt -o bins_v2/bin"

#------------------------------------------------------------------------------
# STEP 8: Taxonomic classification with GTDB-Tk
#------------------------------------------------------------------------------
echo ""
echo "STEP 8: Taxonomic classification with GTDB-Tk..."
echo "----------------------------------------------------------------------"

mkdir -p 08_taxonomy

# Classify HQ and MQ genomes
if command -v gtdbtk &> /dev/null; then
    mamisa run-gtdbtk \
        --selected-dir 04_filtered_genomes/Selected/ \
        --output 08_taxonomy/ \
        --extension fa \
        --cpus "$THREADS" \
        --tiers HQ,MQ
    
    echo "✓ Taxonomic classification completed"
    echo "  Results: 08_taxonomy/HQ/ and 08_taxonomy/MQ/"
else
    echo "⚠ GTDB-Tk not found, skipping..."
    echo "  Please install GTDB-Tk to run this step"
fi

#------------------------------------------------------------------------------
# STEP 9: Generate summary report
#------------------------------------------------------------------------------
echo ""
echo "STEP 9: Generating summary report..."
echo "----------------------------------------------------------------------"

REPORT="workflow_summary.txt"

cat > "$REPORT" <<EOF
====================================================================
  MaMISA Workflow Summary
====================================================================

Workflow completed: $(date)
Working directory: $(pwd)

INPUT:
  Assembly: $ASSEMBLY
  Mapping BAM: $MAPPING_BAM

OUTPUTS:

1. Misassemblies Detection
   Directory: 01_misassemblies/
   Files: MisAsm-clipping.txt, MisAsm-zero_cov.txt

2. Initial Binning
   Directory: 02_initial_bins/
   Bins: $(find 02_initial_bins/ -type f -name "*.fa" 2>/dev/null | wc -l)

3. CheckM2 Quality Assessment
   Directory: 03_checkm2/
   Report: quality_report.tsv

4. Filtered Genomes by Quality
   Directory: 04_filtered_genomes/Selected/
   HQ genomes: $(find 04_filtered_genomes/Selected/HQ/ -type f 2>/dev/null | wc -l)
   MQ genomes: $(find 04_filtered_genomes/Selected/MQ/ -type f 2>/dev/null | wc -l)

5. Assembly without HQ Contigs
   File: 05_assembly_filtered/assembly_no_hq.fa
   Stats: 05_assembly_filtered/removal_stats.txt

6. Clean Assembly (misassemblies filtered)
   File: 06_assembly_clean/assembly_clean.fa
   Stats: 06_assembly_clean/filtering_stats.txt
   Contigs: $(grep -c "^>" 06_assembly_clean/assembly_clean.fa 2>/dev/null || echo "N/A")

8. Taxonomic Classification
   Directory: 08_taxonomy/
   HQ: 08_taxonomy/HQ/
   MQ: 08_taxonomy/MQ/

====================================================================
NEXT STEPS:
  1. Review quality_report.tsv for genome completeness/contamination
  2. Check filtering_stats.txt for assembly filtering details
  3. Examine GTDB-Tk taxonomy in 08_taxonomy/
  4. Consider re-binning with the cleaned assembly
  5. Perform downstream analyses on filtered genomes
====================================================================
EOF

cat "$REPORT"
echo ""
echo "✓ Workflow completed successfully!"
echo "  Summary report: $REPORT"

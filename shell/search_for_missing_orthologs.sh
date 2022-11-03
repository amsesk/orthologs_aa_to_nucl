#!/bin/bash

BASE="/Users/amsesk/DATA/bonitomes"
BINDIR="/Users/amsesk/dev/orthologs_aa_to_nucl"

GFFPATH="${BASE}/gff/"
ASMPATH="${BASE}/asm/"
ISOLATES_PATH="${BASE}/EHB_Isolates.xlsx"
ORTHOTABLE_PATH="${BASE}/all_marker_orthologs_by_ltp.tsv"
OUTPUT_MISSING_PATH="ltp_missing_orthologs_by_orthogroup.tsv"
MULTIMATCH_PROTS_PATH="multimatch_prots.txt"

TARGET=release
#TARGET=debug

currdir=$(pwd)
cd $BINDIR
cargo build --release
cd $currdir

#python ~/dev/scripts/xlsx_to_tsv.py $ISOLATES_PATH | cut -f4,13 | sed 's/\t/.gff\t/' | sed "s#^#$GFFPATH#" | sed "s#\t#\t${ASMPATH}#" | \
#	~/dev/orthologs_aa_to_nucl/target/$TARGET/orthologs_aa_to_nucl -o $ORTHOTABLE_PATH

Rscript --vanilla ${BINDIR}/R/orthologs_missing_from.R -o ${ORTHOTABLE_PATH} -r ${MULTIMATCH_PROTS_PATH} -p ${OUTPUT_MISSING_PATH}

python ${BINDIR}/python/ltp_occupancy_into_blast_commands.py -o ${OUTPUT_MISSING_PATH} -i ${ISOLATES_PATH} -g /nfs4/BPP/Uehling_Lab/data/EHB/genomes -f /nfs4/BPP/Uehling_Lab/amsesk/Bonitomes/phylogenomics_all/round3/phylo_marker_loss/marker_ortholog_nucl_fasta

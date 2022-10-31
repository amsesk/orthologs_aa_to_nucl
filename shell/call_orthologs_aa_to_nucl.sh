#!/bin/bash

BASE="/Users/amsesk/DATA/bonitomes"
BINDIR="/Users/amsesk/dev/orthologs_aa_to_nucl"

GFFPATH="${BASE}/gff/"
ASMPATH="${BASE}/asm/"
ISOLATES_PATH="${BASE}/EHB_Isolates.xlsx"
ORTHOTABLE_PATH="${BASE}/all_marker_orthologs_by_ltp.tsv"

TARGET=release
#TARGET=debug

currdir=$(pwd)
cd $BINDIR
cargo build --release
cd $currdir

python ~/dev/scripts/xlsx_to_tsv.py $ISOLATES_PATH | cut -f4,13 | sed 's/\t/.gff\t/' | sed "s#^#$GFFPATH#" | sed "s#\t#\t${ASMPATH}#" | \
	~/dev/orthologs_aa_to_nucl/target/$TARGET/orthologs_aa_to_nucl -o $ORTHOTABLE_PATH 

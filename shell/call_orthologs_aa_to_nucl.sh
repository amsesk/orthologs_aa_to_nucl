GFFPATH="/home/aimzez/DATA/bonitomes/gff/"
ASMPATH="/home/aimzez/DATA/bonitomes/asm/"
ORTHOTABLE_PATH="/home/aimzez/DATA/bonitomes/all_marker_orthologs_by_ltp.tsv"
currdir=$(pwd)
cd /home/aimzez/dev/orthologs_aa_to_nucl
cargo build 
cd $currdir

python ~/dev/scripts/xlsx_to_tsv.py /home/aimzez/DATA/bonitomes/EHB_Isolates.xlsx | cut -f4,13 | sed 's/\t/.gff\t/' | sed "s#^#$GFFPATH#" | sed "s#\t#\t${ASMPATH}#" | \
	~/dev/orthologs_aa_to_nucl/target/debug/orthologs_aa_to_nucl -o $ORTHOTABLE_PATH 

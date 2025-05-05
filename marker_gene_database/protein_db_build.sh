#!/bin/bash

# set the sensitivity of DIAMOND
while getopts 's:' opt ; do
        case $opt in
                s) SENS=$OPTARG ;;
        esac
done

# skip over the processed options
shift $((OPTIND-1))


#### parameter check
if [ $# -lt 4 ]
then
echo "Usage: $0 protein_name seed_path dataset_path num_threads"
echo ""
echo "run this in the protein_tools conda env"
echo ""
echo "requires 4 arguments:"
echo "protein_name: 	name of the protein you're building the database for"
echo "seed_path: 	path to the fasta file containing 1 or multiple sequences to use as seed for database construction"
echo "dataset_path: 	path to the dataset you are trawling for homologs, in fasta format"
echo "num_threads: 	number of theads to use for DIAMOND"
echo ""
echo "optional argument -s:"
echo "sensitivity:	DIAMOND sensitivity, choose: --sensitive, --very-sensitive, --ultra-sensitive"
echo "if not specified, DIAMOND will run in fast mode"
exit 1
fi

PROT="$1"
SEED="$2"
DB="$3"
THREADS="$4"

# do initial DIAMOND search
mkdir "$PROT"_protein_db
diamond makedb --in $SEED -d "$PROT"_protein_db/"$PROT"_seed_db -p $THREADS
diamond blastp -d "$PROT"_protein_db/"$PROT"_seed_db -q $DB -p $THREADS -o "$PROT"_protein_db/"$PROT"_hits -k 1 --matrix blosum45 --masking 0 --outfmt 6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore score -b 6 -c 2 --min-score 50 --comp-based-stats 0 $SENS

# read lookup
cd "$PROT"_protein_db
/lisc/project/dome/speth/resources/scripts/other/blast_based_read_lookup_new.pl "$PROT"_hits $DB "$PROT"_hits.faa

# 2nd diamond search
/lisc/project/dome/speth/resources/scripts/other/calc_max_score.pl "$PROT"_hits.faa BLOSUM45 "$PROT"_self

# create BSR table for inspection
/lisc/project/dome/speth/resources/scripts/other/bsr_calc.pl "$PROT"_hits "$PROT"_self "$PROT"_bsr

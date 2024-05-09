#! /bin/bash

eval "$(conda shell.bash hook)"
conda deactivate

# set defaults for the options
ASSEMBLY="assembly"
MAP="read_mapping"
ID=98
LEN=80

while getopts 'a:m:i:l:' opt ; do
    case $opt in
        a) ASSEMBLY=$OPTARG ;;
        m) MAP=$OPTARG ;;
        i) ID=$OPTARG ;;
        l) LEN=$OPTARG ;;
    esac
done

# skip over the processed options
shift $((OPTIND-1))

#### parameter check
if [ $# -lt 6 ]
then
echo "Usage: $0 seed_contigs_path fw_reads_path rv_reads_path project_name no_of_iterations threads"
echo ""
echo "requires 6 arguments:"
echo "contigs:		absolute path to fasta file with seed contigs"
echo "fw reads:		absolute path to the forward reads"
echo "rv reads:		absolute path to the reverse reads"
echo "project name:	prefix for files"
echo "num iterations: 	number of iterations to run the loop"
echo "num threads: 	number of threads to use where appropriate"
echo ""
echo "optional arguments:"
echo "-a name of conda env with read mappers, default: read_mapping"
echo "-m name of conda env with spades, default: assembly"
echo "-i value for alignment identity of reads to keep, default 98"
echo "-l value for alignment length fraction of reads to keep, default 80"
exit 1
fi


SEED=$1
FW=$2
RV=$3
NAME=$4
END=$5
THREADS=$6

# create dir for files
mkdir it_1
cd it_1

# map reads to the seed contigs
conda activate ${MAP}
minimap2 -t ${THREADS} -a -x sr --sam-hit-only --no-pair ${SEED} ${FW} ${RV} > ${NAME}_mapped.sam 2> minimap_stderr

# create sorted bam files for processing
samtools view -bShu ${NAME}_mapped.sam | samtools sort -m 97G -@ 3 - -o ${NAME}_mapped_sorted.bam
samtools index ${NAME}_mapped_sorted.bam

# filter mapping to only contain matches above stringent threshold
coverm filter --bam-files ${NAME}_mapped_sorted.bam -o ${NAME}_mapped_sorted_filtered.bam --min-read-percent-identity $ID --min-read-aligned-percent $LEN
samtools sort -o ${NAME}_mapped_sorted_filtered_sort.bam --threads ${THREADS} ${NAME}_mapped_sorted_filtered.bam
samtools fastq -1 ${NAME}_map_98_80_R1.fastq -2 ${NAME}_map_98_80_R2.fastq ${NAME}_mapped_sorted_filtered_sort.bam

# retrieve read pairs with at least one of the reads matching the seed contigs
seqkit seq -n ${NAME}_map_98_80_R1.fastq ${NAME}_map_98_80_R2.fastq | sort | uniq > ${NAME}_map_98_80_ids
seqkit grep -f ${NAME}_map_98_80_ids $FW > ${NAME}_map_98_80_pair1.fastq
seqkit grep -f ${NAME}_map_98_80_ids $RV > ${NAME}_map_98_80_pair2.fastq

conda deactivate

# assemble matching reads
conda activate ${ASSEMBLY}
spades.py --isolate -t ${THREADS} -1 ${NAME}_map_98_80_pair1.fastq -2 ${NAME}_map_98_80_pair2.fastq -k 21,33,55,77,99,111,127 -o spades_${NAME} > spades_stdout 2>&1
conda deactivate

cd spades_${NAME}
fasta_to_gc_cov_length_tab.pl scaffolds.fasta scaf.tab
awk '$3 >= 2000 {print $0}' scaf.tab > scaf_filtered.tab
tab_to_fasta.pl scaf_filtered.tab scaf_filtered.fa 3
cd ../../

START=2
for ITERATION in $(eval echo "{${START}..${END}") ; do
	mkdir it_${ITERATION}
	cp it_$((${ITERATION}-1))/spades_${NAME}/scaf_filtered.fa it_${ITERATION}
	cd it_${ITERATION}
 
	conda activate ${MAP}
  minimap2 -t $THREADS -a -x sr --score-N 2 --sam-hit-only --no-pair scaf_filtered.fa ${FW} ${RV} > ${NAME}_mapped.sam 2> minimap_stderr
	samtools view -bShu ${NAME}_mapped.sam | samtools sort -m 97G -@ 3 - -o ${NAME}_mapped_sorted.bam
	samtools index ${NAME}_mapped_sorted.bam

	coverm filter --bam-files ${NAME}_mapped_sorted.bam -o ${NAME}_mapped_sorted_filtered.bam --min-read-percent-identity 98 --min-read-aligned-percent 80
	samtools sort -o ${NAME}_mapped_sorted_filtered_sort.bam --threads ${THREADS} ${NAME}_mapped_sorted_filtered.bam
	samtools fastq -1 ${NAME}_map_98_80_R1.fastq -2 ${NAME}_map_98_80_R2.fastq ${NAME}_mapped_sorted_filtered_sort.bam

  seqkit seq -n ${NAME}_map_98_80_R1.fastq ${NAME}_map_98_80_R2.fastq | sort | uniq > ${NAME}_map_98_80_ids
  seqkit grep -f ${NAME}_map_98_80_ids ${FW} > ${NAME}_map_98_80_pair1.fastq
  seqkit grep -f ${NAME}_map_98_80_ids ${RV} > ${NAME}_map_98_80_pair2.fastq
	conda deactivate

	conda activate ${ASSEMBLY}
  spades.py -t ${THREADS} --isolate -1 ${NAME}_map_98_80_pair1.fastq -2 ${NAME}_map_98_80_pair2.fastq -k 21,33,55,77,99,111,127 -o spades_${NAME} > spades_stdout 2>&1
  conda deactivate

	cd spades_${NAME}
	fasta_to_gc_cov_length_tab.pl scaffolds.fasta scaf.tab
	awk '$3 >= 2000 {print $0}' scaf.tab > scaf_filtered.tab
	tab_to_fasta.pl scaf_filtered.tab scaf_filtered.fa 3
	cd ../../
done

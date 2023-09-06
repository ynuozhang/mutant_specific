
HOME_LOC=/data/rozen/home/e0833634
PROJECT_LOC=$HOME_LOC/mutant_specific/cd-hit_homology_reduction
FASTA_LOC=$HOME_LOC/mutant_specific/cd-hit_homology_reduction
CD_HIT_2D=$HOME_LOC/packages/cdhit/cd-hit-2d
HOMOLOGY_THRESHOLD=0.9

OUTPUT_LOC=$PROJECT_LOC/cdhit_2d_output/cdhit_2d_output.txt
#if [ -d "$OUTPUT_LOC" ]; then
 # echo "Folder $OUTPUT_LOC exists. Using it as the output location."
#else
#  mkdir "$OUTPUT_LOC"
#  echo "Folder $OUTPUT_LOC created for output."
#fi

FASTA_DB1_WT1=$FASTA_LOC/SKEMPI_wild_seq_1.fasta
FASTA_DB1_WT2=$FASTA_LOC/SKEMPI_wild_seq_2.fasta
FASTA_DB2_WT1=$FASTA_LOC/PPIs_partial_21747_sequences_1.fasta
FASTA_DB2_WT2=$FASTA_LOC/PPIs_partial_21747_sequences_2.fasta

$CD_HIT_2D -i $FASTA_DB1_WT1 -i2 $FASTA_DB2_WT1 -o $OUTPUT_LOC -c $HOMOLOGY_THRESHOLD -n 5 -d 0 -M 0 -T 0 -p 1 -g 1 -s2 0.0 -S2 1500
# -n 5: word_length
# -M: memory limit in MB, 0 for unlimited
# -T: number of threads, 0 for usage of all CPUs
# b: band_width of alignment, left+right how many positions to check, default is 20
# -p: print alignment overlap in .clstr file if set to 1
# -g:  by cd-hit's default algorithm, a sequence is clustered to the first
         # 	cluster that meet the threshold (fast cluster). If set to 1, the program
         # 	will cluster it into the most similar cluster that meet the threshold
         # 	(accurate but slow mode)
         # 	but either 1 or 0 won't change the representatives of final clusters
# by default, only matches where sequences in db2 are not longer than sequences in db1
# -s2: length difference cutoff for db1, default 1.0
    #    by default, seqs in db1 >= seqs in db2 in a same cluster
    #    if set to 0.9, seqs in db1 may just >= 90% seqs in db2
# -S2: length difference cutoff, default 0
    #    by default, seqs in db1 >= seqs in db2 in a same cluster
    #    if set to 60, seqs in db2 may 60aa longer than seqs in db1


SCRIPTS_DIR=$1

for file in `cd fastq/; ls *R1.fq.gz`; 
do 
PREFIX=${file%%_R1.fq.gz}
echo -e "nice $SCRIPTS_DIR\c"
echo $file | awk -v prefix="$PREFIX" -F "_" '{print("/scripts/exb_decoder_parallel_qscore.py -a fastq/"$0" -b fastq/"prefix"_R2.fq.gz -o analysis/exb_stats/"prefix".dec.out -n 40")}'; 

done

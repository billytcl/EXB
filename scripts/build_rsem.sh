SCRIPTS_DIR=$1

#INSTRUCTIONS:
#-change location of RSEM
#-change location of RSEM reference

for file in `cd analysis/dedup_star; ls *toTranscriptome.out.bam`;
do
	prefix=`basename $file .Aligned.toTranscriptome.out.bam`
	echo nice ~/RSEM-1.2.22/rsem-calculate-expression --bam --estimate-rspd --no-bam-output --seed 12345 -p 40 --paired-end analysis/dedup_star/$file /mnt/ix2/research_copy/sep_people/20131229_billy_projects/genomes/rsem/Homo_sapiens.GRCh37.75.RSEM analysis/dedup_rsem/$prefix.rsem --no-qualities
done



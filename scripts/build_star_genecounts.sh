SCRIPTS_DIR=$1

#INSTRUCTIONS:
#
#-change location of STAR reference (--genomeDir)
#-change location of STAR

for file in `cd analysis/dedup/; ls *.fa.gz`; do 
	prefix=`echo $file | awk -F "." '{print $1}'`
	echo nice ~/star --genomeDir /mnt/ix2/research_copy/sep_people/20131229_billy_projects/genomes/star_ensembl_grch37.75_ercc/ --outFileNamePrefix analysis/dedup_star/$prefix. --runThreadN 40 --readFilesIn "<(zcat analysis/dedup/$file | paste - - | awk '{if (NR % 2) print(\$1\"\\n\"\$2);}')" "<(zcat analysis/dedup/$file | paste - - | awk '{if (NR % 2 == 0) print(\$1\"\\n\"\$2);}')" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0
done;

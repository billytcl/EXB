SCRIPTS_DIR=$1

for file in `cd analysis/exb_stats/; ls *.cd.out`; do 
	prefix=`echo $file | awk -F "." '{print $1}'`
	echo -e "nice $SCRIPTS_DIR\c"
	echo "/scripts/exb_extract_decoder.py -i analysis/exb_stats/$file | cutadapt -b CTGTCTCTTATACACATCT - | paste - - - - | awk '{if ((\$2 != \"\") && (\$4 != \"\")) print(\$1\"\\n\"\$2\"\\n\"\$3\"\\n\"\$4)}' | gzip > analysis/dedup/$prefix.fa.gz"
done

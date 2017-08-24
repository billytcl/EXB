SCRIPTS_DIR=$1

for file in `cd analysis/exb_stats/; ls *.dec.out`; 
do 

echo -e "nice $SCRIPTS_DIR\c"

echo $file | awk -F "." '{print("/scripts/exb_consensus_parallel.py -i analysis/exb_stats/"$1".dec.out -o analysis/exb_stats/"$1".cd.out -n 40")}'; 

done;

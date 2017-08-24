#!/bin/bash

#rename this to analysis folder
SCRIPTS_DIR='/mnt/ix2/20170824_BL_EXB_scripts/'

#make analysis directories
mkdir analysis
cd analysis
mkdir exb_stats
mkdir dedup
mkdir dedup_star
mkdir dedup_rsem
cd ..


#build module scripts and parse infinity codes
$SCRIPTS_DIR/scripts/build_decoder_qual.sh $SCRIPTS_DIR > dec.sh
chmod a+x dec.sh
./dec.sh

$SCRIPTS_DIR/scripts/build_consensus_decoder.sh $SCRIPTS_DIR > consensus.sh
chmod a+x consensus.sh
./consensus.sh

$SCRIPTS_DIR/scripts/build_extract_decoder.sh $SCRIPTS_DIR > dedup_decoder.sh
chmod a+x dedup_decoder.sh
./dedup_decoder.sh

$SCRIPTS_DIR/scripts/build_star_genecounts.sh $SCRIPTS_DIR > dedup_star.sh
chmod a+x dedup_star.sh
./dedup_star.sh

$SCRIPTS_DIR/scripts/build_rsem.sh $SCRIPTS_DIR > dedup_rsem.sh
chmod a+x dedup_rsem.sh
./dedup_rsem.sh



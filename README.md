# EXB
## processing scripts for EXB workflow

## Instructions:

- Clone the repository (eg. git clone https://github.com/billytcl/EXB.git).
- Copy `scripts/execute_decoder_qscore.sh` to an analysis folder of your choice.
- Create a folder inside the analysis directory named `fastq`. Copy the fastq files inside, ensuring that they are in the name format `*.fq.gz`. See the `test` folder for the recommended folder structure.
- Edit execute_decoder_qscore.sh and replace `SCRIPTS_DIR` variable to the cloned github folder. Don't include "scripts" in the folder name.
- Edit `build_star_genecounts.sh` so that it points to where STAR is installed. Also change the location to where your STAR reference is located.
- Edit `build_rsem.sh` so that it points to where RSEM is installed. Also change the location to where the RSEM reference is located.


#! /bin/bash
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH --job-name="cleanup"
#SBATCH -t 00:20:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=500Mb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu


cd /projects/b1052/mckenna/obrien_meta/results

rm -r ./megahit_coassembly/intermediate_contigs
rm -r ./metabat_checkm_coassembly/*/bins
rm -r ./metabat_checkm_coassembly/*/storage
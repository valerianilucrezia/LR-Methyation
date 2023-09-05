#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=R
#SBATCH --nodes=1
#SBATCH --mem 100g
#SBATCH --time=2:00:00

module load R/4.2.3

np=/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/input/all.cpg.bed
epic=/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/input/beta_minfi_Noob.csv
out=/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/LR-Methyation/res
name=prova_minfi_Noob


cd /orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/LR-Methyation

Rscript ./scripts/main.R --nanopore $np --epic $epic --output $out --sample $name

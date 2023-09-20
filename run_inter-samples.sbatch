#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=R
#SBATCH --nodes=1
#SBATCH --mem 100g
#SBATCH --time=2:00:00

module load R/4.2.3

samplelist=/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/LR-Methyation/data/sample_list
out=/orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/LR-Methyation/res

cd /orfeo/LTS/LADE/LT_storage/lvaleriani/methylation/prova_GEL/LR-Methyation

Rscript ./scripts/inter_sample.R --samplelist $samplelist --out $out

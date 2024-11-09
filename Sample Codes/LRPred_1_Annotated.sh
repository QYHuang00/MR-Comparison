#! /bin/bash

# This script uses the LDPred method to generate Polygenic Risk Scores (PRS) from PLINK-format genotype data.

# Reference:
# Vilhjálmsson, Bjarni J., Stephan Ripke, Peter Kraft, et al. Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. The American Journal of Human Genetics, vol. 97, no. 4, 2015, pp. 576-592.

# Related GitHub repository:
# https://github.com/bvilhjal/ldpred

# Set the working directory where input and output files are located.
dir="/gpfs/gibbs/project/zhao/qh63/R/PLink_data/output_1"
cd $dir

# Loop through 100 iterations (from 1 to 100), processing each dataset separately.
for i in {1..100}
do

  ldpred coord --gf ./test_1_${i} --ssf ./ssf_1_${i}.txt --N 5000 --out ./output/output_${i} --vbim ./test_1_${i}.bim --vgf ./test_1_${i} --ssf-format LDPRED
  ldpred gibbs --cf ./output/output_${i} --ldr 3 --ldf ./output/ld_1_${i} --out ./output/output_weight_1_${i} --N 5000
  ldpred score --gf ./test_1_${i} --rf ./output/output_weight_1_${i} --out ./output/output_score_1_${i} --only-score
  
  ldpred score --gf ./test_2_${i} --rf ./output/output_weight_1_${i} --out ./output/output_score_2_${i} --only-score

done


# This folder contains results of performing a GWAS using Regenie (both steps 1 and 2).
  - Step 1: Training the model.
      - It consists of building a ridge regression model using genotyped variants to predict phenotypic values.
      - This accounts for polygenic effects, relatedness, and population structure while avoiding issues with large-scale mixed-model approaches. (corrects confoundness)
  - Step 2: Performing the GWAS itself.
      - SNP association with the trait based on the model from step 1.

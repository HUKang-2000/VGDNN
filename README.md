# VGDNN
This package provides the Variants of uncertain significance-guided pretrained deep neural networks (VGDNN) model implementation.

**What we provide:**
- Annotated variant dataset, containing curated features and labels used in this study (Final_annotated.RData).
- Preprocessing scripts for variant filtering, feature preparation, and VUS selection for pretraining (BRCA1_BRCA2_VUS_preprocessing.R, MLH1_MSH2_VUS_preprocessing.R)
- Training code for joint VUS and disease-specific pretraining, followed by gene-specific (BRCA1) or disease-specific (BRCA2, MLH1, MSH2) fine-tuning, as used in this study.

## Copyright
Copyright (c) 2026 Lee, D.-B, Kang, H.-U, and and Hwang, K.-B. All rights reserved. The author retains all rights to this software and documentation

## Contact
If you have any questions, please feel free to contact me: kbhwang@ssu.ac.kr

## Dependencies
VGDNN uses the following software and R packages:
- R (>= 3.6)
- TensorFlow (== 2.5.0)
- Keras (== 2.5.0)
- PRROC (== 1.4)
- caret (== 7.0.1)

## Databases
- ClinVar (GRCh38; clinvar_20220723.vcf.gz)
- Ensembl Variant Effect Predictor version 94
- dbNSFP version 4.2a
- dbSNP build 151
- dbscSNV version 1.1
- gnomAD version 3.1.2
- Korean Variant Archive version 2
- Korean Reference Genome Database phase 2

## Installation
Download the source code and the dataset, then put them in the current working directory.
```sh
wget --no-check-certificate "https://drive.google.com/uc?export=download&id=1lWc2W5TCyj_ujJ81O0dfxx0Sm3LXc0x_" -O VGDNN.tar.gz
```
OR open the save link in any web browser: [https://drive.google.com/file/d/1lWc2W5TCyj_ujJ81O0dfxx0Sm3LXc0x_/view?usp=sharing]

## Usage
1. Extract downloaded folder:
```sh
tar -xvzf VGDNN.tar.gz
cd VGDNN
```

2. Preprocessing for each gene:
```sh
nohup Rscript Preprocessing/BRCA1_BRCA2_VUS_preprocessing.R > logs_BRCA1_BRCA2_VUS_preprocess.txt 2>&1 &
nohup Rscript Preprocessing/MLH1_MSH2_VUS_preprocessing.R  > logs_MLH1_MSH2_VUS_preprocess.txt  2>&1 &
```

3. Run BRCA1 (Gene-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/BRCA1_GS_finetuning_VUS_pretraining/run.sh \
> logs/BRCA1_GS_run.log 2>&1 &
```

4. Run BRCA2 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/BRCA2_DS_finetuning_VUS_pretraining/run.sh \
> logs/BRCA2_DS_run.log 2>&1 &
```

5. Run MLH1 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/MLH1_DS_finetuning_VUS_pretraining/run.sh \
> logs/MLH1_DS_run.log 2>&1 &
```
6. Run MSH2 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/MSH2_DS_finetuning_VUS_pretraining/run.sh \
> logs/MSH2_DS_run.log 2>&1 &
```

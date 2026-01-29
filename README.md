# VGDNN
Copyright (c) 2026 [Author name]. All rights reserved. The author retains all rights to this software and documentation

## Contact
If you have any questions, please feel free to contact me: kbhwang@ssu.ac.kr

## Requirements
- python=3.8.10
- R=3.6.3
- tensorflow=2.5.0
- keras=2.5.0
- PRROC=1.4
- caret=7.0.1

## Quick start guide
1. Download the source code and the dataset, then put them in the current working directory.
```sh
wget --no-check-certificate "https://drive.google.com/uc?export=download&id=1lWc2W5TCyj_ujJ81O0dfxx0Sm3LXc0x_" -O VGDNN.tar.gz
```
OR open the save link in any web browser: [https://drive.google.com/file/d/1lWc2W5TCyj_ujJ81O0dfxx0Sm3LXc0x_/view?usp=sharing]

2. Extract downloaded folder:
```sh
tar -xvzf VGDNN.tar.gz
cd VGDNN
```

3. Preprocessing for each gene:
```sh
nohup Rscript Preprocessing/BRCA1_BRCA2_VUS_preprocessing.R > logs_BRCA1_BRCA2_VUS_preprocess.txt 2>&1 &
nohup Rscript Preprocessing/MLH1_MSH2_VUS_preprocessing.R  > logs_MLH1_MSH2_VUS_preprocess.txt  2>&1 &
```

4. Run BRCA1 (Gene-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/BRCA1_GS_finetuning_VUS_pretraining/run.sh \
> logs/BRCA1_GS_run.log 2>&1 &
```

5. Run BRCA2 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/BRCA2_DS_finetuning_VUS_pretraining/run.sh \
> logs/BRCA2_DS_run.log 2>&1 &
```

6. Run MLH1 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/MLH1_DS_finetuning_VUS_pretraining/run.sh \
> logs/MLH1_DS_run.log 2>&1 &
```
7. Run MSH2 (Disease-specific fine-tuning with VUS & DS pretraining)
```sh
nohup bash VUS_finetuning_pretraining/MSH2_DS_finetuning_VUS_pretraining/run.sh \
> logs/MSH2_DS_run.log 2>&1 &
```

# VGDNN
Copyright (c) 2026 [Author name]. All rights reserved. The author retains all rights to this software and documentation

## Contact
If you have any questions, please feel free to contact me: kbhwang@ssu.ac.kr

## Requirements
- python=3.8.10
- R=4.5.1
- tensorflow=2.5.0
- keras=2.5.0
- numpy=1.19.5
- PRROC=1.4
- caret=7.0.1

## Quick start guide
1. Download the source code and the dataset, then put them in the current working directory.
```sh
wget -c "Download URL"
```
OR open the save link in any web browser: https://xxx.com

2. Extract downloaded folder:
```sh
tar -xvzf Our_Model.tar.gz
cd Our_Model
```

3. Preprocessing for each gene:
```sh
nohup Rscript Preprocessing/BRCA1_BRCA2_VUS_preprocessing.R > logs_BRCA1_BRCA2_VUS_preprocess.txt 2>&1 &
nohup Rscript Preprocessing/MLH1_MSH2_VUS_preprocessing.R  > logs_MLH1_MSH2_VUS_preprocess.txt  2>&1 &
```

5. ㅁㄴㅇ

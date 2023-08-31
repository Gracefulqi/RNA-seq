# RNA-seq test 

## Introduction   
The code in this repository contains the commands and pipeline

## Raw Sequencing Data
Website information


## Requirements

The software packgages needed for the analysis are:

```
Python 3.8.5 
R 
fastp 0.23.4
hisat2
cutadapt 4.4
samtools 1.11
bedtools v2.29.2
```

All of them can be found in the `Packages` directore in this repository.

## 1. RNA-seq pipeline (using RZ1 as an example)
### 1.1 Data quality control
```bash
#PBS N openmpi
##PBS l nodes=1:ppn=12
##PBS j oe
##PBS l walltime=3:00:00
##PBS l mem=10G

fastp -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/rz12d-1_R1.fastq.gz \
      -I /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/rawdata/rz12d-1_R2.fastq.gz \
      -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R1.fastp.fastq.gz \
      -O /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R2.fastp.fastq.gz \
      -w 12 -q 20 -u 20 \
        -j /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1.fastp.json \
        -h /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1.fastp.html
```

### 1.2 Genome mapping
```bash





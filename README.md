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
hisat2 2.1.0
cutadapt 4.4
samtools 1.11
bedtools v2.29.2
StringTie v2.2.1
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
#!/bin/bash
#PBS -N QQtest
#PBS -l nodes=1:ppn=12,mem=10G
#PBS -q batch

hisat2 -p 12
       -x /public/home/zhangqq/Tair10_genome/hisat2_index/TAIR10 \
       --summary-file /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.summary \
       --dta-cufflinks -1 /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R1.fastp.fastq.gz \
       -2 /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/fil/rz12d-1_R2.fastp.fastq.gz | \
samtools view -ShuF 4 -q 20 -f 2 -@ 8 - | \
samtools sort -@ 8 -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.rep1.sorted.bam -
```

### 1.3 Convert the BAM format to bigwig format for tract visualization
```bash
samtools index /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.rep1.sorted.bam \
bamCoverage --bam /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.rep1.sorted.bam \
            -o rz1.deeptools.bw \
            --binSize 10 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 119481543 #拟南芥染色体基因组大小
```

### 1.4 Get the gene expression matrix
### 1.4.1 Assemble the transcripts
```bash
#!/bin/bash
#PBS -N QQtest
#PBS -l nodes=1:ppn=12,mem=10G
#PBS -q batch
#PBS -j oe

stringtie /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/map/rz1.rep1.sorted.bam \ # 此bam是samtools sort处理后的文件
          -G /public/home/zhangqq/Tair10_genome/TAIR10.gff3 \ #参考基因组注释文件
          -l rz1 -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/rz1.transcripts.stringtie.gtf \
          -p 12 -e 
```
### 1.4.2 Merge the transcript samples (处理多个生物学重复样本时需要合并)
```bash
vim mergelist.txt #需要包含之前output.gtf文件的路径
/public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/rz1.transcripts.stringtie.gtf #txt文件内容

stringtie --merge -G /public/home/zhangqq/Tair10_genome/TAIR10.gff3 \
          -F 0.1 -T 0.1 -i -o /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/rz1.stringtie_merged.gtf \ 
          mergelist.txt
``` 
### 1.4.3 Quantitative analysis
'''bash
vim sample_list.txt #需要包含样本名和定量的gtf文件的路径
rz1      /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/rz1.transcripts.stringtie.gtf #txt文件内容，文件名和gtf文件之间用TAB键隔开

source activate python2
python /public/home/zhangqq/software/stringtie-2.2.1/prepDE.py \ #使用python的prepDE.py命令(prepDE.py在stringtie下面，写上prepDE.py的绝对路径)
       -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/sample_list.txt \
       -g /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/gene_count_matrix.csv \
       -t /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/transcript_count_matrix.csv 

   或者用另外一种方式(prepDE.py在base环境下也可以使用)
prepDE.py -i /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/sample_list.txt \
          -g /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/gene_count_matrix.csv \
          -t /public/home/zhangqq/RNA-seq_Col_rz1_FangYJ/gene_expression/transcript_count_matrix.csv 
```





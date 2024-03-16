# Detecting-m5C-sites-using-Nanopore-direct-RNA-sequencing
# Before start
Use the command ```pip install -r requirements.txt``` to download the m5Cnet environment.
# Step I. Basecalling
We must using basecalling software to get the fastq files from fast5 files
```
$ guppy_basecaller -i ./fast5_dir \
-c ./guppy/data/rna_r9.4.1_70bps_hac.cfg \
-s ./guppy_result \
-x "cuda:0" --fast5_out > guppy.log
```
# Step II. Alignment
Alignment basecalled fastqs to reference sequences using minimap2, and sort the bam file by samtools
```
$ minimap2 -ax -splice -uf -k14 ref.fa pass.fastq > pass.sam
$ samtools view -bS pass.sam > pass.bam
$ samtools sort pass.bam -o pass.sort.bam
```
# Step III. Nanopolish index and eventalign
Using Nanopolish to align the electrical signals produced by nanopore sequencing with the reference sequence 
```
$ nanopolish index -d fast5/ pass.fastq
$ nanopolish eventalign --reads pass.fastq --bam pass.sort.bam --genome ref.fa --scale-events --signal index --summary summary.txt > pass.eventalgin.txt
```
# step IV. Generating the features of specific sites
Using dataprep.py to generate the features, the result file included 'data.json', 'data.infor' and 'data.log'. Using command ```$ python dataprep.py --help``` can get the information of usage. This step may use lots of time.
```
$ python data.prep --n_neighbor=1\
--min_reads=20 \
-i pass.eventalgin.txt \
-o dataprep_output_dir_path \
-t 8
```
# step V. Predict m5C sites and get methylation ratio and probability
Using inferece.py to predict the m5Csites. Using command ```$ python inference.py``` can get the information of usage.
```
$ python inference.py -i dataprep_output_dir_path \
-o inference_result_dir_path
```
# Reference:
Hendra C, Pratanwanich PN, Wan YK, Goh WSS, Thiery A, Goke J. Detection of m6A from direct RNA sequencing using a multiple instance learning framework. Nature methods 2022, 19(12): 1590-1598.

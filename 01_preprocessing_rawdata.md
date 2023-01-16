#### Rename raw data files based on sample names
```bash
# Kronja et al., 2014
mv SRR1039771.fastq.gz activatedegg_chx_rpf.fq.gz
mv SRR1039768.fastq.gz activatedegg_nochx_rpf.fq.gz

# Patraquim et al., 2022
mv SRR19387512.fastq.gz embryo_early_rpf.fq.gz
mv SRR19387513.fastq.gz embryo_early_rna.fq.gz

zcat SRR11432001.fastq.gz SRR11432002.fastq.gz | gzip -c > embryo_mid_rna.fq.gz
mv SRR19387520.fastq.gz embryo_mid_rpf.fq.gz

mv SRR19387519.fastq.gz embryo_late_rpf.fq.gz
mv SRR19387518.fastq.gz embryo_late_rna.fq.gz

# Zhang et al., 2018
mv SRR3031124.fastq.gz embryo_2h_rpf.fq.gz
mv SRR5075630.fastq.gz embryo_2h_rna.fq.gz

mv SRR3031125.fastq.gz embryo_6h_rpf.fq.gz
mv SRR5075635.fastq.gz embryo_6h_rna.fq.gz

mv SRR3031126.fastq.gz embryo_12h_rpf.fq.gz
mv SRR5075633.fastq.gz embryo_12h_rna.fq.gz

mv SRR3031127.fastq.gz embryo_24h_rpf.fq.gz
mv SRR5075631.fastq.gz embryo_24h_rna.fq.gz
```

#### trim adaptors
```bash
cutadapt -a TCGTATGCCGTCTTCTGCTTG -u8 -j8 --trim-n -m 18 -o activatedegg_chx_rpf.trim.fa.gz activatedegg_chx_rpf.fq.gz
cutadapt -a TCGTATGCCGTCTTCTGCTTG -u8 -j8 --trim-n -m 18 -o activatedegg_nochx_rpf.trim.fa.gz activatedegg_nochx_rpf.fq.gz

cutadapt -a AGATCGGAAGAGCACACGTC -j8 --trim-n -m 18 -o embryo_mid_rna.trim.fq.gz embryo_mid_rna.fq.gz >embryo_mid_rna.trim.log

cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_2h_rna.trim.fq.gz embryo_2h_rna.fq.gz >embryo_2h_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_2h_rpf.trim.fq.gz embryo_2h_rpf.fq.gz >embryo_2h_rpf.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_6h_rna.trim.fq.gz embryo_6h_rna.fq.gz >embryo_6h_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_6h_rpf.trim.fq.gz embryo_6h_rpf.fq.gz >embryo_6h_rpf.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_12h_rna.trim.fq.gz embryo_12h_rna.fq.gz >embryo_12h_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_12h_rpf.trim.fq.gz embryo_12h_rpf.fq.gz >embryo_12h_rpf.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_24h_rna.trim.fq.gz embryo_24h_rna.fq.gz >embryo_24h_rna.trim.log
cutadapt -a TGGAATTCTCGGGTGCCAAGG -j8 --trim-n -m 18 -o embryo_24h_rpf.trim.fq.gz embryo_24h_rpf.fq.gz >embryo_24h_rpf.trim.log
```

#### filter rRNAs and tRNAs
rRNA and tRNA sequences were downloaded from FlyBase and combined to a single fasta.
```r
library(Biostrings)
trna <- readDNAStringSet('dmel-all-tRNA-r6.47.fasta.gz')
rrna <- readDNAStringSet('dmel-all-miscRNA-r6.47.fasta.gz')
rrna <- rrna[grepl(' type=rRNA;', names(rrna))]
misc <- c(trna, rrna)
writeXStringSet(misc, 'rRNA_tRNA_combined.fa')
```

Build bowtie2 index for rRNAs and tRNAs
```bash
bowtie2-build rRNA_tRNA_combined.fa bt2_index_rtRNA
```

Remove reads mapped to rRNA and tRNAs
```bash
# for files that do not require adaptor trimming
for i in *.fq.gz; do echo "bowtie2 -p8 --local --un-gz ${i%%.*}.clean.fq.gz -x bt2_index_rtRNA -U $i > /dev/null 2>${i%%.*}.clean.log"; done

# for files that require adaptor trimming
for i in *trim.fq.gz; do echo "bowtie2 -p8 --local --un-gz ${i%%.*}.clean.fq.gz -x bt2_index_rtRNA -U $i > /dev/null 2>${i%%.*}.clean.log"; done
```

#### mapping
```bash
for i in *_rna.clean.fq.gz; do echo "STAR --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir /nfs_data/database/ref_genomes/Dmel_em52/STAR_genome --readFilesIn $i --readFilesCommand zcat --outFileNamePrefix ${i%%.clean.fq.gz}_ --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd"; done
```

#### Quantify transcript expression levels using salmon
This step is helpful for choosing the representative transcript isoform for each gene based on highest expression. Activated eggs are sampled 0-2h after egg laying, which is similar to the sampling stage of our previous 0-2h embryos data.
```bash
for i in *_rna.clean.fq.gz; do echo "salmon quant -p8 --seqBias --gcBias --posBias -l A -i /nfs_data/database/ref_genomes/Dmel_em52/salmon_cdna_ncrna -r $i -o ${i%%.clean.fq.gz}"; done

ln -s  embryo_2h_rna activatedegg_chx_rna
ln -s  embryo_2h_rna activatedegg_nochx_rna
```

#### Extract features for developing the PSite model
```bash
for i in *rpf_Aligned.toTranscriptome.out.bam; do echo "python ../src/extract_features.py -n3 /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz $i > ${i%%_Aligned*}.features.3nt.tsv"; done
```

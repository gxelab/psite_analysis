### run PSite pipeline and compute coverage
```bash
# train
psite train -i /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz ../data/embryo_2h_rpf_Aligned.toTranscriptome.out.bam embryo_2h_rpf.psite ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo

# predict
psite predict -i ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa ../data/embryo_2h_rpf_Aligned.sortedByCoord.out.bam embryo_2h_rpf.psite.gbt.pickle embryo_2h_rpf_Aligned.sortedByCoord.psite.bam

psite predict -i ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa ../data/embryo_2h_rpf_Aligned.sortedByCoord.out.bam embryo_2h_rpf.psite.dam.pickle embryo_2h_rpf_Aligned.sortedByCoord.peak.bam

python ../src/setp.py -l 29 -u 35 -n 12 ../data/embryo_2h_rpf_Aligned.sortedByCoord.out.bam embryo_2h_rpf_Aligned.sortedByCoord.fixed.bam

# calculate coverage
psite coverage embryo_2h_rpf_Aligned.sortedByCoord.psite.bam embryo_2h_rpf_Aligned.sortedByCoord.psite
psite coverage embryo_2h_rpf_Aligned.sortedByCoord.peak.bam embryo_2h_rpf_Aligned.sortedByCoord.peak
psite coverage embryo_2h_rpf_Aligned.sortedByCoord.fixed.bam embryo_2h_rpf_Aligned.sortedByCoord.fixed
```

### metagene profile analysis
```bash
python src/gtf.py convert2bed -g Drosophila_melanogaster.BDGP6.32.52.gtf -t exon >Drosophila_melanogaster.BDGP6.32.52.tx.bed12

python ../src/bigwig_covpn.py embryo_2h_rpf_Aligned.sortedByCoord.psite_fw.bw embryo_2h_rpf_Aligned.sortedByCoord.psite_rc.bw ../Drosophila_melanogaster.BDGP6.32.52.tx.bed12 >embryo_2h_rpf_Aligned.sortedByCoord.psite.txcov.tsv

python ../src/bigwig_covpn.py embryo_2h_rpf_Aligned.sortedByCoord.peak_fw.bw embryo_2h_rpf_Aligned.sortedByCoord.peak_rc.bw ../Drosophila_melanogaster.BDGP6.32.52.tx.bed12 >embryo_2h_rpf_Aligned.sortedByCoord.peak.txcov.tsv

python ../src/bigwig_covpn.py embryo_2h_rpf_Aligned.sortedByCoord.fixed_fw.bw embryo_2h_rpf_Aligned.sortedByCoord.fixed_rc.bw ../Drosophila_melanogaster.BDGP6.32.52.tx.bed12 >embryo_2h_rpf_Aligned.sortedByCoord.fixed.txcov.tsv
```

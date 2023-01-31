# Code for PSite manuscript

#### Public data analyzed in the manuscript
Ribo-Seq and RNA-Seq datasets:
- Drosophila Activated eggs, CHX/noCHX, RNase I (GSE52799)
  > Kronja, I., et al. Widespread changes in the posttranscriptional landscape at the Drosophila oocyte-to-embryo transition. Cell Rep 2014;7(5):1495-1508. [PubMed](https://pubmed.ncbi.nlm.nih.gov/24882012/)

- Drosophila embryos early/mid/late, RNase I (GSE204739)
  > Patraquim, P., et al. Translation and natural selection of micropeptides from long non-canonical RNAs. Nat Commun 2022;13(1):6515. [PubMed](https://pubmed.ncbi.nlm.nih.gov/36316320/)

- Drosophila embryos 0-2h/2-6h/6-12h/12-24h, MNase (SRP067542)
  > Zhang, H., et al. Genome-wide maps of ribosomal occupancy provide insights into adaptive evolution and regulatory roles of uORFs during Drosophila development. PLoS Biol 2018;16(7):e2003903. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30028832/)


Sequences and annotations:
- reference genome and annotations:
  - http://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
  - http://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz
  - http://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.32.ncrna.fa.gz
  - http://ftp.ensemblgenomes.org/pub/metazoa/release-52/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.52.gtf.gz
- rRNAs and tRNAs
  - http://ftp.flybase.net/releases/FB2022_04/dmel_r6.47/


#### Software
- [cutadapt](https://github.com/marcelm/cutadapt) v3.7
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.4.5
- [STAR](https://github.com/alexdobin/STAR) 2.7.4a
- [salmon](https://github.com/COMBINE-lab/salmon) v1.8.0
- [PSite](https://github.com/gxelab/psite) v0.1.1
- [RiboCode](https://github.com/xryanglab/RiboCode) v1.2.14
- [Ribo-TISH](https://github.com/zhpn1024/ribotish) dev (2022-08-11)
- python v3.9.10
- R v4.1.3
- Module [ncorf_classifer3.py](https://github.com/gxelab/orftools/blob/main/ncorf_classifier3.py)
- Module [gtf.py](https://github.com/mt1022/GPP/blob/main/gpp/gtf.py)
- Module [bigwig_covpn.py](https://github.com/gxelab/scripts/blob/main/bigwig_covpn.py)

#### Scripts
- 01: raw data processing;
- 02: run the PSite (train + prediction);
- 03: predict translated ORFs with Ribo-TISH and RiboCode using the default pipeline;
- 04: predict translated ORFs with Ribo-TISH and RiboCode using the PSite pipeline;
- 05: metagene profile analysis of P-site coverage around CDS start codons;
- 06：code to generate plots shown in figure 1;

#### Ribo-TISH
```bash
# index bam
for i in *_rpf_Aligned.sortedByCoord.out.bam; do echo "samtools index $i"; done

# estimate offsets
for i in *_rpf_Aligned.sortedByCoord.out.bam; do echo "ribotish quality -p8 -b $i -g ../Drosophila_melanogaster.BDGP6.32.52.gtf"; done

# predict
for i in *_rpf_Aligned.sortedByCoord.out.bam; do echo "ribotish predict -p8 --alt --altcodons CTG,TTG,GTG --framebest -b $i -g ../Drosophila_melanogaster.BDGP6.32.52.gtf -f ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa -o ../orfpred/${i%%_Aligned*}.ribotish_pred.txt --allresult ../orfpred/${i%%_Aligned*}.ribotish_allres.txt"; done

# unified classification of predicted ORFs
for i in *ribotish_pred.txt; do echo "python ../src/ncorf_classifier3.py -k -m ribotish -p ../allorfs/default/${i%%_pred.txt} $i ../Drosophila_melanogaster.BDGP6.32.52.gtf ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo"; done
```

#### RiboCode
```bash
conda activate ribocode
# prepare annotations
ln -s /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa
ln -s /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.52.gtf
prepare_transcripts -g Drosophila_melanogaster.BDGP6.32.52.gtf -f Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa -o ribocode_annot

# estimate P-site offsets by meta-gene profile analysis
metaplots -f0_percent 0.6 -a ../ribocode_annot -r activatedegg_chx_rpf_Aligned.toTranscriptome.out.bam -o activatedegg_chx_rpf
metaplots -f0_percent 0.6 -a ../ribocode_annot -r activatedegg_nochx_rpf_Aligned.toTranscriptome.out.bam -o activatedegg_nochx_rpf

# prediction (default canonical start is ATG and is set by "-s")
RiboCode -a ../ribocode_annot -c activatedegg_chx_rpf_pre_config.txt -l no -g -A CTG,GTG,TTG -o ../orfpred/activatedegg_chx_rpf
RiboCode -a ../ribocode_annot -c activatedegg_nochx_rpf_pre_config.txt -l no -g -A CTG,GTG,TTG -o ../orfpred/activatedegg_nochx_rpf

# unified classification of predicted ORFs
python ../src/ncorf_classifier3.py -k -m ribocode -p ../allorfs/default/activatedegg_chx_rpf.ribocode activatedegg_chx_rpf.txt ../Drosophila_melanogaster.BDGP6.32.52.gtf ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo
python ../src/ncorf_classifier3.py -k -m ribocode -p ../allorfs/default/activatedegg_nochx_rpf.ribocode activatedegg_nochx_rpf.txt ../Drosophila_melanogaster.BDGP6.32.52.gtf ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo
```
NB: other samples have too few frame 0 reads (less than 45%) for `RiboCode` to run.

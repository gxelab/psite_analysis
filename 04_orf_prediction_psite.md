#### Ribo-TISH
```bash
# PSite predict doesn't change alignment order, so that no sorting needed before index
for i in *sortedByCoord.*.tag.bam; do echo "samtools index $i"; done

# predict
for i in *sortedByCoord.om10.tag.bam; do echo "python ../src/ribotish_predict.py -p4 --alt --altcodons CTG,TTG,GTG --framebest -b $i -g ../Drosophila_melanogaster.BDGP6.32.52.gtf -f ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa -o ../orfpred_psite/${i%%_Aligned*}.om10.ribotish_pred.txt --allresult ../orfpred_psite/${i%%_Aligned*}.om10.ribotish_allresult.txt"; done

for i in *sortedByCoord.psite.tag.bam; do echo "python ../src/ribotish_predict.py -p4 --alt --altcodons CTG,TTG,GTG --framebest -b $i -g ../Drosophila_melanogaster.BDGP6.32.52.gtf -f ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa -o ../orfpred_psite/${i%%_Aligned*}.psite.ribotish_pred.txt --allresult ../orfpred_psite/${i%%_Aligned*}.psite.ribotish_allresult.txt"; done
```

#### RiboCode
```bash
# generate pseudo config files
# note: read length filtering is done in prediction step. No need to care it here.
for i in *.psite.tag.bam; do python ../src/ribocode_config.py $i ${i%%_Aligned*}.psite.log > ${i}.pre_config; done
for i in *.om10.tag.bam; do python ../src/ribocode_config.py $i ${i%%_Aligned*}.om10.log > ${i}.pre_config; done

# predict
for i in *om10.tag.bam; do echo "python ../src/RiboCode_psite.py -a ../ribocode_annot -c ${i}.pre_config -l no -g -A CTG,GTG,TTG -o ../orfpred_psite/${i%%_Aligned*}.om10"; done
for i in *psite.tag.bam; do echo "python ..//src/RiboCode_psite.py -a ../ribocode_annot -c ${i}.pre_config -l no -g -A CTG,GTG,TTG -o ../orfpred_psite/${i%%_Aligned*}.psite"; done
```

#### ORF classification
```bash
for i in *_rpf.{psite,om10}.txt; do echo "python ../src/ncorf_classifier3.py -k -m ribocode -p ../allorfs/psite/${i%%.txt}.ribocode $i ../Drosophila_melanogaster.BDGP6.32.52.gtf ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo"; done
for i in *rpf.{om10,psite}.ribotish_pred.txt; do echo "python ../src/ncorf_classifier3.py -k -m ribotish -p ../allorfs/psite/${i%%_pred.txt} $i ../Drosophila_melanogaster.BDGP6.32.52.gtf ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo"; done
```

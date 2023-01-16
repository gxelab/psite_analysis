#### Prepare transcript metadata file
```bash
src/extract_txinfo_ensembl.R Drosophila_melanogaster.BDGP6.32.52.gtf.gz Drosophila_melanogaster.BDGP6.32.52.gtf.info
```

#### Train
Train with parameter `--offset_min 10`.
```bash
for i in *rpf_Aligned.toTranscriptome.out.bam; do echo psite train -i --offset_min 10 /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz $i ${i%%_Aligned*}.om10 ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo; done
```

Train with parameter `--offset_min 11`.
```bash
for i in *rpf_Aligned.toTranscriptome.out.bam; do echo psite train -i /nfs_data/database/ref_genomes/Dmel_em52/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz $i ${i%%_Aligned*}.psite ../Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo; done
```

#### Predict
Get  bam  with "PS" tag using the `predict` submodule and P-site only bams using the `pbam` submodule.
```bash
for i in *rpf_Aligned.toTranscriptome.out.bam; do echo psite predict -i ../Drosophila_melanogaster.BDGP6.32.cdna_plus_ncrna.fa $i ${i%%_Aligned*}.om10.gbt.pickle ${i%%.out.bam}.om10.tag.bam; done
for i in *rpf_Aligned.sortedByCoord.out.bam; do echo psite predict -i ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa $i ${i%%_Aligned*}.om10.gbt.pickle ${i%%.out.bam}.om10.tag.bam; done

for i in *rpf_Aligned.toTranscriptome.out.bam; do echo psite predict -i ../Drosophila_melanogaster.BDGP6.32.cdna_plus_ncrna.fa $i ${i%%_Aligned*}.psite.gbt.pickle ${i%%.out.bam}.psite.tag.bam; done
for i in *rpf_Aligned.sortedByCoord.out.bam; do echo psite predict -i ../Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa $i ${i%%_Aligned*}.psite.gbt.pickle ${i%%.out.bam}.psite.tag.bam; done
```

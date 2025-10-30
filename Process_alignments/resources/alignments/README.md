Tree (mashtree v1.4.6):

```
mashtree --numcpus 12 *.fasta > 8-t2t-apes-2023v2.dnd
```
Reroot the tree with siamang as the outgroup.

Alignment (cactus v2.7.1):

Make a SLURM config (bigger lastz chunks)

```
cp /private/groups/cgl/cactus/cactus-bin-v2.7.1/src/cactus/cactus_progressive_config.xml ./config-slurm.xml
sed -i config-slurm.xml -e 's/blast chunkSize="30000000"/blast chunkSize="90000000"/g'
sed -i config-slurm.xml -e 's/dechunkBatchSize="1000"/dechunkBatchSize="200"/g'
```
Make the alignment using this input
```
(((GCA_028885655.2:0.0017500000000000016,GCA_028885625.2:0.001729999999999999):0.014950000000000001,(GCA_029281585.2:0.00877,((GCA_029289425.2:0.0019999999999999983,GCA_028858775.2:0.0022900000000000004):0.0043300000000000005,(hs1:5.000000000000004E-4,hg38:5.000000000000004E-4):0.005989999999999999):0.0014300000000000007):0.0073599999999999985):0.011345,GCA_028878055.2:0.011345000000000003);

hs1                  https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
hg38                 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
GCF_029281585.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/281/585/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCA_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz
GCF_028858775.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/858/775/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCA_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz
GCF_029289425.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz
GCF_028885655.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/655/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCA_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz
GCF_028885625.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/885/625/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCA_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz
GCF_028878055.2 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/878/055/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri/GCA_028878055.2_NHGRI_mSymSyn1-v2.0_pri_genomic.fna.gz

```
The accessions are:
```
GCF_029281585.2 Gorilla_gorilla
GCF_028858775.2 Pan_troglodytes
GCF_029289425.2 Pan_paniscus
GCF_028885655.2 Pongo_abelii
GCF_028885625.2 Pongo_pygmaeus
GCF_028878055.2 Symphalangus_syndactylus
```
Cactus command

```
TOIL_SLURM_ARGS="--partition=long --time=8000" cactus ./js-8apes ./8-t2t-apes-2023v2.seqfile ./8-t2t-apes-2023v2.hal --batchSystem slurm --caching false --consCores 64 --configFile ./config-slurm.xml --logFile 8-t2t-apes-2023v2.hal.log --batchLogsDir batch-logs-8apes --coordinationDir /data/tmp
```

MAF for each primary assembly

for i in hs1 hg38 GCA_028858775.2 GCA_028885655.2 GCA_028885625.2 GCA_028878055.2 GCA_029281585.2 GCA_029289425.2; do TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-hal2maf ./js_hal2maf8 ./8-t2t-apes-2024v1.hal ./8-t2t-apes-2024v1.${i}.maf.gz --filterGapCausingDupes --refGenome $i --chunkSize 500000 --batchCores 64 --noAncestors --batchCount 16  --batchSystem slurm --caching false --logFile ./8-t2t-apes-2024v1.${i}.gz.log --batchLogsDir batch-logs-8apes --coordinationDir /data/tmp ;done

Then make a bigmaf (below is not run on slurm)

for i in hs1 hg38 GCA_028858775.2 GCA_028885655.2 GCA_028885625.2 GCA_028878055.2 GCA_029281585.2 GCA_029289425.2; do zcat ./8-t2t-apes-2024v1.${i}.maf.gz | mafDuplicateFilter -km - | bgzip > ./8-t2t-apes-2024v1.${i}.single-copy.maf.gz; done

for i in hs1 hg38 GCA_028858775.2 GCA_028885655.2 GCA_028885625.2 GCA_028878055.2 GCA_029281585.2 GCA_029289425.2; do cactus-maf2bigmaf ./js-bb ./8-t2t-apes-2024v1.${i}.single-copy.maf.gz ./8-t2t-apes-2024v1.${i}.bigmaf.bb --refGenome $i --halFile ./8-t2t-apes-2024v1.hal ; done

for i in hs1 hg38 GCA_028858775.2 GCA_028885655.2 GCA_028885625.2 GCA_028878055.2 GCA_029281585.2 GCA_029289425.2; do rm ./8-t2t-apes-2024v1.${i}.single-copy.maf.gz ; done



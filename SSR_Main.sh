#!/bin/bash

### run the rscript independently, make sure it stopped before filterSNPGene 
### Should make it parallele at least for different key words
### Should make it parallele in general
### "immune" is the key word hard-coded here
#Rscript findGO.r immune & 

### have user input the header of the map and ped file
### have user input the identifiers 
plinkHeader="mapped_"
plinkTail='Blocks' 
identifiers=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

for id0 in ${identifiers[@]} 
do 
   plink --file $plinkHeader$id0 --blocks no-pheno-req --out $id0$plinkTail 
done  


### RUN ANNOVAR
### Require downloading annovar package 
### annovar parameters 
ANNOVAR_PATH=`pwd`'/annovar/' 

inputTail=".blocks" 
downdbFlag=False	# user-defined  
MAX_GENE_THREAD=6	# user-defined
THREAD=$MAX_GENE_THREAD 

snpInfoHeader='snp'	# user-defined? 
snpInfoTail='.avinput' 
annotateHeader=$snpInfoHeader 
annotateTail='.variant_function' 

for id0 in ${identifiers[@]} 
do
   snpList0=$annotateHeader$id0'.txt'
   cat $id0$plinkTail$inputTail | cut -d' ' -f2 > $snpList0 
   clist=`echo $clist$annotateHeader$id0" "`
done

clist=($clist) 

### download annotation databases if nneded 
if $downdbFlag
then
   perl $ANNOVAR_PATH'annotate_variation.pl' -downdb -buildver hg19 -webfrom annovar refGene $ANNOVAR_PATH'humandb/'
   perl $ANNOVAR_PATH'annotate_variation.pl' -downdb -buildver hg19 -webfrom annovar snp138 $ANNOVAR_PATH'humandb/' 
fi 

### annotate SNPs
for c0 in ${clist[@]} 
do
   perl $ANNOVAR_PATH'convert2annovar.pl' -format rsid $c0'.txt' -dbsnpfile $ANNOVAR_PATH'humandb/hg19_snp138.txt' > $c0$snpInfoTail    
   perl $ANNOVAR_PATH'annotate_variation.pl' --maxgenethread $MAX_GENE_THREAD --thread $THREAD -out $c0 -build hg19 $c0$snpInfoTail $ANNOVAR_PATH'humandb/' --geneanno -dbtype refGene
done 

### prepare parameters for the filterSNPGene.py

filePath=`pwd` 
mapHeader=$plinkHeader
mapTail=".map" 
pedHeader=$plinkHeader 
pedTail=".ped"  

ids=$(printf "%s," "${identifiers[@]}")  

spark-submit filterSNPGene.py $snpInfoHeader $snpInfoTail $annotateHeader $annotateTail $mapHeader $mapTail $pedHeader $pedTail $ids `pwd` "test"  

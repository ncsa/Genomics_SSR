# Genomics_SSR
Coarse Search Space Reduction Workflow

# Background
Many phenotypes are associated with multiple interacting genomic variants. To model these phenotypes, so many combinations of variant/varint-pairs would be tested that not only it become very inefficient, but also that the multiple-hypothesis correction made the modeling impossible. Therefore, combining biological information with modeling to reduce search space is essential.

# Workflow
## 0. File format 
Any input files should be compatible with PLINK format, so that we could use the PLINK haplotype module to do the haplotype. We choose the ".ped" format with a ".map" file.

## 1. Haplotyping
We use the PLINKv1.90 \[1\] haplotype module to do the haplotyping. The output would include two files, one ".blocks" showing the haplotyped SNPs, each haplotype block in one line. The second one is ".blocks.det", which provides more detailed information about chromosome, bp, length, and number of SNPs.

> plink PLINK --file $FILENAME --blocks no-pheno-req --out $FILENAME"blocks"

In our case, the input data is already haplotyped, therefore there are not many SNPs in one haplotype "block", nor the blocks are very large. Therefore, we use annovar to annotate single SNPs to genes, instead of using R GenomicRanges packages to annotate haplotype blocks to genes. Both annovar and GenomicRanges could annotate SNPs to more than coding regions, but again for simplicity of this prototype, we only use coding regiosn of genes.

## 2. Annotation
We use the package annovar\[2\] (http://annovar.openbioinformatics.org/en/latest/) to annotate SNPs to genes.

## 3. Filtering
First, user would provide key word searching genes associated with the phenotype of interest. The package "rols"\[3\] (http://lgatto.github.com/rols/) from bioconductor finds the list of GO_TERM_ID's associated with the key word. Then the workflow would use the package "biomaRt"\[4\]  select genes that are annotated to the key words. SNPs associated with the selected genes are take into consideration. For this coarse version, only Gene Ontology is supported. We filter the SNPs by pySpark.

## 4. Output
Then, the filtered SNPs are output as a genotype table ready for EpiQuant to use.

# Usage
* **pre-requisits**

spark - 2.4.0

python - 3.7

R - 3.5.1

  bioconductor libraries: rols (2.10.1), biomaRt (2.38.0)
  
PLINK - v1.90b6.7
annovar - 2018Apr16
 

* **Input**

".ped" file for user-measured genotype files.

"keyword" for gene ontology annotation.


The program takes in the keyword and the ped file name and automatically execute the workflow. It might take a long time before the final filtered .epiq file is generated. The .epiq file would be used in GenomicsEpiquant package for 

# Notes
This is a coarse model where only SNPs directly annotated to the selected genes are considered, and coding/non-coding region are not separated. We are adding a ranking procedure to include more genes interacting with select genes and SNPs associated with those genes ranked. Later, the regulatory information would also be incorporated into the workflow.

# References
\[1\] Purcell, Shaun, et al. "PLINK: a tool set for whole-genome association and population-based linkage analyses." The American journal of human genetics 81.3 (2007): 559-575.
\[2\]Wang, Kai, Mingyao Li, and Hakon Hakonarson. "ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data." Nucleic acids research 38.16 (2010): e164-e164.
\[3\]Gatto, Laurent. "rols: an R interface to the Ontology Lookup Service." (2013).
\[4\]Durinck, Steffen, et al. "Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt." Nature protocols 4.8 (2009): 1184.

### add require package installation statement

library("biomaRt")
library("rols")

# read argument
args<-commandArgs(trailingOnly = TRUE)
keyword=args[1]

# starting rols query to get the list of gene ontologies

  # build query, limited to GO
query0<-allRows(OlsSearch(q=keyword, ontology="GO"))

  # search
q<-olsSearch(query0)

  # slow, convert to a list of terms for bioMart Search;
  # GO.db package might be faster because it uses local search
  # rols uses online search.
qtrms<-as(q,"Terms")
queryList<-names(termOntology(unique(qtrms)))

# start biomaRt annotation to genes
### this part needs to be modified, currently hard-coded for human gene dataset
### biomart only have GRCh38.p12 for human genome, but we only need gene names, it's ok.
### need to change the hard-coded part for species
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
genelist = getBM(attributes=c('chromosome_name','external_gene_name'),filters='go',values=queryList, mart=ensembl)

#genelist=unique(genelist)
write.table(genelist,"./selectedGenes.txt",row.names=FALSE,quote=FALSE,sep="\t")

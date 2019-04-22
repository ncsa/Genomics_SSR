from pyspark import SparkContext
from pyspark.sql import SparkSession, SQLContext
from pyspark.sql.functions import col 
from functools import reduce 
import pandas as pd 
import sys

### filter SNP
### annotateFile is the annovar output $FILE.variant_function
### geneSel is a list of genes in "external_gene_name"
### currently only work for exons 

def fileNameCreator(identifiers):
    files={}
    counts=0
    for id0 in identifiers:
        files['annotateFile'+str(counts)]=annotateHeader+id0+annotateTail
        files['mapFile'+str(counts)]=mapHeader+id0+mapTail
        files['snpInfoFile'+str(counts)]=snpInfoHeader+id0+snpInfoTail
        files['pedFile'+str(counts)]=pedHeader+id0+pedTail
        counts+=1 
    return {'files':files,'fileNum':counts} 

def filterSNP(annotateFile,geneSel):
    temp_var = annotateFile.map(lambda k: k.split('\t')) 
    df = temp_var.toDF()
    geneSNP = df.filter(df['_1'].startswith('exonic')).select('_2','_8') 
    filteredGeneSNP = geneSNP.where(col('_2').isin(geneSel)) 
    snpList=[row._8 for row in filteredGeneSNP.select(filteredGeneSNP.columns[1]).dropDuplicates().collect()] 
    return snpList 

def getSNPInfo(mapFile,snpInfoFile,snpList): 
    rsidColMap = mapFile.zipWithIndex().map(lambda k: (k[0].split()+[k[1]]))\
            .map(lambda k: (k[1],k[3])).filter(lambda k: k[0] in snpList) 
    minorMap = snpInfoFile.map(lambda k: k.split('\t'))\
            .map(lambda k: (k[5],k[4]))\
            .filter(lambda k: k[0] in snpList) 
    colMinorMap=rsidColMap.join(minorMap).map(lambda x: x[1]).collectAsMap() 
    return colMinorMap 

def filterPed(pedFile,colMinorMap):
    df = pedFile.map(lambda k: k.split(' ')).map(lambda k: [k[1]]+\
            [2 if(k[6+2*ii]==colMinorMap[ii])\
            else 1 if(k[6+2*ii+1]==colMinorMap[ii])\
            else 0\
            for ii in colMinorMap.keys()]).toDF()
    return df

def joinDF(df1,df2):
    df1 = df1.alias('df1')
    df2 = df2.alias('df2') 
    joined_df = df1.join(df2, col('df1.headerLine')==col('df2.headerLine'),'left')\
            .select(['df1.*']+[col('df2.'+xx) for xx in df2.columns[1:]])
    return joined_df 

def resultAssemble(filePath,fileCollection,geneSel):
    snpList0=[]
    fileNum=fileCollection['fileNum']
    files=fileCollection['files']
    for ii in range(fileNum):
        annotateFile=sc.textFile('file:///'+filePath+'/'+files['annotateFile'+str(ii)]) 
        snpList=filterSNP(annotateFile,geneSel)
        snpList0 = snpList0+snpList

        mapFile=sc.textFile('file:///'+filePath+'/'+files['mapFile'+str(ii)]) 
        snpInfoFile=sc.textFile('file:///'+filePath+'/'+files['snpInfoFile'+str(ii)])
        colMinorMap=getSNPInfo(mapFile,snpInfoFile,snpList)

        pedFile=sc.textFile('file:///'+filePath+'/'+files['pedFile'+str(ii)])
        df0=filterPed(pedFile,colMinorMap)
        oldColumns = df0.schema.names
        newColumns = ["headerLine"] + snpList
        df0 = reduce(lambda data, idx: data.withColumnRenamed(oldColumns[idx],newColumns[idx]),range(len(oldColumns)), df0)

        if ii == 0:
            df = df0
        else:
            df = joinDF(df0,df)

    pdDfOut=df.toPandas().set_index("headerLine").transpose() 
    return pdDfOut 

if __name__ == '__main__': 
    ### take in arguments: 
    ### file names
    snpInfoHeader=sys.argv[1] 
    snpInfoTail=sys.argv[2]
    annotateHeader=sys.argv[3]
    annotateTail=sys.argv[4]
    mapHeader=sys.argv[5]
    mapTail=sys.argv[6] 
    pedHeader=sys.argv[7]
    pedTail=sys.argv[8] 
    identifiers=sys.argv[9]
    identifiers=identifiers[:len(identifiers)-1].split(',') 
    filePath=sys.argv[10] 
    fonm=sys.argv[11]+'.epiq' 


    ### file names: should be changed and passed by arguments 
    #annotateHeader='snp'
    #annotateTail='Blocks.variant_function'
    #mapHeader='mapped_'
    #mapTail='.map'
    #snpInfoHeader='snp'
    #snpInfoTail='BlocksList.avinput'
    #pedHeader='mapped_'
    #pedTail='.ped'
    #identifiers=['chr21','chr22']
    #filePath='/Users/weihaoge/Work/snpAnnotation/testCodes/' 
    #fonm='test'+'.epiq' 

    ### prepare spark session
    sc = SparkContext(appName="FilerSNP")
    spark = SparkSession(sc)

    ### read selected genes
    #print('file:///'+filePath+'/selectedGenes.txt')  
    geneFile = sc.textFile('file:///'+filePath+'/selectedGenes.txt') 
    geneSel=geneFile.map(lambda k: k.split('\t')[1]).collect()[1:] 

    ### combine all map, annotation, and ped files
    fileCollection = fileNameCreator(identifiers) 
    print("==========")  
    print(fileCollection)  
    print("==========")  


    ### filter with selected genes 
    pdDF=resultAssemble(filePath,fileCollection,geneSel) 

    ### write the genotype ".epiq" file 
    pdDF.to_csv(filePath+'/'+fonm,sep='\t',index_label='headerLine') 

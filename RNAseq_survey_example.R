############################################
library(readxl)
dir="" #Path where the data files are located
############################################
### STAR output
fs=list.files(sprintf("%s/bam",dir) )
fs1=sort(grep("Log.final.out$",fs,value=T))
rf=lapply(1:length(fs1),function(fi) apply( read.csv(file=sprintf("%s/bam/%s",dir,fs1[fi])
                                                     ,sep="\t",header=FALSE, stringsAsFactors=FALSE),2, trimws))
rfm1=cbind(Info=rf[[1]][,1], do.call("cbind", lapply(rf,function(m) m[,2])))
rfm1=rfm1[-c(1:3),]
colnames(rfm1)[-1]=sub("\\.(.*)","",fs1)
sm=apply(gsub("%","",rfm1[,-1]),1,as.numeric)
colnames(sm) = rfm1[,1]
cidx=(apply(sm,2,function(v) length(unique(v))) != 1)
sm=sm[,cidx]
cn=c("Number of input reads |" , "Uniquely mapped reads number |" ,  "Number of reads mapped to multiple loci |" ,
     "Number of reads mapped to too many loci |" , "Number of reads unmapped: too short |" ,"Number of reads unmapped: other |"  )  
sm=sm[,match(cn,colnames(sm))]
df=data.frame(total=sm[,1], unique=sm[,2],multi=sm[,3]+sm[,4], unmapped=sm[,5]+sm[,6])
rownames(df)=colnames(rfm1)[-1]

###################################################################
library(GEOquery)
library(AnnotationDbi)
library(ensembldb)
library(EnsDb.Hsapiens.v79)
edb=EnsDb.Hsapiens.v79
gid = keys(edb, keytype = "GENEID")
gdb = ensembldb::select(edb, keys=gid,  columns=c("ENTREZID", "SYMBOL", "GENEBIOTYPE"), keytype="GENEID")
colnames(gdb)= c("Ensembl", "Entrez", "Symbol", "GENEBIOTYPE")
e2s=unique(gdb[,2:3])

############################################
### RSEM -> TPM , expected count
############################################
countToTpm = function(counts, effLen){
  rate = log(counts) - log(effLen)
  er=exp(rate)
  idx= counts ==0 | effLen==0
  denom = log(sum(er[!idx],na.rm=T))
  tpm=exp(rate - denom + log(1e6))
  tpm[idx]=0
  return(tpm)
}
############################################
fs=list.files(sprintf("%s/expression",dir),pattern="genes.results")
m1=read.csv(file=sprintf("%s/expression/%s",dir,fs[1]), sep="\t", header=T, stringsAsFactors = F)	#48795

ginfo=gdb[match(m1$gene_id,gdb$Ensembl),]
ginfo[,1] = m1$gene_id
ginfo[is.na(ginfo)]="-"
pgid = ginfo[ginfo$GENEBIOTYPE == "protein_coding",1]
tpm = rcm = ptpm = matrix(NA, nr=nrow(m1), nc=length(fs))
for(fi in 1:length(fs)){
  mm = read.csv(file=sprintf("%s/expression/%s",dir,fs[fi]), sep="\t", header=T, stringsAsFactors = F)	#48795
  mm = mm[match(m1$gene_id,mm$gene_id),]
  tpm[,fi] = mm$TPM
  rcm[,fi] = mm$expected_count
  pidx=(mm$gene_id %in% pgid)
  pTPM=countToTpm(counts = mm$expected_count[pidx], effLen = mm$effective_length[pidx])
  mm$pTPM=0
  mm$pTPM[pidx] = pTPM
  ptpm[,fi] = mm$pTPM
}
colnames(rcm)= colnames(tpm) = colnames(ptpm) = gsub(".STAR.RSEM.genes.results","",fs)
rownames(rcm)= rownames(tpm) = rownames(ptpm) = m1$gene_id

############################################
sinfo=readxl::read_excel(path=sprintf("%s/sample_info_example.xlsx",dir), sheet =1)
sinfo= as.data.frame(sinfo)
sinfo$ID=gsub("#","",sinfo$ID)
sl=strsplit(sinfo$ID,split="-")
sinfo$Replicate= sapply(sl,function(v) v[2])

fs=list.files(sprintf("%s/bam",dir) )
fs1=sort(grep("Log.final.out$",fs,value=T))
rf=lapply(1:length(fs1),function(fi) apply( read.csv(file=sprintf("%s/bam/%s",dir,fs1[fi])
                                                     ,sep="\t",header=FALSE, stringsAsFactors=FALSE),2, trimws))
rfm1=cbind(Info=rf[[1]][,1], do.call("cbind", lapply(rf,function(m) m[,2])))
rfm1=rfm1[-c(1:3),]
colnames(rfm1)[-1]=sub("\\.(.*)","",fs1)
sm=apply(gsub("%","",rfm1[,-1]),1,as.numeric)
colnames(sm) = rfm1[,1]
cidx=(apply(sm,2,function(v) length(unique(v))) != 1)
sm=sm[,cidx]
cn=c("Number of input reads |" , "Uniquely mapped reads number |" ,  "Number of reads mapped to multiple loci |" ,
     "Number of reads mapped to too many loci |" , "Number of reads unmapped: too short |" ,"Number of reads unmapped: other |"  )  
sm=sm[,match(cn,colnames(sm))]
df=data.frame(total=sm[,1], unique=sm[,2],multi=sm[,3]+sm[,4], unmapped=sm[,5]+sm[,6])
rownames(df)=colnames(rfm1)[-1]

annot= cbind(sinfo,df[match(sinfo$ID,rownames(df)),])

save(rcm, tpm, ptpm, pgid, ginfo, annot, file=sprintf("%s/RNAseq_RSEM_results.Rdata",dir) )

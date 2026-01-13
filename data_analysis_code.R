#######################################################################
### 1. settings
#######################################################################
library(DESeq2)
library(edgeR)
dir = "" #Path where the data files are located
cells=c("A549","HEPG2","HT29")
#######################################################################
# Get functional gene sets from MSigDB (msigdbr 7.5.1, downloaded on 2023-07-19)
#######################################################################
BiocManager::install("msigdbr")
library(msigdbr)
gsets = msigdbr(species = "Homo sapiens")
gdf= as.data.frame(gsets)
save(gdf, file=sprintf("%s/msigdbr_hs.Rdata",dir), version=2 ) 

load( file=sprintf("%s/msigdbr_hs.Rdata",dir ))	#gl
gl= split(x = gdf$gene_symbol, f = gdf$gs_name)
gannot= unique(gdf[,c(1,2,3)])
gannot=gannot[match(names(gl), gannot$gs_name),]
n= sapply(gl,length)
gsidx= gannot$gs_cat=="H"| gannot$gs_subcat %in% c('CP:BIOCARTA','CP:KEGG','CP:REACTOME','CP','CP:PID','CP:WIKIPATHWAYS')
gsidx=gsidx & n>= 5 & n<= 500
gl=gl[gsidx]
gannot=gannot[gsidx,]


#######################################################################
### KIOM - Wortmannin DEG analysis
#######################################################################
load( file=sprintf("%s/RNAseq_results_example.Rdata",dir) ) #rcm, sin, ginfo
sin$Dose[is.na(sin$Dose)]=""
sin$Group=apply(sin[,c(5,6,3,12)],1,paste,collapse="-")
sin$Label=apply(sin[,c(5,6,3,12,4)],1,paste,collapse="-")
grp=unique(sin[,c(5,6,3,12,13)])
ll=lapply(cells, function(cell) data.frame(control=sprintf("Vehicle---%s",cell),case=grp$Group[grp$Cell==cell]))
comp=do.call(rbind,ll)
comp=comp[comp[,1]!=comp[,2],]

degl=list()
for( i in 1:nrow(comp)){
	print(comp[i,])
	con = comp$control[i]
	cas = comp$case[i]
	an=rbind(data.frame(ID=sin[sin$Group == con,"ID2"],fac=1)
		,data.frame(ID=sin[sin$Group == cas,"ID2"],fac=2))
	an=cbind(an,Batch=sin$Seq_Batch[match(an$ID,sin$ID2)])
	stpm= tpm[,match(an$ID,colnames(rcm))]
	crm=cor(stpm[,an$fac==2])
	v=crm[upper.tri(crm)]
	print(v)

	xm= round(rcm[,match(an$ID,colnames(rcm))],digits=0)
	group=as.factor(an$fac)
	batch=as.factor(an$Batch)

	design = model.matrix(~batch + group)
	keep = filterByExpr(xm, design, group=group )
	coldata=data.frame(Condition=group)
	rownames(coldata)= colnames(xm)
	dds=DESeqDataSetFromMatrix(countData = xm,colData = coldata, design = ~Condition)
	featureData = data.frame(gene=rownames(xm))
	mcols(dds) = DataFrame(mcols(dds), featureData)
	ddsMF= DESeq(dds,minReplicatesForReplace=Inf)
	resMF =results(ddsMF, cooksCutoff=F, pAdjustMethod = "fdr",independentFiltering=F)
	res=as.data.frame(resMF)
	res=data.frame(res[match(rownames(xm), rownames(res)),],keep, check.names=F)
	degl[[i]]=res

}
names(degl)=comp$case

save(degl, ginfo, sin, comp,
     file=sprintf("%s/RNAseq_results_example_DESeq_degl.Rdata", dir),
     version=2)
#######################################################################
### Get CMap data from the public repository (https://clue.io/data/CMap2020#)
#######################################################################
#siginfo.txt	Metadata for level 5 signatures
#instinfo.txt	Metadata for levels 3 and 4 profiles
#geneinfo.txt	Metadata for genes (the features of the data matrices)
#cellinfo.txt	Metadata for cell lines
#compoundinfo.txt	Metadata for cell lines (note there is one entry per compound/moa/target/structure combination, so some compounds appear more than once)
#level5_beta_trt_cp_n720216x12328.gctx
#20230719_Repurposing_Hub_export.txt #Informtation on approved drugs in CMap
#####################################################################
# Load CMap profiles (Wortmannin + approved drugs) matched to A549, HEPG2, HT29.
#####################################################################
BiocManager::install("cmapR")
library(cmapR)

cell_info=read.csv(file=sprintf("%s/cellinfo_beta.txt",dir),	header=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors =FALSE)	#240 x 20
comp_info=read.csv(file=sprintf("%s/compoundinfo_beta.txt",dir),	header=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors =FALSE)	#39321 x 7
gene_info=read.csv(file=sprintf("%s/geneinfo_beta.txt",dir),	header=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors =FALSE)	#12328 x 7
inst_info=read.csv(file=sprintf("%s/instinfo_beta.txt",dir),	header=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors =FALSE)	#3027596 x 30
sig_info=read.csv(file=sprintf("%s/siginfo_beta.txt",dir),	header=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors =FALSE)	#1202656  x 35

sig_info$median_recall_rank_wtcs_50[sig_info$median_recall_rank_wtcs_50 < 0]=NA
sig_info$median_recall_rank_spearman[sig_info$median_recall_rank_spearman < 0]=NA
sig_info$pct_self_rank_q25[sig_info$pct_self_rank_q25 <0]=NA
idx1 = sig_info$qc_pass >=1 & (sig_info$median_recall_rank_wtcs_50 <= 5 | sig_info$median_recall_rank_spearman <= 5)
idx2 = sig_info$cc_q75 >= 0.2 & sig_info$pct_self_rank_q25 <= 5 & sig_info$nsample >= 3
gold_idx= idx1 & idx2	#190253
sig_info=sig_info[gold_idx,]

hub=read.csv(file=sprintf("%s/20230719_Repurposing_Hub_export.txt",dir),sep="\t",header=T, check.names=F, stringsAsFactor=F)
idl=strsplit(hub$Id,split=", ")
n=sapply(idl,length)
hub=cbind(hub[rep(1:length(idl),n),!(colnames(hub) %in% c("Id"))],Id=unlist(idl))
hub$Id=substr(hub$Id,1,13)
hub=unique(hub)	#6783 x 7
df=cbind(comp_info,hub[match(comp_info$pert_id, hub$Id),])
idx1= is.na(df$Id) & df$cmap_name %in% hub$Name
df[idx1,9:15]=hub[match(df$cmap_name[idx1],hub$Name),]

fname=sprintf("%s/level5_beta_trt_cp_n720216x12328.gctx",dir)
rid=read_gctx_ids(fname, dim = "row")
cid=read_gctx_ids(fname, dim = "col")
idx=which(sig_info$pert_type == "trt_cp" & sig_info$pert_id %in% df$pert_id & sig_info$cell_iname %in%cells & sig_info$pert_time==24)
icid= intersect(cid, sig_info$sig_id[idx])
cmap_sinfo=sig_info[which(sig_info$sig_id %in% icid),]
sigs=cmap_sinfo$sig_id
print(length(sigs))
cidx = which(cid %in% sigs)
ridx = 1:length(rid)
ds = parse_gctx(fname, rid=ridx, cid=cidx)
m=ds@mat
identical(ds@rid, rownames(m))
identical(ds@cid, colnames(m))
cmap_sinfo= cmap_sinfo[match(ds@cid,cmap_sinfo$sig_id),]
cmap_ginfo=gene_info[match(ds@rid, gene_info$gene_id),]
cinfo=cell_info[cell_info$cell_iname %in% unique(cmap_sinfo$cell_iname),]
cmap_sinfo=cmap_sinfo[,c("sig_id","cmap_name","pert_time","pert_time_unit" ,"pert_dose","pert_dose_unit","cc_q75")]
cmapm=m[,match(cmap_sinfo$sig_id,colnames(m))]
rownames(cmapm)=cmap_ginfo$gene_symbol

#####################################################################
### CMap - Pathway enrichment score calculation
#####################################################################
library(fgsea)
resl=lapply(1:ncol(cmapm), function(i){
	print(i)
	rk= cmapm[,i]
	rk= sort(rk,decreasing=T)
	res = fgsea(pathways = gl, stats = rk, nperm=100000, nproc=1)
	res = as.data.frame(res)
	pes= -log10(res$padj)*sign(res$NES)
	pes[match(names(gl),res$pathway)]
	})
cmap_pesm=do.call("cbind",resl)
cmap_pesm[is.na(cmap_pesm)]=0
colnames(cmap_pesm)= colnames(cmapm)
rownames(cmap_pesm)= names(gl)

#####################################################################
### KIOM - Pathway enrichment score calculation
#####################################################################
load(file=sprintf("%s/RNAseq_results_example_DESeq_degl.Rdata",dir) ) #degl,ginfo,sin,comp

resl=lapply(1:length(degl),function(i){
	print(i)
	m=degl[[i]]
	m=m[m$keep,]
	rk=m$stat
	names(rk)=rownames(m)
	res = fgsea(pathways = gl, stats = rk, minSize=5, maxSize=500)
	res= as.data.frame(res)
	pes= -log10(res$padj)*sign(res$NES)
	pes[match(names(gl),res$pathway)]
})

kmap_pesm=do.call("cbind",resl)
kmap_pesm[is.na(kmap_pesm)]=0
colnames(kmap_pesm)= names(degl)
rownames(kmap_pesm)= names(gl)

#####################################################################
### KIOM vs. CMap comparison
#####################################################################

cl=colorRampPalette(c("gold","white","red"))(20)[c(3,13)]
conc=c(2.22,10,10)

for( i in 1:3){
  cell=cells[i]
  id1=paste("Wortmannin",cell,sep="---")
  kmap_wort=kmap_pesm[,colnames(kmap_pesm)==id1]

  ss=cmap_sinfo[which(cmap_sinfo$cmap_name=="wortmannin" & cmap_sinfo$pert_dose==conc[i]),]
  ds=cmap_sinfo[which(cmap_sinfo$cmap_name!="wortmannin"),]

  idx = colnames(cmap_pesm) %in% ss$sig_id
  idx2= colnames(cmap_pesm) %in% ds$sig_id
  cs=cor(cmap_pesm[,idx],kmap_wort)
  cd=cor(cmap_pesm[,idx2],kmap_wort)

 	# Example export (optional):
  #dd=data.frame(Pathway=names(kmap_wort), KORE_Map_Wortmannin = kmap_wort, CMap_Wortmannin=cmap_pesm[,idx])
  #write.csv(dd, file =sprintf("%s/Datafile_%s_PES.csv",dir,cell),row.names=F)

  cordf=rbind(data.frame(Condition="Same", Correlation=unlist(cs)),
     data.frame(Condition="Different", Correlation=unlist(cd)))

  mu = cordf %>% 
    group_by(Condition) %>%
    summarise(grp.mean = mean(Correlation))
    theme_set(
    theme_classic() +
      theme(legend.position = "top")
    )
  ggplot(cordf, aes(x = Correlation)) + 
  geom_density(aes(fill = Condition), alpha = 0.4) + 
    geom_vline(aes(xintercept = grp.mean, color = Condition), data = mu, linetype = "dashed") +
    scale_color_manual(values = c("#868686FF", "#EFC000FF")) +
    scale_fill_manual(values = c("#868686FF", "#EFC000FF"))
  ggsave(filename =sprintf("%s/Figure_%s_Comparison.tif",dir,cell), width = 7, height = 7,units = 'cm', device='tiff', dpi=300)


}





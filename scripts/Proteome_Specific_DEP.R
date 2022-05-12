rm(list = ls())

library(cowplot)
get_limma_sig = function(x,pv,fc)
{
  df = subset(x,adj.P.Val < pv & abs(logFC) > log2(fc))
  #df = subset(x,p.value < pv & abs(logFC) > log2(fc))
  if (nrow(df) > 0)
  {
    df$SIG.ID = rep(NA,nrow(df))
    df$SIG.ID[df$logFC > 0] = "UP"
    df$SIG.ID[df$logFC < 0] = "DN"
  }
  return(df)
}
source("scripts/R_functions/enrichment_functions.R")
##### load DEP data
depfiles = list.files(path = "Data/DEP",pattern = "DEP.paired.corrected.xls$",full.names = TRUE)
names(depfiles) = gsub("^(.*)/|\\.(.*)$","",depfiles)
print(depfiles)

# read data
dep.lst = lapply(depfiles,function(x) read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE))
for (i in 1:length(dep.lst))
{
  dep.lst[[i]]$cancer = rep(names(dep.lst)[i],nrow(dep.lst[[i]]))
  dep.lst[[i]]$gene.symbol = rownames(dep.lst[[i]])
}

dep.dat = do.call('rbind.data.frame',dep.lst)

##### create volcanoe plots
vobj.list = vector("list",length(dep.lst))
top.n= 5
pval.cutoff = 0.05
fc.cutoff = 2

for (i in 1:length(dep.lst))
{
  dat = dep.lst[[i]]
  vec = rep(NA,nrow(dat))
  vec[order(dat$logFC,decreasing = TRUE)[1:top.n]] = "UP"
  vec[order(dat$logFC,decreasing = FALSE)[1:top.n]] = "DN"
  dat$is.label = vec
  
  vobj = ggplot() + geom_point(data = dat,aes(x = logFC,y = -log10(adj.P.Val))) +
    geom_label_repel(data = subset(dat,!is.na(is.label)),aes(x = logFC,y = -log10(adj.P.Val),label = gene.symbol,colour = is.label),size = 6) + 
    labs(title = names(dep.lst)[i],x = "log2(FC)",y = "-log10(FDR)") + 
    geom_vline(xintercept = c(1,-1)*log2(fc.cutoff),colour = "red") + geom_hline(yintercept = -log10(pval.cutoff),colour = "red") + 
    theme_bw() + scale_colour_manual(values = c("UP" = "red","DN" = "blue")) + 
    guides(colour = "none") + 
    theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),plot.title = element_text(size = 21))
  vobj.list[[i]] = vobj
  rm(dat,vec,vobj)
}

png(file = "DEP.volcano.png",res = 300,width = 3200,height = 5500)
plot_grid(plotlist = vobj.list,ncol = 2)
dev.off()

###### MSigDB enrichments for DEP
if (TRUE)
{
  # create DEP signature list
  dep.sig = lapply(dep.lst,function(x) get_limma_sig(x,pv = pval.cutoff,fc = fc.cutoff))
  lapply(dep.sig,dim)
  
  dep.sig = unlist(lapply(dep.sig,function(x) split(x$gene.symbol,factor(x$SIG.ID))),recursive = FALSE)
  
  # get hallmark set
  hall = lapply(read.geneSet("Data/MSigDB/h.all.v6.2.symbols.gmt"),function(x) x[-1])
  bg = Reduce("union",lapply(dep.lst,function(x) x$gene.symbol))
  # run enrichments
  fet.df = perform.AllPairs.FET(geneSets1 = dep.sig,geneSets2 = hall,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
  fet.df$corrected.FET.pvalue = p.adjust(fet.df$FET_pvalue,"BH")
  fet.df$DEG.ID = gsub("^(.*)\\.","",fet.df$set1_Name)
  fet.df$cancer = gsub("\\.(.*)$","",fet.df$set1_Name)
  fet.df$set2_Name = gsub("_"," ",gsub("^HALLMARK_","",fet.df$set2_Name))
  
  write.table(fet.df,file = "DEP_MSigDB.FET.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  # get score per signature direction
  library(reshape2)
  sigmat = list(UP = acast(data = subset(fet.df,DEG.ID == "UP"),formula = set2_Name ~ cancer,value.var = "corrected.FET.pvalue",fun.aggregate = function(x) min(x,na.rm = T)),
                DN = acast(data = subset(fet.df,DEG.ID == "DN"),formula = set2_Name ~ cancer,value.var = "corrected.FET.pvalue",fun.aggregate = function(x) min(x,na.rm = T)))
  sig.count = lapply(sigmat,function(x) rowSums(x < 0.05,na.rm = TRUE))
  sig.score = lapply(sigmat,function(x) rowSums(-log10(x),na.rm = TRUE))
  
  # get order for hallmark terms and cancers
  deg.score = sig.score$UP - sig.score$DN
  deg.mat = -log10(sigmat$UP) + log10(sigmat$DN)
  hobj = hclust(dist(t(deg.mat)),"complete")
  cancer.tick = hobj$labels[hobj$order]
  term.tick = names(deg.score)[order(deg.score,decreasing = TRUE)]
  
  # finally, generate heatmap
  library(ggplot2)
  heat.obj = ggplot() + facet_grid(. ~ DEG.ID) + 
    geom_tile(data = fet.df,aes(x = cancer,y = set2_Name,fill = -log10(corrected.FET.pvalue)*(c("UP" = 1,"DN" = -1)[DEG.ID]))) + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + 
    geom_text(data = subset(fet.df,corrected.FET.pvalue < 0.05),aes(x = cancer,y = set2_Name,label = signif(enrichment.foldchange,2))) + 
    scale_x_discrete(limits = cancer.tick) + scale_y_discrete(limits = term.tick) + 
    guides(fill = guide_colorbar(title = "-log10(FET FDR)*sign(FC)")) + 
    labs(x = "Cancer",y = "MSigDB Hallmark") + 
    theme_bw() + 
    theme(axis.text.x = element_text(size = 17,angle = 45,vjust = 1,hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 19),
          legend.position = "bottom",legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),strip.text = element_text(size = 17))
  
  png(file = "DEP_Hallmark.FET.png",res = 300,width = 2700,height = 3500)
  print(heat.obj)
  dev.off()
  
  ## Also, make a side plot to complement number of cancers enriched for the signatures
  library(reshape)
  count.df = melt(do.call('cbind',sig.count));colnames(count.df) = c("signature","DEG.ID","counts")
  bar.obj = ggplot(data = count.df) + geom_bar(aes(x = counts,y = signature,fill = DEG.ID),stat = "identity") + 
    scale_fill_manual(values = c("UP" = "red","DN" = "blue")) + 
    labs(x = "#. Cancers enriched",y = "MSigDB Hallmark") + 
    scale_y_discrete(limits = term.tick) + 
    theme_bw() + 
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 17),
          axis.title.x = element_text(size = 19),
          axis.title.y = element_blank(),legend.position = "bottom",legend.direction = "horizontal",
          legend.title = element_text(size = 16), legend.text = element_text(size = 15),
          strip.text = element_text(size = 17))
  
  png(file = "DEP_Hallmark.signif_count.png",res = 300,width = 800,height = 3500)
  print(bar.obj)
  dev.off()
  
}

#################################### DEP and DEG comparison
####### load DEG data
deg.dat = read.delim(file = "Data/DEG/pancancer_DEG.age_gender_race_adjusted.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
colnames(deg.dat)[1:2] = c("cancer","gene.symbol")
deg.dat = subset(deg.dat,cancer %in% c("BRCA","LUAD","UCEC","KIRC","COAD","LIHC"))
deg.dat$cancer = gsub("LIHC","HCC",gsub("COAD","CRC",gsub("KIRC","CCRCC",deg.dat$cancer)))

####### Now, per cancer, add stats
cancer.id = names(depfiles)

pdata = data.frame()
for (i in 1:length(cancer.id))
{
  dep.df = subset(dep.dat,cancer == cancer.id[i])
  deg.df = subset(deg.dat,cancer == cancer.id[i])
  
  # match genes
  genes = intersect(deg.df$gene.symbol,dep.df$gene.symbol)
  dep.df = dep.df[match(genes,dep.df$gene.symbol),]
  deg.df = deg.df[match(genes,deg.df$gene.symbol),]
  
  # create data set
  dep.df = dep.df[,1:6];colnames(dep.df) = paste0("DEP.",colnames(dep.df))
  deg.df = deg.df[,3:8];colnames(deg.df) = paste0("DEG.",colnames(deg.df))
  
  data.df = data.frame(gene.symbol = genes,cancer.id = rep(cancer.id[i],length(genes)),dep.df,deg.df)
  
  # annotate genes
  vec = rep(NA,nrow(data.df))
  vec[data.df$DEP.adj.P.Val < 0.05 & data.df$DEP.logFC > 0] = "DEP"
  vec[data.df$DEG.adj.P.Val < 0.05 & data.df$DEG.logFC > 0] = "DEG"
  vec[data.df$DEP.adj.P.Val < 0.05 & data.df$DEP.logFC > 0 & data.df$DEG.adj.P.Val < 0.05 & data.df$DEG.logFC > 0] = "DEP & DEG"
  vec[data.df$DEP.adj.P.Val < 0.05 & data.df$DEP.logFC > 0 & data.df$DEG.adj.P.Val > 0.2] = "proteome-specific DEP" 
  data.df$gene.class = vec;
  
  # annotate DEG ID
  vec = rep(NA,nrow(data.df))
  
  
  # get top genes
  sig.df = subset(data.df,gene.class == "proteome-specific DEP")
  sig.df = sig.df[order(sig.df$DEP.logFC,decreasing = TRUE),]
  top.genes = sig.df$gene.symbol[1:5]
  
  data.df$is.top.gene = data.df$gene.symbol %in% top.genes
  
  # combine with other data
  pdata = rbind.data.frame(pdata,data.df)
  rm(dep.df,deg.df,genes,data.df)
}

# print out different classes of DEs
table(pdata[,c("gene.class","cancer.id")])


######### Now, check hallmark enrichments

bg = unique(pdata$gene.symbol)

# get hallmark set
hall = lapply(read.geneSet("Data/MSigDB/h.all.v6.2.symbols.gmt"),function(x) x[-1])
dep.specs = subset(pdata,gene.class == "proteome-specific DEP")
dep.specs = split(dep.specs$gene.symbol,factor(paste0(dep.specs$cancer.id)))

# run enrichments
fet.df = perform.AllPairs.FET(geneSets1 = dep.specs,geneSets2 = hall,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
fet.df$corrected.FET.pvalue = p.adjust(fet.df$FET_pvalue,"BH")
fet.df = fet.df[order(fet.df$FET_pvalue),]
fet.df$set2_Name = gsub("_"," ",gsub("^HALLMARK_","",fet.df$set2_Name))

write.table(fet.df,file = "Proteome_specific_DEP_MSigDB.FET.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

####### Now, create plot
library(ggplot2)
library(ggrepel)

pobj = ggplot() + 
  geom_point(data = pdata,aes(x = DEP.logFC,y = DEG.logFC,fill = sign(DEP.logFC)*-log10(DEP.adj.P.Val),colour = gene.class),shape = 21) + 
  geom_text_repel(data = subset(pdata,is.top.gene),aes(x = DEP.logFC,y = DEG.logFC,label = gene.symbol),max.overlaps = 10,size = 6) + 
  scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + 
  labs(x = "DEP FC",y = "DEG FC") + 
  scale_colour_manual(values = c("DEP" = "magenta","DEP & DEG" = "orange","proteome-specific DEP" = "red"),na.value = "grey",
                      limits = c("DEP","DEP & DEG","proteome-specific DEP")) + 
  geom_hline(yintercept = 0,colour = "red") + 
  #facet_grid(. ~ cancer.id,scale= "free_x") + 
  facet_wrap( ~ cancer.id,ncol = 3) + 
  guides(fill = guide_colorbar(title = "sign(FC)x-log10(FDR)",direction = "horizontal"),
         colour = guide_legend(title = "DEP type",direction = "vertical")) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 20),axis.text = element_text(size = 18),
        legend.title = element_text(size = 18),legend.text = element_text(size =16),strip.text = element_text(size = 17))

png(file = "proteome_specific_DEP.volcanoe.png",res = 400,width = 3500,height = 3800)
print(pobj)
dev.off()



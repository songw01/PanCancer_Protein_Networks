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

########################## DEP analysis
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
if (TRUE)
{
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
}

######## Perform super exact test on proteome-specific DEP
if (TRUE)
{
  library(SuperExactTest)
  
  df = subset(pdata,gene.class == "proteome-specific DEP")
  dep.spec.sigs = split(df$gene.symbol,df$cancer.id)
  bg = unique(pdata$gene.symbol)
  super.out = supertest(x = dep.spec.sigs,n = length(bg),degree = 2:4)
  
  png(file = "SuperExactTest_Barplot.ProteomeSpecific_DEP.png",res = 400,width = 4000,height = 2500)
  plot(super.out, Layout = "landscape",degree = c(2:4),sort.by="size")
  dev.off()
  
  write.table(summary(super.out)$Table, 
              file="Proteome_specific_DEP.SET.txt",sep = "\t", row.names=FALSE,col.names = TRUE,quote = FALSE)
  
  ####### Now, concatenate DEPs intersecting from at least 3 sets 
  tbl = summary(super.out)$Table
  tbl$adj.P.value = p.adjust(tbl$P.value,"BH")
  tbl = subset(tbl,Degree >= 3 & adj.P.value < 0.05 & FE > 2) 
  sigs = Reduce("union",strsplit(tbl$Elements,", "))
  
  tbl = summary(super.out)$Table
  tbl$adj.P.value = p.adjust(tbl$P.value,"BH")
  tbl = subset(tbl,Degree >= 4 & adj.P.value < 0.05 & FE > 2) 
  topsigs = Reduce("union",strsplit(tbl$Elements,", "))
  
  ##### now, plot heatmap for DEPs
  # prep data for heatmap
  dep.data = subset(pdata,gene.symbol %in% sigs)[,c(1,2,grep("^DEP",colnames(pdata)))];colnames(dep.data) = gsub("^DEP\\.","",colnames(dep.data));
  dep.data$data.type = rep("DEP",nrow(dep.data))
  
  deg.data = subset(pdata,gene.symbol %in% sigs)[,c(1,2,grep("^DEG",colnames(pdata)))];
  colnames(deg.data) = gsub("^DEG\\.","",colnames(deg.data));
  deg.data$data.type = rep("DEG",nrow(deg.data))
  
  full.data = rbind.data.frame(dep.data,deg.data)
  
  vec = rep(FALSE,nrow(full.data))
  for (i in 1:length(dep.spec.sigs))
  {
    ii = which(full.data$data.type == "DEP" & full.data$gene.symbol %in% dep.spec.sigs[[i]] & full.data$cancer.id == names(dep.spec.sigs)[i])
    vec[ii] = TRUE
  }
  full.data$is.specific.DEP = vec
  
  full.data$gene.symbol = factor(full.data$gene.symbol)
  
  # now, order genes by patterns in DEP fold change
  require(reshape2)
  fcmat = acast(data = dep.data,formula = cancer.id + data.type~ gene.symbol,value.var = "logFC",fun.aggregate = function(x) mean(x,na.rm = TRUE))
  fcmat[is.nan(fcmat)] = 0
  hobj = hclust(dist(t(fcmat)),"complete")
  gene.labs = hobj$labels[hobj$order]
  
  # now, create heatmap 
  genes = gene.labs
  lab.fonts = rep("plain",length(genes))
  lab.fonts[gene.labs %in% topsigs] = "bold"
  lab.col = rep("black",length(genes))
  lab.col[gene.labs %in% topsigs] = "magenta"
  library(ggplot2)
  pobj = ggplot() + geom_tile(data = full.data,aes(x = gene.symbol,y = cancer.id,fill = logFC)) + 
    geom_point(data = subset(full.data,is.specific.DEP),aes(x = gene.symbol,y = cancer.id),colour = "magenta") + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + facet_grid(data.type ~ .) + 
    labs(x = "Gene",y = "Cancer type") + 
    scale_x_discrete(limits = gene.labs) + 
    guides(fill = guide_colorbar(title = "log2(FC)")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 16,face = lab.fonts,colour = lab.col),axis.title = element_text(size = 20),
          axis.text.y = element_text(size = 17),
          legend.title = element_text(size = 18),legend.text = element_text(size = 15),
          strip.text = element_text(size = 20))
  png(file = "proteome_specific_DEP.logFC_heatmap.png",res = 400,width = 12300,height = 3000)
  print(pobj)
  dev.off()
}

####### Now, check hallmark enrichments
if (TRUE)
{
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
}

#################################### Cancer specific DEPs
if (TRUE)
{
  # Now, calculate cancer type-specific DEPs
  cancer.id = names(dep.lst)
  depdat = do.call('rbind.data.frame',dep.lst)
  
  vec = rep(NA,nrow(depdat))
  for (i in 1:length(dep.sig))
  {
    ii = which(depdat$cancer == gsub("\\.(.*)$","",names(dep.sig)[i]) & depdat$gene.symbol %in% dep.sig[[i]])
    val = gsub("^(.*)\\.","",names(dep.sig)[i])
    vec[ii] = val;
  }
  
  depdat$DEP.ID = vec;
  
  # give scores to each direction
  vec = depdat$logFC;
  vec[depdat$logFC < 0] = 0
  depdat$logFC.UP = vec
  
  vec = depdat$logFC;
  vec[depdat$logFC > 0] = 0
  depdat$logFC.DN = vec
  
  ## get up-regulated specific sigs
  library(reshape2)
  library(reshape)
  library(ggplot2)
  # make logFC matrix
  fcmat = acast(data = depdat,formula = gene.symbol ~ cancer,value.var = "logFC.UP") 
  pmat = acast(data = depdat,formula = gene.symbol ~ cancer,value.var = "adj.P.Val")
  mat = fcmat * (pmat < 0.2) 
  
  # look for specific matrix
  is.only = apply(mat,1,function(x) {out = NA;if (sum(x > 1E-320,na.rm = TRUE) == 1) out = names(x)[which(x > 1E-320)];return(out)})
  is.fc = apply(mat,1,function(x) {any(x > log2(fc.cutoff))})
  is.valid = apply(fcmat,1,function(x) sum(is.na(x)) <= 3)
  
  up.spec = split(names(is.only)[is.fc & is.valid],factor(is.only[is.fc & is.valid]))
  sapply(up.spec,length)
  
  # sort DEPs by p-values
  up.spec.dat = vector("list",length(up.spec))
  names(up.spec.dat) = names(up.spec)
  
  vec = rep(NA,nrow(depdat))
  for (i in 1:length(up.spec))
  {
    dat = subset(depdat,gene.symbol %in% up.spec[[i]] & cancer == names(up.spec)[i])
    dat = dat[order(dat$P.Value),]
    dat = subset(dat,adj.P.Val < 0.05)
    up.spec.dat[[i]] = dat;
    
    # update depdat with specificity
    ii = which(depdat$cancer == names(up.spec)[i] & depdat$gene.symbol %in% dat$gene.symbol)
    vec[ii] = paste0(names(up.spec)[i],"-specific DEP")
    rm(dat)
  }
  
  depdat$cancer.specificity = vec
  
  sapply(up.spec.dat,dim)
  write.table(depdat,file = "Cancer_specific_DEP.txt",sep = "\t",
              row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  
  #### Now, make heatmap plots per cancer
  top.n = 5;
  top.genes = do.call('rbind',lapply(up.spec.dat,function(x) cbind(x$gene.symbol[1:min(c(top.n,nrow(x)))],x$cancer[1:min(c(top.n,nrow(x)))])))
  common.genes = c("PLOD2","PYCR1","TOP2A","MCM7","MKI67","MCM5","FEN1","RSL1D1","SMC2","SOAT1","GNL2")
  
  plot.data = subset(depdat,gene.symbol %in% top.genes[,1] | gene.symbol %in% common.genes)
  plot.data$cancer.specificity = top.genes[match(plot.data$gene.symbol,top.genes[,1]),2]
  plot.data$cancer.specificity[plot.data$gene.symbol %in% common.genes] = "COMMON"
  plot.data$cancer.specificity = factor(plot.data$cancer.specificity,levels= c(cancer.id,"COMMON"))
  
  vec = rep(NA,nrow(plot.data))
  vec[plot.data$adj.P.Val < 0.05 & plot.data$logFC > log2(fc.cutoff)] = "UP"
  vec[plot.data$adj.P.Val < 0.05 & plot.data$logFC < -log2(fc.cutoff)] = "DN"
  plot.data$DEG.ID = vec;
  
  fcobj = ggplot() + geom_tile(data = plot.data,aes(x = cancer,y = gene.symbol,fill = logFC)) + facet_grid(cancer.specificity ~ .,scale = "free_y",space = "free_y") +
    geom_point(data = subset(plot.data,adj.P.Val < 0.05),aes(x = cancer,y = gene.symbol,colour = DEG.ID),shape = 16) + 
    scale_fill_gradient2(na.value = "white",low = "blue",mid = "white",high = "red",midpoint = 0) + theme_bw() + 
    scale_colour_manual(values = c("UP" = "red","DN" = "blue"),limits = c("UP","DN")) + 
    guides(fill = guide_colorbar(title = "log2(FC)"),colour = guide_legend(title = "DEP")) + 
    labs(x = "Cancer",y = "Cancer type DEPs") + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 15),axis.text.y = element_text(size= 15),
          axis.title = element_text(size = 17),legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),strip.text = element_text(size = 17))
  
  # plot out cancer specificity
  png(file = "Cancer_specific_DEP.png",res = 300,width = 1800,height = 3900)
  print(fcobj)
  dev.off()
}

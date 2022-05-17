rm(list = ls())

library(igraph)
library(RCy3)
library(reshape2)
library(gtable)
library(ggplot2)
library(grid)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#root.dir = "C:/Users/songw01/Documents/CPTAC_Proteome_v2";setwd(root.dir)
out.dir = "Network_Connectivity_Analysis";dir.create(out.dir)

common.can = c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")
################# Load networks
## load proteome data
if (TRUE)
{
  # load MEGENA networks
  if (TRUE)
  {
    require(igraph)
    net = read.delim(file = "Data/MEGENA/MEGENA_Network.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    net.lst = lapply(split(1:nrow(net),factor(net$group)),function(ii,x) {out = x[ii,1:3];out[[1]] = gsub("\\|(.*)$","",out[[1]]);out[[2]] = gsub("\\|(.*)$","",out[[2]]);subset(out,row != col)},x = net)
    
    networks <- lapply(net.lst,function(x) {graph.data.frame(x)})
    networks <- networks[sapply(networks,vcount) > 500]
    networks = networks[common.can]
  }
  
  # load hub results
  if (TRUE)
  {
    filenames = list.files(path = "Data/MEGENA/hubs/proteome",pattern = "^(.*)\\.global_hubs\\.txt",full.names = TRUE)
    names(filenames) = gsub("^(.*)/|\\.(.*)$","",filenames)
    hub.lst = lapply(filenames,function(x) subset(read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE),hub.p.value <= 1))
    hub.lst = hub.lst[common.can]
    hub.lst = lapply(hub.lst,function(node.tbl) {
      # add node stratification data 
      hub.type = rep("Non-hub",nrow(node.tbl))
      hub.type[node.tbl$hub.p.value < 0.1] = "Local hub"
      hub.type[node.tbl$hub.p.value < 0.05] = "Sub-global hub"
      hub.type[node.tbl$hub.p.value < 0.01] = "Global hub"
      node.tbl$hub.type = hub.type
      return(node.tbl)
    })
  }
}

## load transcriptome data
if (TRUE)
{
  # load hub results
  if (TRUE)
  {
    filenames = list.files(path = "Data/MEGENA/hubs/transcriptome",pattern = "^(.*)\\.global_hubs\\.txt",full.names = TRUE)
    names(filenames) = gsub("^(.*)/|\\.(.*)$","",filenames)
    hub.lst.tx = lapply(filenames,function(x) subset(read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE),hub.p.value <= 1))
    hub.lst.tx= hub.lst.tx[common.can]
    hub.lst.tx = lapply(hub.lst.tx,function(x) {out = x;out[[1]] = gsub("\\|(.*)$","",out[[1]]);return(out)})
    hub.lst.tx = lapply(hub.lst.tx,function(node.tbl) {
      # add node stratification data 
      hub.type = rep("Non-hub",nrow(node.tbl))
      hub.type[node.tbl$hub.p.value < 0.1] = "Local hub"
      hub.type[node.tbl$hub.p.value < 0.05] = "Sub-global hub"
      hub.type[node.tbl$hub.p.value < 0.01] = "Global hub"
      node.tbl$hub.type = hub.type
      return(node.tbl)
    })
  }
}

######### classify PR vs TX hubs
if (TRUE)
{
  
  library(ggplot2)
  library(ggrepel)
  all.dat = data.frame()
  for (i in 1:length(common.can))
  {
    tx = hub.lst.tx[[i]]
    pr = hub.lst[[i]]
    
    common.gene = intersect(tx[[1]],pr[[1]])
    
    tx = tx[match(common.gene,tx[[1]]),]
    pr = pr[match(common.gene,pr[[1]]),]
    
    vec = rep(NA,length(common.gene))
    vec[which(tx$hub.type %in% c("Global hub","Sub-global hub","Local hub") & pr$hub.type %in% c("Global hub","Sub-global hub","Local hub"))] = "shared"
    vec[which(tx$hub.type %in% c("Global hub","Sub-global hub","Local hub") & pr$hub.p.value > 0.2)] = "TX-specific"
    vec[which(pr$hub.type %in% c("Global hub","Sub-global hub","Local hub") & tx$hub.p.value > 0.2)] = "PR-specific"
    
    df.pr = pr[,-1];colnames(df.pr) = paste0("PR.",colnames(df.pr))
    df.tx = tx[,-1];colnames(df.tx) = paste0("TX.",colnames(df.tx))
    
    dat = data.frame(gene = common.gene,df.pr,df.tx,comparison.type = vec,cancer = rep(common.can[i],nrow(df.pr)))
    
    
    # mark top 5 genes per category
    is.top = rep(FALSE,nrow(dat))
    
    ii = which(dat$comparison.type == "shared")
    score = -log10(dat$TX.hub.p.value) + -log10(dat$PR.hub.p.value);
    score = score[ii];k = order(score,decreasing = T);ii = ii[k];ii = ii[1:5];is.top[ii] = TRUE
    
    ii = which(dat$comparison.type == "TX-specific")
    score = -log10(dat$TX.hub.p.value) 
    score = score[ii];k = order(score,decreasing = T);ii = ii[k];ii = ii[1:5];is.top[ii] = TRUE
    
    ii = which(dat$comparison.type == "PR-specific")
    score = -log10(dat$PR.hub.p.value) 
    score = score[ii];k = order(score,decreasing = T);ii = ii[k];ii = ii[1:5];is.top[ii] = TRUE
    
    dat$is.top = is.top
    
    all.dat = rbind.data.frame(all.dat,dat)
    
    rm(tx,pr,common.gene,vec,df.pr,df.tx)
  }
  
  pobj = ggplot() + geom_point(data = all.dat, aes(x = -log10(PR.hub.p.value),y = -log10(TX.hub.p.value),
                                                   fill = comparison.type,alpha = comparison.type,
                                                   size = -log10(PR.hub.p.value) + -log10(TX.hub.p.value)),shape = 21,stroke = 0) + 
    labs(x = "-log10(Hub p-value: Proteome)",y = "-log10(Hub p-value: Transcriptome)") + 
    geom_text_repel(data = subset(all.dat,is.top),aes(x = -log10(PR.hub.p.value),y = -log10(TX.hub.p.value),label = gene),size = 6,fontface = "bold") + 
    guides(colour = guide_legend(title = "Hub type by\nsignificance"),
           fill = guide_legend(title = "TX & PR\n Comparison",override.aes = list(size = 4)),size = "none",alpha = "none") + 
    geom_hline(data = data.frame(threshold = c(0.01,0.05,0.1),type = rev(c("Local","Sub-global","Global"))),mapping = aes(yintercept = -log10(threshold),colour = type)) +
    geom_vline(data = data.frame(threshold = c(0.01,0.05,0.1),type = rev(c("Local","Sub-global","Global"))),mapping = aes(xintercept = -log10(threshold),colour = type)) +
    scale_colour_manual(values = c("Local" = "yellow","Sub-global" = "pink","Global" = "red")) + 
    scale_alpha_manual(values = c("shared" = 1,"TX-specific" = 1,"PR-specific" = 1),na.value = 0.01) + 
    scale_fill_manual(values = c("shared" = "magenta","TX-specific" = "green","PR-specific" = "red")) + 
    facet_wrap( ~ cancer) + theme_bw() + 
    theme(axis.title = element_text(size = 24),axis.text = element_text(size = 20),legend.title = element_text(size = 17),legend.text = element_text(size = 14),
          strip.text = element_text(size = 19))
  #pobj
  
  png(file = paste0(out.dir,"/TX_vs_PR.hub_pvalue_comparison.png"),res = 500,width = 8000,height = 8000)
  print(pobj)
  dev.off()
  
  ######## Perform pathway enrichments on different hubs
  if (TRUE)
  {
    source("scripts/R_functions/enrichment_functions.R")
    bg = unique(all.dat$gene)
    
    # get hallmark set
    hall = lapply(read.geneSet("Data/MSigDB/h.all.v6.2.symbols.gmt"),function(x) x[-1])
    
    # get connectivity signatures
    df= subset(all.dat,!is.na(comparison.type))
    con.sig=split(as.character(df$gene),factor(paste0(df$cancer,"_",df$comparison.type)))
    
    # run enrichments
    fet.df = perform.AllPairs.FET(geneSets1 = con.sig,geneSets2 = hall,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    fet.df$corrected.FET.pvalue = p.adjust(fet.df$FET_pvalue,"BH")
    fet.df = fet.df[order(fet.df$FET_pvalue),]
    fet.df$set2_Name = gsub("_"," ",gsub("^HALLMARK_","",fet.df$set2_Name))
    fet.df$signature.type = gsub("^(.*)_","",fet.df$set1_Name)
    fet.df$cancer = gsub("_(.*)$","",fet.df$set1_Name)
    
    # order hallmark signatures by enrichment patterns
    pdata = subset(fet.df,set2_Name %in% subset(fet.df,corrected.FET.pvalue < 0.001)$set2_Name)
    mat = acast(data = pdata,formula = set2_Name ~ cancer + signature.type,value.var = "corrected.FET.pvalue")
    hobj = hclust(as.dist(sqrt(2*(1-cor(t(-log10(mat)))))))
    htick = hobj$label[hobj$order]
    
    pobj = ggplot() + geom_tile(data = pdata,aes(x = signature.type,y = set2_Name,fill = -log10(corrected.FET.pvalue))) + 
      geom_point(data = subset(pdata,corrected.FET.pvalue < 0.05),mapping = aes(x = signature.type,y = set2_Name,size = enrichment.foldchange),colour = "black") + 
      scale_y_discrete(limits = htick) + 
      scale_x_discrete(limits = c("PR-specific","TX-specific","shared")) + 
      facet_grid(. ~ cancer) + 
      scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = -log10(0.05)) + 
      theme_bw() + 
      guides(fill = guide_colorbar(title = "-log10(FET FDR)"),size = guide_legend(title = "EFC")) + 
      labs(x = "Connectivity Signature",y = "Hallmark") + 
      theme(axis.text.x = element_text(angle = 45,vjust =1 ,hjust = 1,size = 17),axis.text.y = element_text(size = 17),
            axis.title = element_text(size = 20),strip.text =element_text(size = 18),
            legend.title = element_text(size = 18),legend.text = element_text(size = 16),
            legend.position = "bottom",legend.direction = "horizontal")
    png(file = paste0(out.dir,"/TX_vs_PR.hub_signature_hallmark.png"),res = 500,width = 6000,height = 4500)
    print(pobj)
    dev.off()
    
  }

}

######### summarize hub statistics
common.gene = Reduce("intersect",lapply(hub.lst,function(x) x$gene))

hub.pmat = do.call('cbind',lapply(hub.lst,function(x) x$hub.p.value[match(common.gene,x[[1]])]))
rownames(hub.pmat) = common.gene

# subset for varying connectivity
hub.pmat = hub.pmat[which(rowSums(hub.pmat < 0.05,na.rm = TRUE) > 0),]

# extract genes to label for multiplicity
top.genes = rownames(hub.pmat)[rowSums(hub.pmat < 0.01,na.rm = TRUE) >= 3]

# transform into -log10(p) matrix
hub.neglogmat = -log10(hub.pmat)
hobj = hclust(dist(hub.neglogmat),"complete")

######### summarize pan-cancer mutation hub driver connectome: DDX21, RSL1D1, SMC2
if (TRUE)
{
  
  source("scripts/R_functions/enrichment_functions.R")
  
  ######## Load DEP results
  depfiles = list.files(path = "Data/DEP",pattern = "DEP.paired.corrected.xls$",full.names = TRUE)
  names(depfiles) = gsub("^(.*)/|\\.(.*)$","",depfiles)
  print(depfiles)
  
  # extract signatures
  depres = lapply(depfiles,function(x) {
    df = read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE);
    colnames(df)[4:5] = c("p.value","adj.p.value")
    df = data.frame(gene.symbol = rownames(df),df,stringsAsFactors = FALSE)
    df = df[order(df$p.value),]
    return(df)
  })
  
  bg = Reduce("union",lapply(depres,function(x) x$gene.symbol))
  
  # set up MSigDB signatures
  gmt.folder <- "Data/MSigDB";
  gmt.files <- list.files(path = gmt.folder,pattern = "\\.gmt$",full.names = T)
  names(gmt.files) = gsub("^(.*)/|\\.symbols\\.gmt","",gmt.files)
  
  # get genes
  genes = c("DDX21","RSL1D1","SMC2")
  
  fet.res = data.frame()
  for (gi in 1:length(genes))
  {
    neigh = lapply(net.lst,function(x) Reduce("union",subset(x,row == genes[gi] | col == genes[gi])[,c(1,2)]))
    for (fi in 1:length(gmt.files))
    {
      # get hallmark set
      hall = lapply(read.geneSet(gmt.files[fi]),function(x) x[-1])
      
      # run enrichments
      fet.df = perform.AllPairs.FET(geneSets1 = neigh,geneSets2 = hall,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
      fet.df$corrected.FET.pvalue = p.adjust(fet.df$FET_pvalue,"BH")
      fet.df$set2_Name = gsub("_"," ",gsub("^HALLMARK_","",fet.df$set2_Name))
      fet.df$database = rep(names(gmt.files)[fi],nrow(fet.df))
      fet.df$gene = rep(genes[gi],nrow(fet.df))
      fet.res = rbind.data.frame(fet.res,fet.df)
      rm(fet.df)
    }
    
  }
  
  write.table(fet.res,file = paste0(out.dir,"/pancancer_hubs_immediate_pathways.FET.txt"),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
}


######### load surfaceosome data
library(readxl)
surf.tbl = read_xlsx(path = "Data/surfaceosome/table_S3_surfaceome.xlsx",
                     sheet = 1,
                     range = NULL,
                     col_names = TRUE,
                     col_types = NULL,
                     na = "",
                     trim_ws = TRUE,
                     skip = 1
)

surf.tbl = surf.tbl[which(surf.tbl$"MachineLearning FPR class (1=1%, 2=5%, 3=15%)" %in% c(1) & surf.tbl$`Surfaceome Label` == "surface" & surf.tbl$`Membranome Almen main-class` %in% c("Enzymes","Receptors")),]

# add in LINCS validation data
lincs.dat = subset(data.res,padj < 0.05 & layer < 3)
lincs.lst = lapply(split(lincs.dat$target.gene,factor(lincs.dat$cancer)),unique)
lincs.lst.surface = lapply(lincs.lst,function(x) intersect(x,surf.tbl$`UniProt gene`))

######### Summarize and plot LINCS and CRisPRi hit results
if (TRUE)
{
  cancers = c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")
  
  ######## summarize LINCS evaluation
  res = rbind.data.frame(
    read.delim(file = "Data/KD_stratification/KD_Stratification_Evaluaton.Proteome.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE),
    read.delim(file = "Data/KD_stratification/KD_Stratification_Evaluaton.Transcriptome.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  )
  
  # summarize LINCS validation
  require(ggplot2)
  plot.data = subset(res,group.type == "HUB_CALL" & validation.type == "LINCS")
  plot.data$group = factor(plot.data$group,levels = c("P_0.01","P_0.05","P_0.1","non-hub"))
  rate.obj = ggplot(data = plot.data,aes(x= cancer,y = Rate.true,fill = group)) + geom_bar(stat = "identity",position = "dodge") + 
    scale_fill_manual(labels = c("P_0.01" = "Global hub (Hub p < 0.01)","P_0.05" = "Sub-global hub (0.01 < Hub p < 0.05)","P_0.1" = "Local hub (0.05 < Hub p < 0.1)",
                                 "non-hub" = "Non-hub (Hub p > 0.1)",
                                 "global.KD" = "Global KD","local.KD" = "Local KD"),
                      values = c("P_0.01" = "red","P_0.05" = "coral","P_0.1" = "darkgoldenrod1","non-hub" = "grey")) + 
    scale_y_continuous(labels=scales::percent) +
    facet_grid(. ~ data.type) + theme_bw() + labs(y = "Validation Rate(%)") + 
    guides(fill = guide_legend(title = "Node\nCategory",ncol = 2)) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 17),axis.title.x = element_blank(),axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 19),strip.text = element_text(size = 16),legend.position = "bottom",legend.direction = "horizontal",
          legend.text = element_text(size = 13),legend.title = element_text(size = 17)) + guides(fill = FALSE)
  
  
  # summarize LINCS validation
  require(ggplot2)
  plot.data = subset(res,group.type == "HUB_CALL" & validation.type == "CRISPRi")
  plot.data$group = factor(plot.data$group,levels = c("P_0.01","P_0.05","P_0.1","non-hub"))
  FET.obj = ggplot(data = plot.data,aes(x = cancer,y = Rate.true,fill = group)) + 
    geom_bar(stat = "identity",position = "dodge") + 
    scale_fill_manual(labels = c("P_0.01" = "Global hub (Hub p < 0.01)","P_0.05" = "Sub-global hub (0.01 < Hub p < 0.05)","P_0.1" = "Local hub (0.05 < Hub p < 0.1)",
                                 "non-hub" = "Non-hub (Hub p > 0.1)",
                                 "global.KD" = "Global KD","local.KD" = "Local KD"),
                      limits = c("P_0.01","P_0.05","P_0.1","non-hub"),
                      values = c("P_0.01" = "red","P_0.05" = "coral","P_0.1" = "darkgoldenrod1","non-hub" = "grey")) + 
    theme_bw() + 
    facet_grid(. ~ data.type) + 
    scale_y_continuous(labels=scales::percent) +
    labs(y = "Validation Rate (%)") + 
    guides(fill = guide_legend(title = "Node\nCategory",ncol = 2)) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 17),
          axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 19),strip.text = element_text(size = 15),legend.position = "bottom",legend.direction = "horizontal",
          legend.text = element_text(size = 16),legend.title = element_text(size = 17))
  
  ### extract legend for plotting separately
  
  
  mylegend <- g_legend(FET.obj)
  FET.obj = FET.obj + guides(fill = "none")
  
  library(cowplot)
  pdf(file = paste0(out.dir,"/Node_Stratification_Validation.CRISPRi_and_LINCS.pdf"),width = 9,height = 9)
  plot_grid(rate.obj,FET.obj,mylegend,ncol = 1,nrow = 3,rel_heights = c(0.92,1,0.25))
  dev.off()
  
}

######## Make genewise immediate neighborhood pathway plots
if (TRUE)
{
  library(ggraph)
  library(igraph)
  library(scatterpie)
  library(reshape2)
  library(ggrepel)
  library(cowplot)
  
  cancer.col = c("BRCA" = "red","CCRCC" = "orange","CRC" = "thistle4","HCC" = "darkgreen",
                 "LUAD" = "blue","MB" = "deepskyblue",
                 "STAD" = "slateblue","UCEC" = "darkturquoise","PANCAN" = "black");
  
  net.list = net.lst[c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")]
  
  genes = c("DDX21","RSL1D1","SMC2")
  plst = vector("list",length(genes));names(plst) = genes
  for (i in 1:length(genes))
  {
    gene = genes[i]
    net.el =kvec =  data.frame()
    for (l in 1:length(net.list)) 
    {
      # retrieve neighborhood topology
      v = Reduce("union",subset(net.list[[l]],row == gene | col == gene)[,c("row","col")])
      el = subset(net.list[[l]],row %in% v | col %in% v)
      el$cancer = rep(names(net.list)[l],nrow(el))
      net.el = rbind.data.frame(net.el,el)
      
      # get cancer-wise connectivity
      cg = graph.data.frame(el,directed = FALSE)
      k = igraph::degree(cg);
      kvec = rbind.data.frame(kvec,data.frame(gene = names(k),degree = k,cancer = rep(names(net.list)[l],length(k))))
      
      rm(v,el,k,cg)
    }
    
    
    # get neighborhood interactions
    g = graph.data.frame(net.el,directed = FALSE)
    V(g)$degree <- igraph::degree(g)
    
    # draw main network
    pobj = ggraph(graph = g) + 
      geom_edge_arc(aes(edge_colour = factor(cancer)),alpha = 0.1,strength = 0.3,show.legend = FALSE) + 
      scale_edge_color_manual(values = cancer.col) + 
      #geom_node_point(aes(size = degree)) + 
      theme_void() 
    
    # summarize cancer-wise connectivity, and add as pie chars
    kmat = acast(data = kvec,formula = gene ~ cancer,value.var = "degree",fun.aggregate = function(x) mean(x,na.rm = T))
    kmat[is.nan(kmat)] = 0
    kdata = data.frame(gene = rownames(kmat),as.data.frame(kmat),pobj$data[match(rownames(kmat),pobj$data$name),c("x","y")],total.degree = igraph::degree(g)[rownames(kmat)])
    vec = log(kdata$total.degree);vec = vec/max(vec,na.rm = T);vec[vec < 1E-320] = 0.5*min(vec[vec >1E-320])
    vec = vec * 0.1
    kdata$radius = vec;
    
    pobj = pobj + geom_scatterpie(aes(x=x, y=y,r = radius), data=kdata,
                                  cols=names(net.list),color = NA,alpha = 0.4) + scale_fill_manual(values = cancer.col[names(net.list)])  + 
      guides(colour = "none",fill = guide_legend(title = "Cancer type",ncol = 1)) + 
      labs(title = paste0(gene,"-centered Network")) + 
      theme(legend.position = "right",legend.direction = "vertical",plot.title = element_text(hjust = 0.5,size = 17))
    
    # add labels for top connected nodes
    txtdata = subset(kdata,total.degree >= quantile(kdata$total.degree,0.95));
    txtdata$is.center = txtdata$gene == gene
    pobj = pobj + geom_text_repel(data = txtdata,aes(x = x,y = y,label = gene,size = log(total.degree)*1.4,colour = is.center)) + 
      scale_colour_manual(values = c("TRUE" = "red","FALSE" = "black")) + 
      guides(size = "none",colour = "none") 
    plst[[i]] = pobj
  }
  gl = g_legend(plst[[1]])
  plst = lapply(plst,function(x) x + guides(fill = "none"))
  plst = c(plst,list(gl))
 
  pdf(file = paste0(out.dir,"/top_driver_networks.pdf"),width = 17,height = 6.1)
  plot_grid(plotlist = plst,ncol = 4,rel_widths = c(1,1,1,0.2))
  dev.off()
}

############################### Hub heatmap analysis
if (TRUE)
{
  ######### identify cancer specific hubs
  cancer.col = c("BRCA" = "red","CCRCC" = "orange","CRC" = "thistle4","HCC" = "darkgreen",
                 "LUAD" = "blue","MB" = "deepskyblue",
                 "STAD" = "slateblue","UCEC" = "darkturquoise","PANCAN" = "black");
  m = hub.pmat[which(rowSums(hub.pmat < 0.05) == 1),]
  
  lst = apply(m,2,function(x) names(x)[x < 0.05])
  vec = rep(NA,nrow(hub.neglogmat))
  for (i in 1:length(lst)) vec[rownames(hub.neglogmat) %in% lst[[i]]] = names(lst)[i]
  vec[rownames(hub.neglogmat) %in% top.genes] = "PANCAN"
  
  ######### load PANCAN mutation data
  library(reshape2)
  mut.data = read.delim(file = "Data/MutationalDrivers/mutation_driver_list.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  mut.data$Cancer.type = gsub("LIHC","HCC",gsub("COADREAD","CRC",gsub("KIRC","CCRCC",mut.data$Cancer.type)))
  mut.data = subset(mut.data,Cancer.type %in% c(colnames(hub.pmat),"PANCAN"))
  
  mut.matrix = acast(data = mut.data,formula = Gene ~ Cancer.type,value.var = "Pathway")
  mut.matrix = mut.matrix[match(rownames(hub.neglogmat),rownames(mut.matrix)),]
  mut.matrix = mut.matrix[,c(colnames(hub.neglogmat),"PANCAN")]
  
  ######### Perform k-means to split into groups
  set.seed(1234)
  kmout = kmeans(x = hub.neglogmat,iter.max = 20,centers = 9)
  kmres = kmout$cluster[rownames(hub.neglogmat)]
  
  ######### Plot out heatmap
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  
  # generate markers for driver mutation intersection
  clus.map = c("6" = "CCRCC","5" = "HCC","8" = "LUAD","1" = "BRCA","9" = "UCEC","4" = "STAD","2" = "CRC","7" = "PANCAN","3" = "PANCAN")
  at = labs = c()
  for (i in 1:length(clus.map))
  {
    ii = which(as.character(kmres) == names(clus.map)[i] & rownames(hub.neglogmat) %in% subset(mut.data,Cancer.type == clus.map[i])[[1]])
    at = c(at,ii)
    labs = c(labs,rownames(hub.neglogmat)[ii])
  }
  
  # re-label gene clusters
  clus.map = c("6" = "CCRCC","5" = "HCC","8" = "LUAD","1" = "BRCA","9" = "UCEC","4" = "STAD","2" = "CRC","7" = "PANCAN-I","3" = "PANCAN-II")
  kmres.rename = clus.map[match(as.character(kmres),names(clus.map))]
  
  # generate cancer specific annotation 
  ln = unique(mut.data$Pathway)
  lncol = colorRampPalette(brewer.pal(9,"Set1"))(length(ln));names(lncol) = ln
  hb = rowAnnotation('Cancer-Specificity' = vec,
                     'Driver Mutation' = mut.matrix,
                     'Driver Marks' = anno_mark(at = at,labels = labs),
                     col = list('Cancer-Specificity' = cancer.col,'Driver Mutation' = lncol))
  
  # generate surface protein & top genbe annotation
  surf.gene = rownames(hub.neglogmat)[rownames(hub.neglogmat) %in% surf.tbl$`UniProt gene`]
  surf.mark = surf.gene
  
  genes = c()
  for (i in 1:length(lincs.lst.surface))
  {
    ii = rownames(hub.neglogmat)[which(hub.neglogmat[,names(lincs.lst.surface)[i]] > -log10(0.05) & rownames(hub.neglogmat) %in% lincs.lst.surface[[i]])]
    genes = c(genes,ii)
    rm(ii)
  }
  surf.mark[surf.gene %in% genes] = paste0(surf.mark[surf.gene %in% genes],"*")
  
  ha = rowAnnotation(
    Surfaceosome = anno_mark(at = which(rownames(hub.neglogmat) %in% surf.gene), labels = surf.mark))
  
  # generate main heatmap
  col_fun = colorRamp2(c(0,-log10(0.05),-log10(0.01),3.5),c("green","white","pink","red"))
  ht = Heatmap(matrix = hub.neglogmat,col = col_fun,name = "-log10(Hub p-value)",
               use_raster = TRUE,show_row_names = FALSE,right_annotation = ha,left_annotation = hb,row_split = factor(kmres.rename))
  
  pdf(file = paste0(out.dir,"/hub_pvalue.pdf"),width = 7,height = 13)
  draw(ht,show_heatmap_legend = FALSE,show_annotation_legend = FALSE)
  dev.off()

}

############################### Pan-cancer regulator nomination: Heatmap plotting
if (TRUE)
{
  library(ComplexHeatmap)
  library(circlize)
 
  if (TRUE)
  {
    ##### load KDA results
    pancan.res = read.delim(file = "Data/KD_stratification/summarized_KDA_results.txt",sep= "\t",header = TRUE,stringsAsFactors = FALSE)
    pancan.res = pancan.res[order(pancan.res$total.count.score,decreasing = TRUE),]
    
    # call unique driver
    pmat = as.matrix(pancan.res[,grep("cFETp_ALL",colnames(pancan.res))]);colnames(pmat) = gsub("cFETp_ALL\\.","",colnames(pmat))
    hmat = as.matrix(pancan.res[,grep("hubcall",colnames(pancan.res))])
    comb.mat = (pmat < 0.05) * (hmat != "");colnames(comb.mat) = colnames(pmat)
    pancan.res$is.unique = (rowSums(pmat < 0.05) == 1 & rowSums(hmat != "") == 1 & rowSums(comb.mat) == 1)
    pancan.res$is.unique.cancer.driver = rep(NA,nrow(pancan.res))
    vec = apply(comb.mat,1,function(x) paste(names(x)[which(x != 0)],collapse = ","))
    pancan.res$is.unique.cancer.driver[pancan.res$is.unique] = vec[pancan.res$is.unique]
    
    # get top results
    top.res = subset(pancan.res,hub.count > 3 & (UP.hit_count > 4 | DN.hit_count > 4) & is.signature.count > 3 & sig.count > 3)
    #top.res = subset(pancan.res,hub.count > 3 & sig.count > 3)
    print(dim(top.res))
   
    # extract enrichment heatmaps
    FET.mat = as.matrix(top.res[,c(grep("cFETp_UP",colnames(top.res)),grep("cFETp_DN",colnames(top.res)))]);
    rownames(FET.mat) = top.res$gene
    FET.mat[FET.mat < 1E-320] = 1E-320
    neglog.mat = -log10(FET.mat);
    neglog.up = neglog.mat[,grep("cFETp_UP",colnames(neglog.mat))];
    colnames(neglog.up) = gsub("cFETp_UP\\.","",colnames(neglog.up))
    neglog.dn = neglog.mat[,grep("cFETp_DN",colnames(neglog.mat))];
    colnames(neglog.dn) = gsub("cFETp_DN\\.","",colnames(neglog.dn))
    
    # annotation data
    col.annt = data.frame(DEP.ID = gsub("^(.*)_|\\.(.*)$","",colnames(FET.mat)),cancer = gsub("^(.*)\\.","",colnames(FET.mat)))
    row.annt = top.res[,c(grep("CRISPRi_hit_rates",colnames(top.res)))]
    colnames(row.annt) = gsub("_hit_rates","",colnames(row.annt))
    colnames(row.annt) = gsub("CRISPRi\\.","",colnames(row.annt))
    hit.annt = top.res[,c("UP.hit_count","DN.hit_count")];colnames(hit.annt) = c("DEP-UP","DEP-DN")
    hub.annt = top.res[,grep("hubcall",colnames(top.res))];colnames(hub.annt) = gsub("hubcall\\.","",colnames(hub.annt))
    hub.annt = as.matrix(hub.annt)
    hub.annt = t(apply(hub.annt,1,function(x) c(sum(x == "global",na.rm = TRUE),sum(x == "sub-global",na.rm = TRUE),sum(x == "local",na.rm = TRUE))))
    ##### fabricate protein annotation bar
    # specify colormap: CRISPRi features
    gene.annot = rowAnnotation(CRISPRi = as.matrix(row.annt),col = list(CRISPRi = colorRamp2(breaks = c(0,1),colors = c("white","darkorchid2"))))
    
    # KD hit by signatures
    hit.annot = rowAnnotation('KD' = anno_barplot(hit.annt,gp = gpar(fill = c("red","blue"))),
                              'Hub' = anno_barplot(hub.annt,gp = gpar(fill = c("red","coral1","darkgoldenrod1"))))
    
    # cluster genes
    rowh = hclust(dist(neglog.mat),"complete")
    
    ### make main heatmap
    heat.up = Heatmap(matrix = neglog.up,name = "-log10[cFET p(DEP-UP)]",
                      col = colorRamp2(breaks = c(0,2,5,10,max(neglog.mat,na.rm = TRUE)*1.01),colors = c("white","white","yellow","red","red")),
                      cluster_rows = rowh,left_annotation = hit.annot
    )
    heat.dn = Heatmap(matrix = neglog.dn,name = "-log10[cFET p(DEP-DN)]",
                      col = colorRamp2(breaks = c(0,2,5,10,max(neglog.mat,na.rm = TRUE)*1.01),colors = c("white","white","cyan","blue","blue")),
                      cluster_rows = rowh,
                      right_annotation = gene.annot)
    heat.comb= heat.up + heat.dn
    
    # make legends
    up.lgd = Legend(col_fun = colorRamp2(breaks = c(0,2,5,10),colors = c("white","white","yellow","red")),
                    title = "-log10[cFET p(DEP-UP)]",direction = "horizontal",labels_rot = 45)
    dn.lgd = Legend(col_fun = colorRamp2(breaks = c(0,2,5,10),colors = c("white","white","cyan","blue")),
                    title = "-log10[cFET p(DEP-DN)]",direction = "horizontal",labels_rot = 45)
    rate.lgd = Legend(col_fun = colorRamp2(breaks = c(0,1),colors = c("white","darkorchid2")),
                      title = "CRISPRi hit rate",direction = "horizontal",labels_rot = 45)
    hub.lgd = Legend(labels = c("global","sub-global","local"),title = "Hub Type",legend_gp = gpar(fill = c("red","coral1","darkgoldenrod1")),ncol = 1)
    KD.lgd = Legend(labels = c("DEP-UP","DEP-DN"),title = "KD Type",legend_gp = gpar(fill = c("red","blue")),ncol = 1)
    
    pdf(file = paste0(out.dir,"/PanCancer_Drivers.pdf"),width = 8,height = 13)
    draw(heat.comb,show_heatmap_legend = FALSE,
         annotation_legend_list = list(up.lgd,dn.lgd,rate.lgd,hub.lgd,KD.lgd),annotation_legend_side = "bottom")
    dev.off()
  }
}

library(VennDiagram)
library(gridExtra)
library(lattice)
library(gtable)

require(reshape2)
require(ggplot2)
require(ggrepel)

load_preservation_data <- function(fd)
{
  prvfiles = list.files(path = fd,pattern = "GEX\\.results\\.txt",full.names = TRUE)
  names(prvfiles) = gsub("_(.*)$","",gsub("^(.*)/","",prvfiles))
  #prvfiles = prvfiles[cancers]
  prv.lst = lapply(prvfiles,function(x) read.delim(file = x,sep = "\t",header = TRUE,stringsAsFactors = FALSE))
  return(prv.lst)
}

generate_pres_volcano <- function(tbl1,tbl2,top.n = 5)
{
  mid = intersect(tbl1[[1]],tbl2[[1]])
  tbl1 = tbl1[match(mid,tbl1[[1]]),];tbl2 = tbl2[match(mid,tbl2[[1]]),]
  tbl1[[1]] = gsub("c1_","M",tbl1[[1]]);tbl2[[1]] = gsub("c1_","M",tbl2[[1]])
  tbl = rbind.data.frame(tbl1,tbl2)
  valid.modules = unique(subset(tbl,log.p.Bonfsummary.pres <= -10)$Module)
  spec.modules = intersect(subset(tbl1,Zsummary.pres <= 2)$Module,subset(tbl2,Zsummary.pres <= 2)$Module)
  cancer = unique(tbl$cancer)
  
  zmat = acast(data = tbl,formula = Module ~ cohort,value.var = "Zsummary.pres",fun.aggregate = function(x) max(x,na.rm = TRUE) )
  zmat[is.infinite(zmat)] = 0
  vec = rep(NA,nrow(zmat));vec[rownames(zmat) %in% valid.modules] = "Preserved";vec[rownames(zmat) %in% spec.modules] = "Proteome-specific";
  top.mods = rownames(zmat)[order(rowSums(zmat,na.rm = TRUE),decreasing = TRUE)[1:top.n]]
  top.spec = rownames(zmat)[order(rowSums(zmat,na.rm = TRUE),decreasing = FALSE)[1:top.n]]
  zdata = data.frame(Module = rownames(zmat),as.data.frame(zmat),is.significant = vec)
  pobj = ggplot() + 
    geom_point(data = zdata,aes_string(x = colnames(zmat)[1],y = colnames(zmat)[2],colour = "is.significant"),alpha = 0.5) + 
    geom_text_repel(data = subset(zdata,Module %in% c(top.mods)),aes_string(x = colnames(zmat)[1],y = colnames(zmat)[2],label = "Module"),size = 8) + 
    scale_colour_manual(values = c("Preserved" = "brown","Proteome-specific" = "chartreuse1"),na.value = "grey") + 
    guides(colour = guide_legend(title = "-log10(Bonf.P) < 10")) + 
    labs(title = cancer,x = paste0("Z-score:",colnames(zmat)[1]),y = paste0("Z-score:",colnames(zmat)[2])) + 
    theme_bw() +
    scale_x_sqrt() + scale_y_sqrt() + 
    theme(axis.text = element_text(size = 17),axis.title = element_text(size = 18),plot.title = element_text(size = 20,hjust = 0.5),
          legend.title = element_text(size = 16),legend.text = element_text(size = 14),legend.position = "bottom")
  return(pobj)
}

make_venn_diagrams <- function(tbl1,tbl2)
{
  # format tables
  mid = intersect(tbl1[[1]],tbl2[[1]])
  tbl1 = tbl1[match(mid,tbl1[[1]]),];tbl2 = tbl2[match(mid,tbl2[[1]]),]
  tbl1[[1]] = gsub("c1_","M",tbl1[[1]]);tbl2[[1]] = gsub("c1_","M",tbl2[[1]])
  tbl = rbind.data.frame(tbl1,tbl2)
  
  # call significant calls
  vec = rep(NA,nrow(tbl));vec[tbl$log.p.Bonfsummary.pres <= -10] = "Preserved";vec[tbl$Zsummary.pres <= 2] = "Proteome-specific"
  tbl$Preservation.call = vec;
  
  # make preservation diagram
  cohorts = unique(tbl$cohort)
  prv.lst = lapply(cohorts,function(x,y) subset(y,Preservation.call == "Preserved" & cohort == x)$Module,y = tbl)
  #names(prv.lst) = cohorts
  names(prv.lst) = rep("",length(cohorts))
  spec.lst = lapply(cohorts,function(x,y) subset(y,Preservation.call == "Proteome-specific" & cohort == x)$Module,y = tbl)
  #names(spec.lst) = cohorts
  names(spec.lst) = rep("",length(cohorts))
  
  prv.grob = venn.diagram(x = prv.lst,filename = NULL,cat.cex = 2,cex = 2,fill = c("red","chartreuse1"),alpha = c(0.2,0.2))
  spec.grob = venn.diagram(x = spec.lst,filename = NULL,cat.cex = 2,cex = 2,fill = c("red","chartreuse1"),alpha = c(0.2,0.2))
  return(list("preserved" = prv.grob,"specific" = spec.grob))
}


library(MEGENA)
generate_sunburst <- function(tbl1,tbl2,mtbl)
{
  ### process preservation table
  if (!is.null(tbl1) & !is.null(tbl2))
  {
    mid = intersect(tbl1[[1]],tbl2[[1]])
    tbl1 = tbl1[match(mid,tbl1[[1]]),];tbl2 = tbl2[match(mid,tbl2[[1]]),]
    tbl1[[1]] = gsub("c1_","M",tbl1[[1]]);tbl2[[1]] = gsub("c1_","M",tbl2[[1]])
    tbl = rbind.data.frame(tbl1,tbl2)
    valid.modules = intersect(subset(tbl1,log.p.Bonfsummary.pres <= -10)$Module,subset(tbl2,log.p.Bonfsummary.pres <= -10)$Module)
    spec.modules = intersect(subset(tbl1,Zsummary.pres <= 2)$Module,subset(tbl2,Zsummary.pres <= 2)$Module)
    cancer = unique(tbl$cancer)
    
    zmat = acast(data = tbl,formula = Module ~ cohort,value.var = "Zsummary.pres",fun.aggregate = function(x) max(x,na.rm = TRUE) )
    zmat[is.infinite(zmat)] = 0
    vec = rep(NA,nrow(zmat));vec[rownames(zmat) %in% valid.modules] = "Preserved";vec[rownames(zmat) %in% spec.modules] = "Proteome-specific";
    zdata = data.frame(Module = rownames(zmat),as.data.frame(zmat),is.significant = vec)
  }else{
    tbl1[[1]] = gsub("c1_","M",tbl1[[1]]);
    vec = rep(NA,nrow(tbl1));
    vec[tbl1$log.p.Bonfsummary.pres <= -10] = "Preserved";
    vec[tbl1$Zsummary.pres <= 2] = "Proteome-specific"
    zdata = tbl1[,c("Module","Zsummary.pres","log.p.Bonfsummary.pres")]
    zdata$is.significant = vec
  }
  
  # update module table
  mtbl$GEX.Preservation.class = factor(zdata$is.significant[match(mtbl$id,zdata$Module)])
  cancer.id = unique(mtbl$group)
  # get sunburst 
  sbobj = draw_sunburst_wt_fill(module.df = mtbl,
                                parent.col = "module.parent",id.col = "id",
                                min.angle = Inf,
                                feat.col = "GEX.Preservation.class",
                                fill.type = "discrete",log.transform = FALSE,
                                fill.scale = scale_fill_manual(values = c("Preserved" = "brown","Proteome-specific" = "chartreuse1"),na.value = alpha("grey",0.5)),
                                theme.adjust = NULL
  )
  sbobj = sbobj + 
    geom_text(x = 0,y = 0,label = cancer.id,size = 4) +
    theme(plot.title = element_text(hjust = 0.5,size = 18),
          panel.border = element_blank(),panel.background = element_blank()
    ) + guides(fill = FALSE)
  return(list(sb = sbobj,mtbl = mtbl))
}
rm(list = ls())

source("scripts/R_functions/preservation_analysis_functions.v1.R")
#################

out.dir = "Transcriptome_Presevation_Results";dir.create(out.dir)

###### load results
res = list(
  "PanCanAtlas" = load_preservation_data("Data/Transcriptome_Preservation/Proteome_vs_Transcriptome_PanCan"),
  "CPTAC" = load_preservation_data("Data/Transcriptome_Preservation/Proteome_vs_Transcriptome_CPTAC")
)

for (i in 1:length(res))
{
  for (j in 1:length(res[[i]]))
  {
    res[[i]][[j]]$cohort = names(res)[i]
    res[[i]][[j]]$cancer = names(res[[i]])[j]
  }
}

# match cancer types
#cancers = Reduce("intersect",lapply(res,names))
cancers = c("BRCA","CCRCC","CRC","LUAD","UCEC")
res = lapply(res,function(x,y) x[y],y = cancers)

###### Cohort comparison: See if CPTAC and TCGA transcriptome yield similar results via the volcanoe plot
require(reshape2)
require(ggplot2)
require(ggrepel)

vlst = mapply(FUN = function(x,y) generate_pres_volcano(tbl1 = x,tbl2 = y),x = res[[1]],y = res[[2]],SIMPLIFY = FALSE)
for (i in 1:length(vlst)) vlst[[i]] = vlst[[i]] + guides(colour = FALSE)

# output .png
require(grid)
png(file = paste0(out.dir,"/Cohort_Comparison.png"),res = 350,width = 6400,height = 1500)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,5)))
for (i in 1:length(vlst))
{
  print(vlst[[i]], vp = viewport(layout.pos.row = 1,layout.pos.col = i))
}
dev.off()

# output .pdf
pdf(file = paste0(out.dir,"/Cohort_Comparison.pdf"),width = 20,height = 4)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,5)))
for (i in 1:length(vlst))
{
  print(vlst[[i]] + theme(), vp = viewport(layout.pos.row = 1,layout.pos.col = i))
}
dev.off()

####### Now, generate venn diagram
library(VennDiagram)
library(gridExtra)
library(lattice)
library(gtable)

# create /png figues
if (TRUE)
{
  for (i in 1:length(cancers))
  {
    tbl1 = res[[1]][[i]]
    tbl2 = res[[2]][[i]]
    
    vres = make_venn_diagrams(tbl1,tbl2)
    png(file = paste0(out.dir,"/",cancers[i],".venn_diagram.png"),res = 350,width = 1400,height = 2000)
    grid.newpage()
    grid.arrange(gTree(children=vres[[1]]), gTree(children=vres[[2]]),ncol=1,nrow = 2)
    dev.off()
  }
  
  # get legend
  png(file = paste0(out.dir,"/venn_legend.png"),res = 350,width = 2500,height = 1500)
  plot.new()
  legend("center",legend = c("PanCanAtlas TX","CPTAC TX"),fill = alpha(c("red","chartreuse1"),0.3),cex = 3)
  dev.off()
}

# create .pdf figure
if (TRUE)
{
  vdat = list()
  for (i in 1:length(cancers))
  {
    tbl1 = res[[1]][[i]]
    tbl2 = res[[2]][[i]]
    
    vres = make_venn_diagrams(tbl1,tbl2);
    vdat = c(vdat,vres)
    rm(vres)
  }
  pdf(file = paste0(out.dir,"/venn_diagram.pdf"),width = 20,height = 8)
  #grid.newpage()
  grid.arrange(gTree(children=vdat[[1]]),gTree(children=vdat[[3]]),gTree(children=vdat[[5]]), gTree(children=vdat[[7]]), gTree(children=vdat[[9]]),
               gTree(children=vdat[[2]]),gTree(children=vdat[[4]]), gTree(children=vdat[[6]]),gTree(children=vdat[[8]]), gTree(children=vdat[[10]]),
               ncol=5,nrow = 2)
  dev.off()
  
  pdf(file = paste0(out.dir,"/venn_legend.pdf"),width = 6,height = 3)
  plot.new()
  legend("center",legend = c("PanCanAtlas TX","CPTAC TX"),fill = alpha(c("red","chartreuse1"),0.3),cex = 3,border = "white",bty = "n")
  dev.off()
}


######## Now, generate sunburst plots
###### load results
res = list(
  "PanCanAtlas" = load_preservation_data("Data/Transcriptome_Preservation/Proteome_vs_Transcriptome_PanCan"),
  "CPTAC" = load_preservation_data("Data/Transcriptome_Preservation/Proteome_vs_Transcriptome_CPTAC")
)

for (i in 1:length(res))
{
  for (j in 1:length(res[[i]]))
  {
    res[[i]][[j]]$cohort = names(res)[i]
    res[[i]][[j]]$cancer = names(res[[i]])[j]
  }
}

# match cancer types
cancers = c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")
res = lapply(res,function(x,y) x[y],y = cancers)

# get module table per cancer
modtbl = read.delim(file = "Data/MEGENA/multiscale_module_summary.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
modtbl = subset(modtbl,group %in% cancers)
modtbl$id = gsub("c1_","M",modtbl$id)
modtbl$module.parent = gsub("c1_","M",modtbl$module.parent)
modtbl.lst = lapply(split(1:nrow(modtbl),factor(modtbl$group)),function(x,y) y[x,],y = modtbl)
for (i in 1:length(modtbl.lst))
{
  modtbl.lst[[i]]$id = gsub("\\|(.*)$","",modtbl.lst[[i]]$id)
  modtbl.lst[[i]]$module.parent = gsub("\\|(.*)$","",modtbl.lst[[i]]$module.parent)
}

sb.res = mapply(FUN = function(tbl1,tbl2,mtbl) generate_sunburst(tbl1,tbl2,mtbl),tbl1 = res[[1]],tbl2 = res[[2]],mtbl = modtbl.lst,SIMPLIFY = FALSE)
plot.list = lapply(sb.res,function(x) x$sb)

# plot sunbursts
require(grid)
png(file = paste0(out.dir,"/Transcriptome_Preservation_Sunburst.png"),res = 350,width = 3000,height = 1500)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,4)))

for (i in 1:length(plot.list))
{
  col = ceiling(i/4);row = i - 4*(col-1)
  print(plot.list[[i]], vp = viewport(layout.pos.row = col,layout.pos.col = row))
}

dev.off()

pdf(file = paste0(out.dir,"/Transcriptome_Preservation_Sunburst.pdf"),width = 11,height = 5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,4)))

for (i in 1:length(plot.list))
{
  col = ceiling(i/4);row = i - 4*(col-1)
  print(plot.list[[i]], vp = viewport(layout.pos.row = col,layout.pos.col = row))
}

dev.off()


# get barplot summarizing proportions
mat = matrix(0,nrow = 3,ncol = length(cancers));colnames(mat) = cancers;rownames(mat) = c("Preserved","Proteome-specific","NA")
for (i in 1:length(sb.res))
{
  vec = sb.res[[i]]$mtbl$GEX.Preservation.class
  mat[1,i] = sum(vec == "Preserved",na.rm = TRUE)
  mat[2,i] = sum(vec == "Proteome-specific",na.rm = TRUE)
  mat[3,i] = sum(is.na(vec))
  rm(vec)
}

# plot fractions
require(reshape)
fmat = apply(mat,2,function(x) x/sum(x));rownames(fmat)[3] = "Unassigned"
pdata = melt(fmat[1:3,]);colnames(pdata) = c("Category","cancer","proportion")
pobj = ggplot(data = pdata,aes(x = cancer,y = proportion,fill = Category)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values =  c("Preserved" = "brown","Proteome-specific" = "chartreuse1","Unassigned" = "grey")) + theme_minimal() + 
  guides(fill = guide_legend(nrow = 2)) + 
  theme(legend.position = "bottom",legend.direction = "horizontal",axis.text = element_text(size = 16),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 18),
        axis.title = element_text(size = 22),legend.title = element_text(size = 17),legend.text = element_text(size = 13))

png(file = paste0(out.dir,"/Transcriptome_Preservation_Barplot.png"),res = 350,width = 2200,height = 1700)
print(pobj)
dev.off()

pdf(file = paste0(out.dir,"/Transcriptome_Preservation_Barplot.pdf"),width = 6,height = 5)
print(pobj + guides(fill = "none") + theme(axis.title.x = element_blank()))
dev.off()

pdf(file = paste0(out.dir,"/barplot_legend.pdf"),width = 21,height = 3)
plot.new()
legend("center",legend = c("Preserved","Proteome-specific","Unassigned"),fill = c("brown","chartreuse1","grey"),cex = 3,border = "white",bty = "n",ncol=3)
dev.off()

################################# Now summarize MSigDB terms per category
cancers = c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")

# get enriched term
get_mod_term <- function(modtbl)
{
  term.mat = as.matrix(modtbl[,grep("^SigTerms",colnames(modtbl))]);rownames(term.mat) = modtbl$id;
  term.mx = apply(term.mat,1,function(x) {
    dat = do.call('rbind',strsplit(do.call('c',strsplit(x,",")),"/"))
    dat[which.min(as.numeric(dat[,2])),1]
  })
  return(term.mx)
}

rename_terms = function(x)
{
  out = gsub("^BIOCARTA_|^KEGG_|^REACTOME_|^HALLMARK_","",x)
  out = gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","EMT",out)
  out = gsub("ENDOPLASMIC_RETICULUM","ER",out)
  out = gsub("EXTRACELLULAR_MATRIX","ECM",out)
  out = gsub("EXTRACELLULAR","EC",out)
  return(out)
}

format_name = function(x,nmx = 20) 
{
  str = strsplit(x,"_")[[1]]
  nvec = cumsum(sapply(str,nchar))
  nvec[1:(length(nvec)-1)] = nvec[1:(length(nvec)-1)] + 1
  nline = ceiling(nvec/nmx)
  str.n = split(str,factor(nline))
  str.n = sapply(str.n,function(x) paste(x,collapse = "_"))
  if (length(str.n) > 1) str.n[1:(length(str.n)-1)] = paste(str.n[1:(length(str.n)-1)],"\n",sep = "")
  paste(str.n,collapse = "")
}

##############################
filenames = list.files(path = "Data/MEGENA/ModuleRanking",pattern = "\\.ranked_module_table\\.txt$",full.names = TRUE)
names(filenames) = gsub("^(.*)/|\\.(.*)$","",filenames)
filenames = filenames[cancers]

## collect list of interesting modules
mtbl.list = vector("list",length(filenames));names(mtbl.list) = names(filenames)
sig.list = list("DEP-UP" = c("DEP_UP","EBV_UP"),
                "DEP-DN" = c("DEP_DN","EBV_DN"))
for (i in 1:length(filenames))
{
  tbl = read.delim(file = filenames[i],sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  tbl$module.term = rename_terms(get_mod_term(tbl))
  
  # get signatue score
  tbl$UP.score = apply(cbind(tbl[,colnames(tbl) %in% sig.list$`DEP-UP`]),1,function(x) sum(-log10(x)))
  tbl$DN.score = apply(cbind(tbl[,colnames(tbl) %in% sig.list$`DEP-DN`]),1,function(x) sum(-log10(x)))
  tbl$id = gsub("\\|(.*)$","",gsub("^c1_","M",tbl$id))
  tbl$module.parent = gsub("\\|(.*)$","",gsub("^c1_","M",tbl$module.parent))
  
  tbl$GEX.Preservation.class = sb.res[[i]]$mtbl$GEX.Preservation.class[match(tbl$id,sb.res[[i]]$mtbl$id)]
  
  stbl = subset(tbl,!is.na(GEX.Preservation.class))
  mtbl.list[[i]] = stbl
  rm(tbl,stbl)
}

## collect MSigDB results
msigdb.dirs = list.dirs(path = "Data/MEGENA/Individual_MEGENA_Proteome")
msigdb.dirs = msigdb.dirs[grep("/MSigDB",msigdb.dirs)]
names(msigdb.dirs) = gsub("\\.(.*)$","",gsub("^(.*)__","",msigdb.dirs))
msigdb.dirs = msigdb.dirs[cancers]

all.df = data.frame()
for (i in 1:length(msigdb.dirs))
{
  fetfiles = list.files(path = msigdb.dirs[i],pattern = "_FET-Table.txt$",full.names = TRUE)
  names(fetfiles) = gsub("^(.*)/|_FET-Table.txt$","",fetfiles)
  #fetfiles = fetfiles[c("c2.cp.kegg.v5.0.symbols","c5.bp.v5.0.symbols","c5.mf.v5.0.symbols","h.all.v5.0.symbols")]
  #fetfiles = fetfiles["h.all.v5.0.symbols"]
  for (j in 1:length(fetfiles))
  {
    df = read.delim(file = fetfiles[j],sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    df$database = rep(names(fetfiles)[j],nrow(df))
    df$cancer = rep(names(msigdb.dirs)[i],nrow(df))
    df$module.id = gsub("^c1_","M",df$set1_Name)
    df = subset(df,module.id %in% mtbl.list[[i]]$id)
    df$conservation.call = mtbl.list[[i]]$GEX.Preservation.call[match(df$module.id,mtbl.list[[i]][[1]])]
    all.df = rbind.data.frame(all.df,df)
    rm(df)
  }
}

######### Create ranking of MSigDB terms best represented in a group of modules
require(reshape2)
require(ggrepel)
require(ggbeeswarm)
top.n = 2

all.term.data = all.pathways = data.frame();
for (preserve.call in c("Preserved","Proteome-specific"))
{
  ### summarize up-regulated DEPs
  if (TRUE)
  {
    term.data = data.frame();
    for (mi in 1:length(mtbl.list))
    {
      mid = subset(mtbl.list[[mi]],GEX.Preservation.class == preserve.call & (UP.score > -log10(0.05) ))[[1]]
      
      if (length(mid) > 0)
      {
        mod.df = subset(all.df,module.id %in% mid & cancer == names(mtbl.list)[mi])
        mod.df = mod.df[order(mod.df$FET_pvalue),]
        mod.df$preservation.call = rep(preserve.call,nrow(mod.df))
        mod.df$DEP.call = rep("DEP-UP",nrow(mod.df))
        all.pathways = rbind.data.frame(all.pathways,mod.df)
        
        term.mat = acast(data = mod.df,formula = module.id ~ set2_Name,value.var = "corrected.FET.pvalue",fun.aggregate = function(x) min(x,na.rm = TRUE))
        #term.mat = acast(data = mod.df,formula = module.id ~ set2_Name,value.var = "FET_pvalue",fun.aggregate = function(x) min(x,na.rm = TRUE))
        term.mat[is.infinite(term.mat)] = 1
        
        #term.score = colSums(rbind(-log10(term.mat)))
        term.score = apply(rbind(-log10(term.mat)),2,function(x) max(x,na.rm = TRUE))
        term.score = term.score[order(term.score,decreasing = TRUE)]
        is.top = rep(FALSE,length(term.score))
        is.top[1:top.n] = TRUE
        is.top[term.score < -log10(0.05)] = FALSE
        
        term.data = rbind.data.frame(term.data,data.frame(term.name = names(term.score),term.score = term.score,cancer = rep(names(mtbl.list)[mi],length(term.score)),
                                                          is.top = is.top))
        
      }
      
      rm(mid)
    }
    term.data$preservation.call = rep(preserve.call,nrow(term.data))
    term.data$DEP.call = rep("DEP-UP",nrow(term.data))
    term.data$term.short = rename_terms(term.data$term.name) 
  }
  all.term.data = rbind.data.frame(all.term.data,term.data)
  
  ### summarize down-regulated DEPs
  if (TRUE)
  {
    term.data = data.frame();
    for (mi in 1:length(mtbl.list))
    {
      mid = subset(mtbl.list[[mi]],GEX.Preservation.class == preserve.call & (DN.score > -log10(0.05) ))[[1]]
      
      if (length(mid) > 0)
      {
        mod.df = subset(all.df,module.id %in% mid  & cancer == names(mtbl.list)[mi])
        mod.df = mod.df[order(mod.df$FET_pvalue),]
        mod.df$preservation.call = rep(preserve.call,nrow(mod.df))
        mod.df$DEP.call = rep("DEP-DN",nrow(mod.df))
        all.pathways = rbind.data.frame(all.pathways,mod.df)
        
        term.mat = acast(data = mod.df,formula = module.id ~ set2_Name,value.var = "corrected.FET.pvalue",fun.aggregate = function(x) min(x,na.rm = TRUE))
        #term.mat = acast(data = mod.df,formula = module.id ~ set2_Name,value.var = "FET_pvalue",fun.aggregate = function(x) min(x,na.rm = TRUE))
        term.mat[is.infinite(term.mat)] = 1
        
        #term.score = colSums(rbind(-log10(term.mat)))
        term.score = apply(rbind(-log10(term.mat)),2,function(x) max(x,na.rm = TRUE))
        term.score = term.score[order(term.score,decreasing = TRUE)]
        is.top = rep(FALSE,length(term.score))
        is.top[1:top.n] = TRUE
        is.top[term.score < -log10(0.05)] = FALSE
        
        term.data = rbind.data.frame(term.data,data.frame(term.name = names(term.score),term.score = term.score,cancer = rep(names(mtbl.list)[mi],length(term.score)),
                                                          is.top = is.top))
        
      }
      
      rm(mid)
    }
    term.data$preservation.call = rep(preserve.call,nrow(term.data))
    term.data$DEP.call = rep("DEP-DN",nrow(term.data))
    term.data$term.short = rename_terms(term.data$term.name) 
  }
  all.term.data = rbind.data.frame(all.term.data,term.data)
  
}

all.term.data$term.reformat = sapply(all.term.data$term.short,format_name)
#all.pathways$term.reformat = sapply(all.pathways$set2_Name,format_name)

# plot out prioritized pathways: preliminary plot
cancer.col = c("BRCA" = "red","CCRCC" = "orange","CRC" = "thistle4","HCC" = "darkgreen",
               "LUAD" = "blue","MB" = "deepskyblue",
               "STAD" = "slateblue","UCEC" = "darkturquoise","PANCAN" = "black");

pobj = ggplot() + geom_quasirandom(data = all.term.data,aes(x = cancer,y = term.score,colour = cancer),alpha = 0.5) + 
  geom_label_repel(data = subset(all.term.data,is.top),aes(x = cancer,y = term.score,label = term.reformat,colour = cancer),alpha = 0.8,angle = 45) + 
  labs(x = "Cancers",y = "-log10(FET FDR)") + facet_grid(DEP.call + preservation.call ~ .,scale = "free_y") + 
  geom_hline(yintercept = -log10(0.05),colour = "red") + 
  scale_colour_manual(values = cancer.col[!(names(cancer.col) %in% c("MB","PANCAN"))]) +
  theme_bw() + theme(axis.title = element_text(size = 25),axis.text = element_text(size = 20),
                     strip.text = element_text(size = 20),
                     legend.position = "bottom",legend.direction = "horizontal")

pdf(file = paste0(out.dir, "/Transcriptome_Comparison_MSigDB_Terms.pdf"),width = 16,height = 14)
print(pobj)
dev.off()



rm(list = ls())

source("scripts/R_functions/MoM_functions.v1.R")

####################

out.dir = "Data/CrossCancer_Module_Preservation"
res.dir = "Conserved_Module_Networks";dir.create(res.dir)
thresh.val = -10;

modfile = "Data/MEGENA/multiscale_significant.modules.txt"
netfile = "Data/MEGENA/MEGENA_Network.txt"
outfile = "CrossCancer_Module_Preservation_Summary.txt"

cancers = c("BRCA","CCRCC","CRC","LUAD","HCC","STAD","UCEC")

############################### get Preservation data 
if (TRUE)
{
  resfiles = list.files(path = out.dir,pattern = "\\.results\\.txt$",full.names = TRUE)
  names(resfiles) = gsub("^(.*)/|\\.results\\.txt$","",resfiles)
  
  all.res = data.frame()
  for (i in 1:length(resfiles))
  {
    res = read.delim(file = resfiles[i],sep = "\t",header = TRUE,stringsAsFactors = FALSE)
    res$comparison.id = rep(names(resfiles)[i],nrow(res))
    res$case.id = rep(gsub("_vs_(.*)$","",names(resfiles)[i]),nrow(res))
    res$ctrl.id = rep(gsub("^(.*)_vs_","",names(resfiles)[i]),nrow(res))
    
    # call preserved modules across cancers
    res$is.conserved = res$log.p.Bonfsummary.pres <= thresh.val
    res$Module = paste(res$Module,res$case.id,sep = "|")
    
    # call unpreserved 
    res$is.not.conserved = res$Zsummary.pres < 2 #& res$log.p.Bonfsummary.pres > thresh.val & res$Zconnectivity.pres < 2 & res$Zsummary.qual < 2
    all.res = rbind.data.frame(all.res,res);
    rm(res)
  }
  
  write.table(all.res,file = outfile,
              sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
}else{
  all.res = read.delim(file = outfile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
}

# subset for seven cancers
all.res = subset(all.res,case.id %in% cancers & ctrl.id %in% cancers)

############################### get module overlap
fetfile = "Data/CrossCancer_Module_Overlap/Overlap_Results.signifFET.txt"
all.FET.call = read.delim(file = fetfile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
all.FET.call = subset(all.FET.call,corrected.FET.pvalue < 0.05)

####### test FET thresholds to call preserved pairs of modules
modules = MEGENA::read.geneSet(modfile);
modules = modules[gsub("^(.*)\\|","",names(modules)) %in% cancers]
max.size = 3000;min.size = 10;
overlap.cutoff = seq(0.1,0.55,0.01)
efc.cut = 1;
pval.cut = 0.05
excl.score = tot.sizes = rep(NA,length(overlap.cutoff))
mod.sizes = vector("list",length(overlap.cutoff))
for (i in 1:length(overlap.cutoff))
{
  cat(paste("processing:",overlap.cutoff[i],"\n",sep = ""))
  ### step 0: tag pairs by cutoff
  all.signif.call = assign_overlap_call(atbl = all.FET.call,pval.cut = pval.cut,efc.cut = efc.cut,overlap.cut = overlap.cutoff[i])
  all.signif.call$set1_Name.full = paste(all.signif.call$set1_Name,all.signif.call$cancer.i,sep = "|")
  all.signif.call$set2_Name.full = paste(all.signif.call$set2_Name,all.signif.call$cancer.j,sep = "|")
  
  ### step 1. Update conservation call from Module Presrvation analysis
  require(reshape2)
  cons.mat = acast(data = all.res,formula = Module ~ ctrl.id,value.var = "is.conserved",fun.aggregate = function(x) any(x))
  
  # update preservation call 
  ij = cbind(match(all.signif.call$set1_Name.full,rownames(cons.mat)),match(all.signif.call$cancer.j,colnames(cons.mat)))
  cons.call = apply(ij,1,function(x,m) {out = NA;if (all(!is.na(x))) out = m[x[1],x[2]];return(out)},m = cons.mat)
  all.signif.call$set1_preserved = cons.call
  
  ij = cbind(match(all.signif.call$set2_Name.full,rownames(cons.mat)),match(all.signif.call$cancer.i,colnames(cons.mat)))
  cons.call = apply(ij,1,function(x,m) {out = NA;if (all(!is.na(x))) out = m[x[1],x[2]];return(out)},m = cons.mat)
  all.signif.call$set2_preserved = cons.call

  ### step 2. get finalized pairs of preserved moduels
  all.preserved = subset(all.signif.call,set1_preserved & set2_preserved & conservation.class %in% c("both.conserved"))
  #all.preserved = subset(all.signif.call,set1_preserved & set2_preserved & overall_overlap > 0.2 )
  print(dim(all.preserved))
  
  if (nrow(all.preserved) > 1)
  {
    ### step 3. identify modules of modules (MoMs) through the conserved relationships
    mom.res = identify_MoM(all.preserved,modules)
    cls.list = split(as.character(mom.res$MoM.table$module.id),factor(mom.res$MoM.table$cluster));names(cls.list) = paste("MoM",names(cls.list),sep = "")
    cls.list = cls.list[sapply(cls.list,length) > 1]
    raw.modules = lapply(cls.list,function(x,y) Reduce("union",y[x]),y = modules)
    raw.modules = raw.modules[sapply(raw.modules,length) <= max.size & sapply(raw.modules,length) >= min.size]
    
    ######## evaluation
    modlen = sapply(raw.modules,length)
    tot = length(Reduce("union",raw.modules))
    int.num = rep(NA,length(raw.modules))
    
    for (m in 1:length(raw.modules))
    {
      int.num[m] = length(intersect(raw.modules[[m]],Reduce("union",raw.modules[-m])))
    }
    
    excl.score[i] = mean((modlen-int.num)/modlen)
    tot.sizes[i] = tot
    mod.sizes[[i]] = modlen
  }
  
}

require(ggplot2)
require(ggrepel)
plot.data = data.frame(x = tot.sizes,y = excl.score,overlap.cutoff = overlap.cutoff)

pobj = ggplot(data = plot.data,aes(x = x,y = y)) + geom_point() + 
  geom_text_repel(data = subset(plot.data,excl.score > 0.6),aes(x = x,y = y,label = overlap.cutoff)) + labs(x = "Total union size",y = "# of MoM-specific gene")

pdf(file = paste(res.dir,"/Overlap_Evaluation.pdf",sep = ""))
print(pobj) # gives overlap cutoff = 0.44 as optimal cutoff giving the least overlap across unique MoMs
dev.off()

############################### Now, apply the following criteria: FET FDR p < 0.05 & pairs are preserved in both direction
if (TRUE)
{
  ### step 0: tag pairs by cutoff
  all.signif.call = assign_overlap_call(atbl = all.FET.call,pval.cut = pval.cut,efc.cut = efc.cut,overlap.cut = 0.44)
  all.signif.call$set1_Name.full = paste(all.signif.call$set1_Name,all.signif.call$cancer.i,sep = "|")
  all.signif.call$set2_Name.full = paste(all.signif.call$set2_Name,all.signif.call$cancer.j,sep = "|")
  
  write.table(all.signif.call,paste(res.dir,"/Overlap_Results.annotated_pairs.txt",sep = ""),sep = "\t",row.names = FALSE,
              col.names = TRUE,quote = FALSE)
  
  ### step 1. Update conservation call from Module Presrvation analysis
  require(reshape2)
  cons.mat = acast(data = all.res,formula = Module ~ ctrl.id,value.var = "is.conserved",fun.aggregate = function(x) any(x))
  
  # update preservation call 
  ij = cbind(match(all.signif.call$set1_Name.full,rownames(cons.mat)),match(all.signif.call$cancer.j,colnames(cons.mat)))
  cons.call = apply(ij,1,function(x,m) {out = NA;if (all(!is.na(x))) out = m[x[1],x[2]];return(out)},m = cons.mat)
  all.signif.call$set1_preserved = cons.call
  
  ij = cbind(match(all.signif.call$set2_Name.full,rownames(cons.mat)),match(all.signif.call$cancer.i,colnames(cons.mat)))
  cons.call = apply(ij,1,function(x,m) {out = NA;if (all(!is.na(x))) out = m[x[1],x[2]];return(out)},m = cons.mat)
  all.signif.call$set2_preserved = cons.call
  
  ### step 2. get finalized pairs of preserved moduels
  all.preserved = subset(all.signif.call,set1_preserved & set2_preserved & conservation.class %in% c("both.conserved"))
  
  # identify modules of modules (MoMs) through the conserved relationships
  mom.res = identify_MoM(all.preserved,modules)
  cls.list = split(as.character(mom.res$MoM.table$module.id),factor(mom.res$MoM.table$cluster));names(cls.list) = paste("MoM",names(cls.list),sep = "")
  cls.list = cls.list[sapply(cls.list,length) > 1]
  
  ####
  raw.modules = lapply(cls.list,function(x,y) Reduce("union",y[x]),y = modules)
  sapply(raw.modules,length)
  overlap.res = check_overlap_to_merge(raw.modules)
  
  
  write.table(mom.res$MoM.table,file = paste(res.dir,"/Overlap_Results.Module_Membership.txt",sep = ""),sep = "\t",
              row.names = FALSE,col.names= TRUE,quote = FALSE)
  
}

############################### Now, apply the following criteria: FET p > 0.05 & pairs are not preserved in any other cancers
require(reshape2)
conserve.mat = acast(data = all.res,formula= Module ~ ctrl.id,value.var = "is.conserved",fun.aggregate = function(x) sum(x,na.rm = TRUE))
no.conserve.mat = acast(data = all.res,formula= Module ~ ctrl.id,value.var = "is.not.conserved",fun.aggregate = function(x) sum(x,na.rm = TRUE))

unique.mods = rownames(no.conserve.mat)[rowSums(no.conserve.mat > 0,na.rm = TRUE) == (ncol(no.conserve.mat)-1)]
MEGENA::output.geneSet.file(list("Cancer.Specific.Modules" = unique.mods),gsub("\\.txt","\\.Cancer_Specific_Modules\\.txt",outfile))

############################### Now, run MoM analysis
### Params

out.dir = res.dir

### load modules
modules.o = MEGENA::read.geneSet(modfile)
modules = lapply(modules.o,function(x) gsub("\\|(.*)$","",x))

##### construct MoM-wise network
net.el = read.delim(file = netfile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
cancers = unique(net.el$group)

if (TRUE)
{
  mom.netres = identify_MoM_Network(net.el,cls.list,modules,out.dir = out.dir,core.cutoff = 0.95,plot.network = TRUE)
  sapply(mom.netres$mom.cores,length)
  sapply(mapply(FUN = function(x,y) intersect(x,y),x = mom.netres$mom.cores,y = mom.netres$mom.overlap,SIMPLIFY = FALSE),length)
  mom.netres$mom.final = lapply(mom.netres$node.stats,function(x) subset(x,overlap.count.normalized > 0.5 & is.core)$gene)
  sapply(mom.netres$mom.final,length)
 
  #mean(core.overlap)
  
  # merge those that show intersection larger 
  
  if (!all(sapply(mom.netres$plotlist,is.null)))
  {
    pdf(file = paste(out.dir,"/Merged_MoM_Network.pdf",sep = ""),width = 17,height = 15)
    for (i in 1:length(mom.netres$plotlist)) print(mom.netres$plotlist[[i]])
    dev.off()
  }
  
  MEGENA::output.geneSet.file(mom.netres$mom.modules,paste(out.dir,"/MoM.modules.txt",sep = ""))
  MEGENA::output.geneSet.file(mom.netres$mom.cores,paste(out.dir,"/MoM.cores.txt",sep = ""))
  MEGENA::output.geneSet.file(cls.list,paste(out.dir,"/MoM.module_list.txt",sep = ""))
  write.table(data.frame(MoM = rownames(mom.netres$link.count),as.data.frame(mom.netres$link.count),total.link = rowSums(mom.netres$link.count,na.rm = TRUE)),
              file = paste(out.dir,"/MoM.link_counts.txt",sep = ""),
              sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  saveRDS(mom.netres,file = paste(out.dir,"/MoM_Network_Results.RDS",sep = ""))
}else{
  mom.netres = readRDS(paste(out.dir,"/MoM_Network_Results.RDS",sep = ""))
}
###### check MoM enrichment with protein signatures
if (TRUE)
{
  ## run FET
  if (FALSE)
  {
    source("scripts/R_functions/enrichment_functions.R")
    ### set up functions and parameters for signature extraction
    pval.cutoff = 0.05
    fc.cutoff = 2
    min.size = 5
    max.size = 3000
    
    bg = union(net.el[[1]],net.el[[2]])
    
    # functions
    if (TRUE)
    {
      get_limma_sig = function(x,pv,fc)
      {
        df = subset(x,adj.p.value < pv & abs(logFC) > log2(fc))
        #df = subset(x,p.value < pv & abs(logFC) > log2(fc))
        if (nrow(df) > 0)
        {
          df$SIG.ID = rep(NA,nrow(df))
          df$SIG.ID[df$logFC > 0] = "UP"
          df$SIG.ID[df$logFC < 0] = "DN"
        }
        return(df)
      }
      
      get_JT_sig = function(x,pv)
      {
        df = subset(x,adj.p.value < pv)
        if (nrow(df) > 0)
        {
          df$SIG.ID = rep(NA,nrow(df))
          df$SIG.ID[df$p.value.inc < df$p.value.dec] = "INC"
          df$SIG.ID[df$p.value.inc > df$p.value.dec] = "DEC"
        }
        return(df)
      }
      
      get_COR_sig = function(x,pv)
      {
        df = subset(x,adj.p.value < pv)
        if (nrow(df) > 0)
        {
          df$SIG.ID = rep(NA,nrow(df))
          df$SIG.ID[df$statistic.rho < 0] = "NEG"
          df$SIG.ID[df$statistic.rho > 0] = "POS"
        }
        return(df)
      }
      
      extract_sig_results = function(stat.results,sig.categs,sig.fx.lst,pval.cutoff = 0.05,fc.cutoff = 1.2)
      {
        # get signatures per type
        sig.results = vector("list",length(stat.results))
        names(sig.results) = names(stat.results)
        for (i in 1:length(stat.results))
        {
          sig.lst = vector("list",length(stat.results[[i]]))
          names(sig.lst) = names(stat.results[[i]])
          for (j in 1:length(stat.results[[i]]))
          {
            sig.type = names(sig.categs)[sapply(sig.categs,function(x,y) any(y %in% x),y = names(stat.results[[i]])[j])]
            sig.func = sig.fx.lst[sig.type][[1]]
            if (sig.type == "limma")
            {
              sig.df = sig.func(x = stat.results[[i]][[j]],pv = pval.cutoff,fc = fc.cutoff)
            }else{
              sig.df = sig.func(x = stat.results[[i]][[j]],pv = pval.cutoff)
            }
            
            sig.lst[[j]] = sig.df
            rm(sig.type,sig.func,sig.df)
          }
          sig.results[[i]] = sig.lst
          rm(sig.lst)
        }
        return(sig.results)
      }
      
      make_signature_list = function(sig.results)
      {
        sig.list = vector("list",length(sig.results))
        names(sig.list) = names(sig.results)
        for (i in 1:length(sig.results))
        {
          sig.lst = list()
          for (j in 1:length(sig.results[[i]]))
          {
            if (nrow(sig.results[[i]][[j]]) > 0)
            {
              sig = split(sig.results[[i]][[j]][[1]],factor(sig.results[[i]][[j]]$SIG.ID))
              names(sig) = paste(names(sig.results[[i]])[j],names(sig),sep = "_")
              sig.lst = c(sig.lst,sig)
              rm(sig)
            }
          }
          sig.list[[i]] = sig.lst;
          rm(sig.lst)
        }
        return(sig.list)
      }
      
      sig.categs = list(limma = c("node.status","disease.stage","DEP","Immune.infiltration","overall.survival","MSI","EBV","cirrhosis","lymph.node.metastasis"),
                        JT = c("overall.stage","T.stage","N.stage","M.stage","grade","ISUP"),
                        COR = c("PSA"))
      
      sig.fx.lst = list(limma = get_limma_sig,
                        JT = get_JT_sig,
                        COR = get_COR_sig)
    }
    
    # load stat 
    stat.results = readRDS("Data/PRO_Clinical_Signatures.RDS")
    
    # Now, update with DEP results
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
    
    # attach to relevant cancer data
    for (i in 1:length(depres))
    {
      ii = which(names(stat.results) == names(depres)[i])
      stat.results[[ii]] = c(stat.results[[ii]],list(DEP = depres[[i]]))
    }
    
    
    ### extract signatures
    if (TRUE)
    {
      # extract signatures with pvalue cutoff
      sig.results = extract_sig_results(stat.results,sig.categs,sig.fx.lst,pval.cutoff = pval.cutoff,fc.cutoff = fc.cutoff)
      sig.list = make_signature_list(sig.results)
      sig.list = lapply(sig.list,function(x,mn) x[sapply(x,length) >= mn],mn = min.size)
      sapply(sig.list,names)
      sapply(sig.list,function(x) sapply(x,length))
      
      # add pan-cancer mutation signatures to each
      mut.tbl = read.delim(file = "Data/MutationalDrivers/mutation_driver_list.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
      mut.tbl$Cancer.type[mut.tbl$Cancer.type == "COADREAD"] = "CRC"
      mut.tbl$Cancer.type[mut.tbl$Cancer.type == "KIRC"] = "CCRCC"
      mut.tbl$Cancer.type[mut.tbl$Cancer.type == "LIHC"] = "HCC"
      
      mut.sigs = split(mut.tbl[[1]],factor(mut.tbl$Cancer.type))
      names(mut.sigs)[names(mut.sigs) %in% names(sig.list)]
      
      for (i in 1:length(sig.list))
      {
        if (any(names(mut.sigs) == names(sig.list)[i])) sig.list[[i]] = c(sig.list[[i]],list(mutation.driver_CAN = mut.sigs[names(sig.list)[i]][[1]]))
        #sig.list[[i]] = c(sig.list[[i]])
      }
      
      # collapse for enrichment testing
      sig.list.collapse = list()
      for (n in 1:length(sig.list)) 
      {
        sig = sig.list[[n]]
        if (length(sig) > 0)
        {
          names(sig) = paste(names(sig.list)[n],names(sig),sep = "_")
          sig.list.collapse = c(sig.list.collapse,sig)
          
        }
        rm(sig)
      }
      sig.list.collapse = sig.list.collapse[sapply(sig.list.collapse,length) >= 5]
      sig.list.collapse = c(sig.list.collapse,list(PANCAN_mutation.driver_ALL = mut.sigs$PANCAN))
    }
    sapply(sig.list.collapse,length)
    
    core.sig.FET = perform.AllPairs.FET(geneSets1 = mom.netres$mom.cores,geneSets2 = sig.list.collapse,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    core.sig.FET$corrected.FET.pvalue = p.adjust(core.sig.FET$FET_pvalue,"BH")
    
    sig.a = do.call('rbind',strsplit(as.character(core.sig.FET$set2_Name),"_"))
    core.sig.FET$cancer = sig.a[,1]
    core.sig.FET$signature.id = paste(sig.a[,2],sig.a[,3],sep = "_")
    core.sig.FET = core.sig.FET[order(core.sig.FET$FET_pvalue),]
    
    write.table(core.sig.FET,file = paste(out.dir,"/MoM.cores.Clinical_FET.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  }else{
    core.sig.FET = read.delim(file = paste(out.dir,"/MoM.cores.Clinical_FET.txt",sep = ""),sep = "\t",header =TRUE,stringsAsFactors = FALSE)
  }
  
  ## make summary heatmap
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  uniq.sigs = unique(core.sig.FET$signature.id)
  
  # get signatures hit by the cores
  sig.terms = subset(core.sig.FET,corrected.FET.pvalue < 0.05)$set2_Name
  plot.data = subset(core.sig.FET,set2_Name %in% sig.terms)
  
  # order by hclust
  require(reshape2)
  pmat = acast(data = plot.data,formula = set1_Name ~ set2_Name,value.var = "corrected.FET.pvalue",fun.aggregate = function(x) min(x,na.rm = TRUE))
  pmat[is.infinite(pmat)] = 1;
  hobj = hclust(dist(-log10(pmat)),"complete")
  mom.lab = hobj$labels[hobj$order]
  # re-label signature names
  labmap = c("DEP_DN" = "tumorDEP:DN","DEP_UP" = "tumorDEP:UP","EBV_UP"= "EBV:UP","grade_INC" = "Tumor\nGrade:\nINC","overall.stage_INC" = "Cancer\nstage:\nINC")
  plot.data$signature.id = labmap[plot.data$signature.id]
  
  # make heatmap
  heatobj = ggplot(data = plot.data,aes(x = cancer,y = set1_Name,fill = -log10(corrected.FET.pvalue))) + geom_tile() + 
    scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = -log10(0.05)) + 
    scale_y_discrete(limits = mom.lab) + 
    geom_text(data = subset(plot.data,corrected.FET.pvalue < 0.05),aes(x = cancer,y = set1_Name,label = signif(enrichment.foldchange,2))) + 
    facet_grid(. ~ signature.id,scale = "free_x",space = "free_x") + theme_bw() + 
    guides(fill = guide_colorbar(title = "-log10(FDR)",title.position = "top")) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 16),axis.text.y = element_text(size = 16),
          axis.title.y = element_blank(),axis.title.x = element_text(size = 19),
          legend.position = "bottom",legend.direction = "horizontal",legend.text = element_text(size = 13),legend.title = element_text(size = 16),
          strip.text = element_text(size = 16)
          )
  png(file = paste(out.dir,"/MoM.cores.Clinical_FET.png",sep = ""),res = 300,width = 2300,height = 1900)
  print(heatobj)
  dev.off()
}

###### check overlaps with real modules
if (TRUE)
{
  source("scripts/R_functions/enrichment_functions.R")
  bg = Reduce("union",modules)
  
  mom.FET = perform.AllPairs.FET(geneSets1 = mom.netres$mom.cores,geneSets2 = modules,background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
  mom.FET$corrected.FET.pvalue = p.adjust(mom.FET$FET_pvalue,"BH")
  mom.FET$cancer = gsub("^(.*)\\|","",as.character(mom.FET$set2_Name))
  mom.FET = mom.FET[order(mom.FET$FET_pvalue),]
  mom.FET$Proportion = mom.FET$actual.overlap/mom.FET$set1_size
  
  # make the call 
  mom.FET$call = rep(NA,nrow(mom.FET))
  mom.FET$call[which(mom.FET$corrected.FET.pvalue < 0.05 & mom.FET$Proportion > 0.5)] = "Overlap"
  mom.FET$call[which(mom.FET$Proportion > 0.95 & mom.FET$corrected.FET.pvalue < 0.05 & mom.FET$enrichment.foldchange >= 2)] = "Preserved"
  mom.FET$call[which(mom.FET$Proportion > 0.95 & mom.FET$corrected.FET.pvalue < 0.05 & (mom.FET$enrichment.foldchange < 2 & mom.FET$enrichment.foldchange >= 1))] = "Inclusive"
  mom.FET$call[which(mom.FET$corrected.FET.pvalue > 0.05)] = "Not Preserved"
  table(mom.FET$call)
  
  require(reshape2)
  require(ggplot2)
  #sig.mom.FET = subset(mom.FET,corrected.FET.pvalue < 0.05)
  sig.mom.FET = mom.FET
  sig.mom.ii = sapply(split(1:nrow(sig.mom.FET),factor(paste(sig.mom.FET$set1_Name,sig.mom.FET$cancer,sep = "_"))),function(ii,x) ii[which.min(x$FET_pvalue[ii])],x = sig.mom.FET)
  
  plot.data = sig.mom.FET[sig.mom.ii,]
  core.pobj = ggplot(data = plot.data,aes(x = cancer,y = set1_Name,fill = -log10(corrected.FET.pvalue))) + geom_tile() + 
    scale_x_discrete(limits = c("BRCA","CCRCC","CRC","HCC","LUAD","STAD","UCEC")) + 
    scale_y_discrete(limits = rev(names(mom.netres$mom.cores))) + 
    scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = -log10(0.05)) + 
    geom_point(data = subset(plot.data,!is.na(call)),aes(x = cancer,y = set1_Name,colour = call,size = Proportion)) + 
    scale_colour_manual(values = c("Preserved" = "chartreuse3","Not Preserved" = "black","Overlap" = "darkorange"),
                        labels = c("Preserved" = "Preserved\n(Prop > 95%,\nFDR < 0.05, EFC > 2)",
                                   "Not Preserved" = "Not Preserved\n(FDR > 0.05)",
                                   "Overlap" = "Overlap\n(FDR < 0.05,\nProp > 50% & < 95%)")) + 
    guides(colour = guide_legend(title = "Status",ncol = 1,override.aes = list(size = 6)),
           fill = guide_colourbar(title = "-log10(FDR)",title.position = "top"),
           size = guide_legend(ncol = 2,title.position = "top")) + 
    labs(y = "") + 
    theme_bw() + 
    theme(legend.position = "right",axis.text.y = element_text(size = 16),axis.text.x = element_text(angle = 45,vjust = 1,hjust =1,size = 16),
          axis.title = element_text(size = 20),legend.title = element_text(size = 18),legend.text = element_text(size = 15))
  core.pobj
  call.mat = acast(data = mom.FET,formula = set1_Name ~ cancer,value.var = "call",fun.aggregate = function(x) any(x == "Preserved",na.rm = TRUE))
  
  
  write.table(mom.FET,file = paste(out.dir,"/mom_core_preservation.FET.txt",sep = ""),sep = "\t",row.names= FALSE,col.names = TRUE,quote = FALSE)
  
}

############# evaluate inter-MoM interactions in each cancer
if (TRUE)
{
  library(GSVA)
  data.files = list.files(path = "Data/QCed_Proteome",pattern = "BatchAgeGender_Corrected.txt$",full.names= TRUE)
  names(data.files) = gsub("^Data/QCed_Proteome/|\\.BatchAgeGenderRace_Corrected.txt$","",data.files)
  min.size = 3;
  
  mom.cor = data.frame()
  for (i in 1:length(data.files))
  {
    data.mat = as.matrix(read.delim(file = data.files[i],sep = "\t",header = TRUE,stringsAsFactors = FALSE))
    
    # skip the union sets
    if (FALSE)
    {
      mom.rho = check_mom_corr(data.mat,mods.matched = mom.netres$mom.modules)
      mom.rho$cancer = rep(names(data.files)[i],nrow(mom.rho))
      mom.rho$type = rep("MoM",nrow(mom.rho))
      mom.cor = rbind.data.frame(mom.cor,mom.rho)
    }
    
    core.rho = check_mom_corr(data.mat,mods.matched = mom.netres$mom.cores)
    core.rho$cancer = rep(names(data.files)[i],nrow(core.rho))
    core.rho$type = rep("core",nrow(core.rho))
    
    mom.cor = rbind.data.frame(mom.cor,core.rho)
    rm(rhores)
  }
  
  write.table(mom.cor,file = paste(out.dir,"/inter_MoM_correlation.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  
  ### make heatmap summary
  input.cor = subset(mom.cor,type == "core")
  # decide order
  gen_mom_heatmap = function(input.cor)
  {
    require(reshape2)
    rmat = acast(data = input.cor,formula = mi ~ mj,value.var = "rho",fun.aggregate = function(x) median(x,na.rm = TRUE))
    hobj = hclust(dist(rmat),"complete")
    mtick = hobj$labels[hobj$order]
    
    # get tile plot
    hplot = ggplot(data = input.cor,aes(x = mi,y = mj,fill = -log10(adj.p.value)*sign(rho))) + 
      geom_tile() + facet_wrap( ~ cancer ) + scale_x_discrete(limits = mtick) + scale_y_discrete(limits = mtick) + 
      scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0,na.value = "White") + theme_bw() + 
      guides(fill = guide_colourbar(title = "-log10(FDR)xsign(rho)")) + 
      labs(x = "MoM cores",y = "") +
      theme(legend.position = "bottom",legend.direction = "horizontal",axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
    return(hplot)
  }
  
  # per cancer heatmap
  library(ComplexHeatmap)
  library(circlize)
  heat.col = colorRamp2(breaks = c(320,20,10,2,0,-2,-10,-20,-320),colors = c("red","red","deeppink","pink","white","cadetblue1","cyan","blue","blue"))
  
  mom.cor.res = mom.cor
  mom.cor.res$cancer = gsub("^(.*)/|\\.(.*)$","",mom.cor.res$cancer)
  mom.cor.res = subset(mom.cor.res,type == "core" & cancer %in% c("BRCA","CCRCC","CRC","LUAD","HCC","STAD","UCEC"))
  
  cancers = c("BRCA","CCRCC","CRC","LUAD","HCC","STAD","UCEC")
  ha.list = vector("list",length(cancers));names(ha.list) = cancers
  for (i in 1:length(cancers))
  {
    mom.can = subset(mom.cor.res,cancer == cancers[i])
    mom.can$adj.p.value[mom.can$adj.p.value < 1E-320] = 1E-320
    mom.can$composite = sign(mom.can$rho) * -log10(mom.can$adj.p.value)
    
    # get heatmap values
    rmat = acast(data = mom.can,formula = mi ~ mj,value.var = "composite",fun.aggregate = function(x) median(x,na.rm = TRUE))
    rmat = rmat[match( names(mom.netres$mom.cores),rownames(rmat)),match( names(mom.netres$mom.cores),colnames(rmat))]
    rmat[is.na(rmat)] = 0
    rownames(rmat) = colnames(rmat) = names(mom.netres$mom.cores)
    rmat = rmat + t(rmat);
    
    # order by dendrogram
    hobj = hclust(as.dist(sqrt(2*(1-cor(rmat,method = "spearman")))),"complete")
    ha.list[[i]] = Heatmap(rmat,col = heat.col,cluster_rows = hobj,cluster_columns = hobj,show_heatmap_legend = FALSE,column_title = cancers[i])
  }
  
  # get legend 
  lgd.col = colorRamp2(breaks = c(20,10,2,0,-2,-10,-20),colors = c("red","deeppink","pink","white","cadetblue1","cyan","blue"))
  heat.lgd = Legend(col_fun = lgd.col,title = "sign(rho) * -log10(FDR)",direction = "horizontal",
                    legend_width = unit(4,"cm"),labels_gp = gpar(size = 8))
  #draw(heat.lgd)
  
  # output pdf
  pdf(file = paste(out.dir,"/inter_MoM_core_correlation_per_cancer.pdf",sep = ""),width = 5,height = 5.5)
  for (i in 1:length(ha.list)) draw(ha.list[[i]])
  dev.off()
  
  ###### get co-occurrence
  mom.cor.res$pair.id = paste(mom.cor.res$mi,mom.cor.res$mj,sep = "_")
  vec = rep(NA,nrow(mom.cor.res))
  vec[mom.cor.res$rho > 0 & mom.cor.res$adj.p.value < 0.05] = "POS"
  vec[mom.cor.res$rho < 0 & mom.cor.res$adj.p.value < 0.05] = "NEG"
  mom.cor.res$cor.id = vec
  
  co.tbl = table(mom.cor.res$pair.id,mom.cor.res$cor.id)
  co.tbl = rbind.data.frame(data.frame(mi = gsub("_(.*)$","",rownames(co.tbl)),mj = gsub("^(.*)_","",rownames(co.tbl)),count = co.tbl[,"POS"]),
                            data.frame(mi = gsub("^(.*)_","",rownames(co.tbl)),mj = gsub("_(.*)$","",rownames(co.tbl)),count = -co.tbl[,"NEG"]))
  
  co.mat = acast(data = co.tbl,formula = mi ~ mj,value.var = "count")
  co.mat[is.na(co.mat)] = 0
  co.mat = co.mat[names(mom.netres$mom.cores),names(mom.netres$mom.cores)]
  hobj = hclust(dist(co.mat),"complete")
  co.mat = co.mat[hobj$labels[hobj$order],hobj$labels[hobj$order]]
  
  for (i in 1:(nrow(co.mat)-1))
  {
    for (j in (i+1):nrow(co.mat))
    {
      if (co.mat[i,j] < 0) 
      {
        pos.val = co.mat[j,i];neg.val = co.mat[i,j]
        co.mat[i,j] = pos.val;co.mat[j,i] = neg.val
      }
      if (co.mat[j,i] > 0)
      {
        pos.val = co.mat[j,i]
        neg.val = co.mat[i,j]
        co.mat[i,j] = pos.val
        co.mat[j,i] = neg.val
      }
    }
  }
  
  co.tbl.align = melt(co.mat)
  colnames(co.tbl.align) = c("mi","mj","count")
  
  # get co-occurence matrix
  map.obj = ggplot(data = co.tbl.align,aes(x = mi,y = mj,fill = count)) + geom_tile() + 
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) + theme_minimal() + 
    scale_x_discrete(limits = rownames(co.mat)) + scale_y_discrete(limits = colnames(co.mat)) + 
    geom_text(data = subset(co.tbl.align,count != 0),aes(x = mi,y = mj,label = abs(count))) + 
    #scale_x_discrete(limits= hobj$labels[hobj$order]) + scale_y_discrete(limits= hobj$labels[hobj$order]) + 
    guides(fill = guide_legend(title = "co-occurence")) + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 14),axis.title = element_blank(),
          axis.text.y = element_text(size = 14),legend.title = element_text(size = 16),legend.text = element_text(size = 13),
          legend.position = "bottom",legend.direction = "horizontal")
  map.obj
  
  # get MoM core network
  co.tbl = table(mom.cor.res$pair.id,mom.cor.res$cor.id)
  mx.tbl = co.tbl[co.tbl[,2] >= 4 | co.tbl[,1] >= 4,]
  mx.df = melt(mx.tbl);colnames(mx.df) = c("pair.id","cor.id","count")
  mx.df = subset(mx.df,count > 0)
  mx.df$mi = gsub("_(.*)$","",mx.df$pair.id)
  mx.df$mj = gsub("^(.*)_","",mx.df$pair.id)
  mx.df = mx.df[,c(4,5,3,2)]
  
  # add functional annotation 
  mom.annot = c("HEME metabolism","Proteasome Accessory","NADH Dehydrogenase","Translation Termination","INF-A,B signaling","Mitochodrial Translation",
                "Integrator Complex","Proteasome","NADH Dehydrogenase","Oxidoreductase Activity","Transcription Export","Telomerase RNA Localization","ER Lumen",
                "Immunoglobulin","Mitochondrial Translation","N Acetyltransferase Activity","Golgi Vesicle Transport","Defense Response","Golgi Transport","Blood Microparticle","Ribosome")
  names(mom.annot) = paste("MoM",1:21,sep = "")
  
  library(ggraph)
  library(ggrepel)
  graph <- graph_from_data_frame(mx.df);k = igraph::degree(graph);V(graph)$degree = k
  gobj = ggraph(graph,layout = "kk") + scale_edge_colour_manual(values = c("NEG" = "blue","POS" = "red")) + 
    geom_edge_link(aes(colour = factor(cor.id))) + 
    geom_node_point(aes(size = degree),colour = alpha("black",0.6)) + theme_minimal() + 
    guides(edge_colour = guide_legend(title = "Correlation"),size = FALSE) + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_blank(),axis.title = element_blank(),
          legend.position = "bottom",legend.title = element_text(size = 15),legend.text = element_text(size = 13))
  gobj$data$label = paste(gobj$data$name,"\n",mom.annot[gobj$data$name],sep = "")
  gobj = gobj + geom_text_repel(data = gobj$data,aes(x = x,y = y,label = label),size = 5)
  gobj
}

######################## MSigDB enrichments
## Now, check MoM enrichments for MSigDB
source("scripts/R_functions/enrichment_functions.R")
if (TRUE)
{
  # get background
  bg = Reduce("union",modules)
  print(length(bg))
  
  ## get GMT files
  gmt.folder = "Data/MSigDB"
  
  # get respective gmt files
  gmt.files <- list.files(path = gmt.folder,full.names = TRUE,pattern = "\\.gmt$")
  #gmt.files <- gmt.files[-grep("c6",gmt.files)]
  names(gmt.files) <- gsub("^(.*)/|\\.gmt$","",gmt.files)
  gmt.sets <- lapply(gmt.files,function(x) lapply(MEGENA::read.geneSet(x),function(y) y[-1]))
  
  #### run MSigDB enrichments: MoMs
  res.df <- data.frame()
  for (i in 1:length(gmt.sets))
  {
    df <- perform.AllPairs.FET(geneSets1 = mom.netres$mom.modules,geneSets2 = gmt.sets[[i]],background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    df$database = rep(names(gmt.sets)[i],nrow(df))
    df$corrected.FET.pvalue <- p.adjust(df$FET_pvalue,"BH")
    res.df <- rbind.data.frame(res.df,df);rm(df)
  }
  res.df$corrected.FET.pvalue = p.adjust(res.df$FET_pvalue,"BH")
  res.df = res.df[order(res.df$FET_pvalue),]
  
  # write results
  write.table(res.df,file = paste(out.dir,"/MoM.MSigDB_FET.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,
              quote = FALSE)
  
  # produce summary plot
  require(reshape2)
  top.n = 10
  sdf = subset(res.df,corrected.FET.pvalue < 0.05 & enrichment.foldchange > 2)
  mom.ids = sort(as.character(unique(sdf$set1_Name)))
  mom.obj = vector("list",length(mom.ids))
  mom.topterm = rep(NA,length(mom.ids))
  for (i in 1:length(mom.ids))
  {
    mom.sdf = subset(sdf,set1_Name == mom.ids[i])
    mom.sdf = mom.sdf[order(mom.sdf$FET_pvalue),]
    mom.sdf = mom.sdf[1:min(c(top.n,nrow(mom.sdf))),]
    mom.topterm[i] = as.character(mom.sdf$set2_Name[1])
    mom.obj[[i]] = ggplot(data = mom.sdf,aes(x = set2_Name,y = -log10(corrected.FET.pvalue))) + geom_bar(stat = "identity") + 
      scale_x_discrete(limits = mom.sdf$set2_Name) + geom_hline(yintercept = -log10(0.05),colour = "red") + 
      coord_flip() + labs(x = "-log10(FDR FET P)",y = "Enriched Terms",title = mom.ids[i]) + theme_bw()
    rm(mom.sdf)
  }
  pdf(file =paste(out.dir,"/MoM.MSigDB_FET.pdf",sep = ""),width = 14,height = 8)
  for (i in 1:length(mom.obj)) print(mom.obj[[i]])
  dev.off()
  write.table(data.frame(MoM.id = mom.ids,top.term = mom.topterm),file = paste(out.dir,"/MoM.top_MSigDB.txt",sep = ""),
              sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  ####### MoM cores
  res.df <- data.frame()
  for (i in 1:length(gmt.sets))
  {
    df <- perform.AllPairs.FET(geneSets1 = mom.netres$mom.cores,geneSets2 = gmt.sets[[i]],background = bg,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    df$database = rep(names(gmt.sets)[i],nrow(df))
    df$corrected.FET.pvalue <- p.adjust(df$FET_pvalue,"BH")
    res.df <- rbind.data.frame(res.df,df);rm(df)
  }
  res.df$corrected.FET.pvalue = p.adjust(res.df$FET_pvalue,"BH")
  res.df = res.df[order(res.df$FET_pvalue),]
  
  # write results
  write.table(res.df,file = paste(out.dir,"/MoM_Cores.MSigDB_FET.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,
              quote = FALSE)
  
  # produce summary plot
  require(reshape2)
  top.n = 10
  sdf = subset(res.df,corrected.FET.pvalue < 0.05 & enrichment.foldchange > 2)
  mom.ids = sort(as.character(unique(sdf$set1_Name)))
  mom.obj = vector("list",length(mom.ids))
  mom.topterm = rep(NA,length(mom.ids))
  for (i in 1:length(mom.ids))
  {
    mom.sdf = subset(sdf,set1_Name == mom.ids[i])
    mom.sdf = mom.sdf[order(mom.sdf$FET_pvalue),]
    mom.sdf = mom.sdf[1:min(c(top.n,nrow(mom.sdf))),]
    mom.topterm[i] = as.character(mom.sdf$set2_Name[1])
    mom.obj[[i]] = ggplot(data = mom.sdf,aes(x = set2_Name,y = -log10(corrected.FET.pvalue))) + geom_bar(stat = "identity") + 
      scale_x_discrete(limits = mom.sdf$set2_Name) + geom_hline(yintercept = -log10(0.05),colour = "red") + 
      coord_flip() + labs(x = "-log10(FDR FET P)",y = "Enriched Terms",title = mom.ids[i]) + theme_bw()
    rm(mom.sdf)
  }
  pdf(file =paste(out.dir,"/MoM_Cores.MSigDB_FET.pdf",sep = ""),width = 14,height = 8)
  for (i in 1:length(mom.obj)) print(mom.obj[[i]])
  dev.off()
  write.table(data.frame(MoM.id = mom.ids,top.term = mom.topterm),file = paste(out.dir,"/MoM_Cores.top_MSigDB.txt",sep = ""),
              sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  
}

# plot top terms in the cores
plot_best_MSigDB <- function(data.res)
{
  require(ggplot2)
  # show barplots of top ranked terms per signatures
  idx = sapply(split(1:nrow(data.res),factor(data.res$set1_Name)),function(x,y) x[which.min(y$FET_pvalue[x])],y = data.res)
  plot.data = data.res[idx,]
  bar.max = ggplot() + 
    geom_bar(data = plot.data,aes(x = set1_Name,y = -log10(corrected.FET.pvalue)),stat = "identity",fill = "cornflowerblue") + 
    geom_label(data = subset(plot.data,corrected.FET.pvalue < 0.05),
               y = max(-log10(plot.data$corrected.FET.pvalue))/2,
               aes(x = set1_Name,label = set2_Name),alpha = 0.5) + 
    geom_text(data = subset(plot.data,corrected.FET.pvalue < 0.05),
               y = -log10(0.05),
               aes(x = set1_Name,label = signif(enrichment.foldchange,2)),colour = "black",hjust = -0.2) + 
    theme_bw() + 
    labs(x = "MoM",y = "-log10(FET FDR)") +
    geom_hline(yintercept = -log10(0.05),colour = "red") + 
    coord_flip() + theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),strip.text = element_text(size = 15))
  bar.max
}

data.res = read.delim(file = paste(out.dir,"/MoM_Cores.MSigDB_FET.txt",sep = ""),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
bar.max = plot_best_MSigDB(data.res)
bar.max = bar.max + scale_x_discrete(limits = rev(names(mom.netres$mom.cores)))
bar.max

# output list of core proteins
tbl = data.frame(MoM = names(mom.netres$mom.cores),Core.proteins = sapply(mom.netres$mom.cores,function(x) paste(x,collapse = ",")))
write.table(tbl,file = paste(out.dir,"/MoM_Cores.table.txt",sep = ""),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
###### Now, see how MoMs interact within each network 
if (TRUE)
{
  # MoM interactions
  cancers = unique(net.el$group)
  mom.modules = mom.netres$mom.modules
  mom.inter = summarize_mom_interaction(gene.lst = mom.modules,net.el)
  mom.inter$weight = mom.inter$ne/mom.inter$nij 
  rownames(mom.inter) = NULL
  mom.all = gen_mom_network(df = subset(mom.inter, weight > 0.1),sigs = mom.modules)
  all.mom = c(lapply(cancers,function(x,y,z) gen_mom_network(df = subset(y,cancer.id == x & weight > 0.1),sigs = z),y = mom.inter,z = mom.modules),
              list(all = gen_mom_network(df = subset(mom.inter, weight > 0.1),sigs = mom.modules)))
  all.mom = all.mom[c(which(sapply(cancers,function(x,y) nrow(subset(y,cancer.id == x & weight > 0.1)),y = mom.inter) > 0),length(all.mom))]
  pdf(file = paste(out.dir,"/MoM_Interaction.pdf",sep = ""),width = 8,height = 6)
  for (i in 1:length(all.mom)) print(all.mom[[i]])
  dev.off()
  
  # MoM core interactions
  mom.cores = mom.netres$mom.cores
  mom.core.inter = summarize_mom_interaction(gene.lst = mom.cores,net.el)
  mom.core.inter$weight = mom.core.inter$ne/mom.core.inter$nij 
  rownames(mom.core.inter) = NULL
  
  all.cores = c(lapply(cancers,function(x,y,z) gen_mom_network(df = subset(y,cancer.id == x & weight > 0.1),sigs = z),y = mom.core.inter,z = mom.cores),
                list(all = gen_mom_network(df = subset(mom.core.inter, weight > 0.1),sigs = mom.cores)))
  all.cores = all.cores[c(which(sapply(cancers,function(x,y) nrow(subset(y,cancer.id == x & weight > 0.1)),y = mom.core.inter) > 0),length(all.cores))]
  
  pdf(file = paste(out.dir,"/MoM_Cores_Interaction.pdf",sep = ""),width = 8,height = 6)
  for (i in 1:length(all.mom)) print(all.cores[[i]])
  dev.off()
}


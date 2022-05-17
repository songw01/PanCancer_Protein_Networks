
assign_overlap_call <- function(atbl,pval.cut = 0.05,efc.cut = 1.5,overlap.cut = 0.3)
{
  atbl$set1_overlap = atbl$actual.overlap/atbl$set1_size
  atbl$set2_overlap = atbl$actual.overlap/atbl$set2_size
  atbl$overall_overlap = atbl$actual.overlap/(atbl$set1_size + atbl$set2_size - atbl$actual.overlap)
  # get significant conservation 
  atbl$is.set1.conserved = atbl$corrected.FET.pvalue < pval.cut & atbl$enrichment.foldchange >= efc.cut & atbl$set1_overlap > overlap.cut
  atbl$is.set2.conserved = atbl$corrected.FET.pvalue < pval.cut & atbl$enrichment.foldchange >= efc.cut & atbl$set2_overlap > overlap.cut
  
  # categorize based on enrichments
  atbl$conservation.class = rep(NA,nrow(atbl))
  atbl$conservation.class[atbl$corrected.FET.pvalue < pval.cut] = "enriched"
  atbl$conservation.class[atbl$corrected.FET.pvalue > 0.2] = "not.enriched"
  atbl$conservation.class[atbl$corrected.FET.pvalue < pval.cut & atbl$enrichment.foldchange > efc.cut & (atbl$is.set1.conserved)] = "set1.conserved"
  atbl$conservation.class[atbl$corrected.FET.pvalue < pval.cut & atbl$enrichment.foldchange > efc.cut & (atbl$is.set2.conserved)] = "set2.conserved"
  atbl$conservation.class[atbl$corrected.FET.pvalue < pval.cut & atbl$enrichment.foldchange > efc.cut & (atbl$is.set1.conserved & atbl$is.set2.conserved)] = "both.conserved"
  atbl = atbl[order(atbl$FET_pvalue),]
  
  return(atbl)
}

### Mom functions
summarize_mom_interaction = function(gene.lst,net.el)
{
  inter.res = data.frame()
  for (i in 1:(length(gene.lst)-1))
  {
    for (j in (i+1):length(gene.lst))
    {
      el = subset(net.el,(row %in% gene.lst[[i]] & col %in% gene.lst[[j]]) | (col %in% gene.lst[[i]] & row %in% gene.lst[[j]]))
      res = do.call('rbind',lapply(split(el$weight,factor(el$group)),function(x) c(sum(x,na.rm = TRUE),length(x))))
      if (!is.null(res))
      {
        res = data.frame(cancer.id = rownames(res),strength = res[,1],ne = res[,2],
                         mi = rep(names(gene.lst)[i],nrow(res)),mj = rep(names(gene.lst)[j],nrow(res)),
                         ni = rep(length(gene.lst[[i]]),nrow(res)),nj = rep(length(gene.lst[[j]]),nrow(res)),
                         nij = rep(length(union(gene.lst[[i]],gene.lst[[j]])),nrow(res)))
        inter.res = rbind.data.frame(inter.res,res)
      }
      
      rm(res,el)
    }
  }
  return(inter.res)
}

gen_mom_network <- function(df,sigs)
{
  mom.core.g = graph.data.frame(df[,c(4,5,9,setdiff(1:ncol(df),c(4,5,9)))],directed = FALSE,
                                vertices = data.frame(id = names(sigs),module.size = sapply(sigs,length)))
  
  # get cores
  comp.res = components(mom.core.g)
  comp.id = which(comp.res$csize > 2)
  core.mom = do.call('c',lapply(comp.id,function(x,g,mem) {sg = subgraph(g,v = names(mem)[mem == x]);core.mom = coreness(sg);names(core.mom)[core.mom == max(core.mom)]},g = mom.core.g,mem = comp.res$membership))
  vertex_attr(mom.core.g)$is.core = vertex_attr(mom.core.g)$name %in% core.mom
  mom.core.pobj = ggraph(graph = mom.core.g,layout = "kk") + geom_edge_fan(aes(colour = cancer.id,edge_width = weight),strength = 1.5,alpha= 0.2) + 
    scale_edge_width_continuous(range = c(0.5,1.5)) + 
    geom_node_point(aes(size = module.size,colour = is.core),stroke = 0.2,alpha = 0.7)
  txt.dat = mom.core.pobj$data
  mom.core.pobj = mom.core.pobj + 
    geom_text_repel(data = txt.dat,
                    aes(x = x,y = y,label = name,colour = is.core),
                    fontface = "bold",
                    min.segment.length = 0.1,size = 5) + guides(colour = guide_legend(alpha = 1)) + 
    scale_colour_manual(values = c("TRUE" = "red","FALSE" = "grey"))
  return(mom.core.pobj)
}

identify_MoM <- function(conserve.table,modules,min.comp.size = 5)
{
  require(igraph)
  # conserve.table contains list of pairs of preserved modules. Contain set1_Name, set2_Name, cancer.i and cancer.j column
  conserve.table$set1_id = paste(conserve.table$set1_Name,conserve.table$cancer.i,sep = "|")
  conserve.table$set2_id = paste(conserve.table$set2_Name,conserve.table$cancer.j,sep = "|")
  
  el = data.frame(mi = conserve.table$set1_id,mj = conserve.table$set2_id,weight = conserve.table$overall_overlap)
  g = graph.data.frame(el,directed = FALSE)
  
  #gp = ggraph(g) + 
  #  geom_edge_link() + 
  #  geom_node_point()
  #### cluster per connected component
  ml = cluster_walktrap(graph = g,weights = E(g)$weight,steps = 2)
  table(membership(ml))
  
  # summarize into data.frame table
  ml.res = membership(ml)
  k = igraph::degree(g)
  df = data.frame(module.id = names(ml.res),cluster = as.vector(ml.res),cancer.id = gsub("^(.*)\\|","",names(ml.res)),k = k[names(ml.res)],
                  module.size = sapply(modules,length)[names(ml.res)])
  
  output = list(MoM.table = df,graph = g)
  return(output)
}

identify_MoM_Network <- function(net.el,cls.list,modules,out.dir,subset.cutoff = 0.8,core.cutoff = 0.95,overlap.cutoff = 0.95,plot.networks = FALSE)
{
  #################################
  cancers = unique(net.el$group)
  
  ### extract MoM-wise network
  # MoM bucket
  mom.modules = mom.cores = mom.overlap = mom.occur = mom.kcoreval = vector("list",length(cls.list))
  names(mom.modules) = names(mom.cores) = names(mom.overlap) = names(mom.occur)= names(mom.kcoreval) =  names(cls.list)
  
  # other params
  plst = vector("list",length(cls.list))
  names(plst) = names(cls.list)
  link.count = matrix(0,nrow = length(cls.list),ncol = length(cancers));
  rownames(link.count) = names(cls.list);colnames(link.count) = cancers
  require(reshape2)
  for (i in 1:length(cls.list))
  {
    cat(paste("processing:",names(cls.list)[i],"\n",sep = ""))
    # extract genes
    mod.match = modules[match(cls.list[[i]],names(modules))]
    cls.genes = Reduce("union",mod.match);
    mom.modules[[i]] = cls.genes
    co.occur = do.call('cbind',lapply(mod.match,function(x,y) y %in% x,y = cls.genes))
    rownames(co.occur) = cls.genes
    co.occur = co.occur[order(rowSums(co.occur),decreasing = TRUE),]
    head(co.occur)
    mom.occur[[i]] = co.occur
    
    ovr.vec = rowSums(co.occur)
    ovr.count = table(ovr.vec)
    top.ovr.val = as.integer(names(ovr.count)[which(cumsum(ovr.count)/length(ovr.vec) >= overlap.cutoff)])
    mom.overlap[[i]] = rownames(co.occur)[rowSums(co.occur) >= min(top.ovr.val)]
    #mom.overlap[[i]] = rownames(co.occur)[rowSums(co.occur) > ceiling(ncol(co.occur)/2)]
    # extract links
    cls.cid = unique(gsub("^(.*)\\|","",cls.list[[i]]))
    net.cls = subset(net.el,row %in% cls.genes & col %in% cls.genes & group %in% cls.cid)
    
    ### extract cores by centrality 
    # obtain cancer specific degree
    cids = unique(net.cls$group)
    
    # get degree
    kmat = do.call('cbind',lapply(cids,function(x,y,z) {g = graph.data.frame(subset(y,group == x),directed = FALSE);out = igraph::degree(g);out[match(z,names(out))]},y = net.cls,z = cls.genes))
    colnames(kmat) = paste("k",cids,sep = ".");rownames(kmat) = cls.genes
    kmat[is.na(kmat)] = 0;
    
    # get cores index score
    g = graph.data.frame(net.cls,directed = FALSE)
    k.overall = igraph::degree(g)
    core.score = igraph::coreness(g);
    core.count = table(core.score)
    top.core.val = as.integer(names(core.count)[which(cumsum(core.count)/length(core.score) >= core.cutoff)])
    mx.cores = names(core.score)[which(core.score %in% top.core.val)]
    mom.cores[[i]] = mx.cores
    
    # get closeness
    #close = igraph::closeness(graph = g,vids = V(g),mode = "all")
    
    # generate network figure
    ndf = data.frame(gene = rownames(kmat),as.data.frame(kmat),k.overall = k.overall[rownames(kmat)],
                     #closeness = close[rownames(kmat)],
                     coreness = core.score[rownames(kmat)],overlap.count = rowSums(co.occur)[rownames(kmat)],overlap.count.normalized = rowSums(co.occur)[rownames(kmat)]/length(cls.list[[i]]))
    ndf$is.core = ndf$gene %in% mx.cores
    ndf$show.label = rep(TRUE,nrow(ndf))
    mom.kcoreval[[i]] = ndf
    
    write.table(ndf,file = paste(out.dir,"/",names(cls.list)[i],".node.txt",sep = ""),sep = "\t",
                row.names = FALSE,col.names = TRUE,quote= FALSE)
    write.table(data.frame(gene = rownames(co.occur),as.data.frame(co.occur)),
                file = paste(out.dir,"/",names(cls.list)[i],".co_occurrence.txt",sep = ""),sep = "\t",
                row.names = FALSE,col.names = TRUE,quote= FALSE)
    write.table(net.cls,
                file = paste(out.dir,"/",names(cls.list)[i],".edges.txt",sep = ""),sep = "\t",
                row.names = FALSE,col.names = TRUE,quote= FALSE)
    
    ### 
    
    # add network plot
    if (plot.networks)
    {
      cls.g = graph.data.frame(net.cls,directed = FALSE,vertices = ndf)
      require(ggraph)
      require(ggrepel)
      pobj = ggraph(graph = cls.g,layout = "kk") + geom_edge_fan(aes(colour = group,edge_width = weight),strength = 1.5,alpha= 0.2) + 
        scale_edge_width_continuous(range = c(0.5,1.5)) + 
        geom_node_point(aes(size = k.overall),stroke = 0.2,alpha = 0.7) 
      txt.dat = pobj$data;
      pobj = pobj + geom_text_repel(data = subset(txt.dat,is.core),
                                    aes(x = x,y = y,label = name),
                                    fontface = "bold",
                                    min.segment.length = 0.1,size = 5) +
        labs(title = names(cls.list)[i]) + guides(colour = guide_legend(alpha = 1)) 
      plst[[i]] = pobj
    }
    
    ### summarize links
    ln = sapply(split(1:nrow(net.cls),factor(net.cls$group)),length)
    link.count[i,match(names(ln),colnames(link.count))] = ln
    
  }
  
  output = list(mom.modules = mom.modules,mom.cores = mom.cores,mom.overlap = mom.overlap,co.occurence = mom.occur,
                node.stats = mom.kcoreval,
                plotlist = plst,link.count = link.count)
  
  return(output)
}

check_mom_corr <- function(data.mat,mods.matched)
{
  # mods.matched = mom.netres$mom.modules
  mods.matched = lapply(mods.matched,function(x,y) intersect(x,y),y = rownames(data.mat))
  mods.matched = mods.matched[sapply(mods.matched,length) >= min.size]
  print(sapply(mods.matched,length))
  gsva.mat = gsva(expr = data.mat,gset.idx.list = mods.matched,method = "gsva")
  
  x <- y <- rho <- pval <- c();
  for (mi in 1:(nrow(gsva.mat)-1))
  {
    for (mj in (mi+1):nrow(gsva.mat))
    {
      out = cor.test(x = gsva.mat[mi,],y = gsva.mat[mj,],method = "spearman")
      x = c(x,rownames(gsva.mat)[mi])
      y = c(y,rownames(gsva.mat)[mj])
      rho = c(rho,out$estimate)
      pval = c(pval,out$p.value)
      rm(out)
    }
  }
  rhores = data.frame(mi = x,mj = y,rho = rho,p.value = pval,adj.p.value = p.adjust(pval,"BH"))
  output = list(gsva.matrix = gsva.mat,rho.results = rhores)
  return(rhores)
}

check_overlap_to_merge <- function(raw.modules,subset.cutoff = 0.9)
{
  ### screen out redundant moms
  raw.modules = lapply(cls.list,function(x,y) Reduce("union",y[x]),y = modules)
  
  # check if cores overlap
  core.overlap = matrix(0,length(raw.modules),length(raw.modules))
  for (i in 1:length(raw.modules))
  {
    core.overlap[i,-i] = sapply(raw.modules[-i],function(x,y) length(intersect(x,y)),y = raw.modules[[i]])
  }
  
  core.overlap.rate = core.overlap
  for (i in 1:length(raw.modules))
  {
    core.overlap.rate[i,] = core.overlap.rate[i,]/length(raw.modules[[i]])
  }
  
  overlap.ij = which(core.overlap.rate >= subset.cutoff,arr.ind = TRUE)
  list(overlap.matrix = core.overlap,overlap.rate.matrix = core.overlap.rate,overlap.pair = data.frame(from = names(raw.modules)[overlap.ij[,1]],to = names(raw.modules)[overlap.ij[,2]]))
}

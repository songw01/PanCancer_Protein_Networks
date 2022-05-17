# functions for signatures extraction
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
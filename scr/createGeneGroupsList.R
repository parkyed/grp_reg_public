createGeneGroupsList <- function(pathway.lst, gene.lst){
  
  #' Generate named list of gene groups
  #' 
  #' Create a named list of gene groups. Pathway names are the group names, and items are the member genes
  #' Filter the groups on a list of groups (e.g. all or enriched only) and list of genes (e.g. all or DE only)
  #' Create a orphan group for each gene not in any other group, named by that gene
  
  #' @param pathway.lst the list of pathways to include in the output
  #' @param gene.lst the list of genes to include in the output
  
  # filter pathway list on the genes in the gene list, and remove any empty groups
  pathway.lst.fil <- lapply(pathway.lst, function(x){x[which(x %in% gene.lst)]}) %>% purrr::compact() 
  
  # remove gene groups with single genes
  gene.grp.lst.pw <- pathway.lst.fil[which(!lapply(pathway.lst.fil, length) %>% unlist == 1)]
  
  # expand the groups list to include a group for each gene in the gene list not included in any other group - i.e. they have their own groups
  (genes.own.group <- gene.lst[!gene.lst %in% (gene.grp.lst.pw %>% unlist %>% unname %>% unique)])
  
  names(genes.own.group) <- genes.own.group
  
  (gene.grp.lst <- c(gene.grp.lst.pw, genes.own.group))
  
  # check the gene lists match
  if(length(dplyr::symdiff(gene.grp.lst %>% unlist %>% unname %>% unique, gene.lst)) != 0){stop()}
  
  return(list('gene.grp.lst' = gene.grp.lst, 'gene.grp.lst.pw' = gene.grp.lst.pw))
  
}
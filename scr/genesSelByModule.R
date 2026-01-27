genesSelByModule <- function(feat.names, feat.groups.idx, gl){
  
  #' Tally the number of genes selected for each module
  #' 
  #' Takes a list of selected genes, and an index of the selected groups in the full gene group list gl
  #' Filters the gene group list on the selected groups, and then filters on the selected genes
  #' Result is a list of selected gene groups, with the constituent genes as members
  #' Function then counts the number selected in each group
  #'
  #' @param feat.names names of all genes selected, regardless of group
  #' @param feat.groups.idx index of the groups selected by a given model
  #' @param gl named list of gene groups, names as group names, values as gene names

  # names of selected groups with the member genes
  (feat.groups <- gl[unlist(unname(feat.groups.idx))])
  
  # filter the groups on the selected genes only
  (feat.names.groups <- lapply(feat.groups, FUN = function(x){x[x %in% feat.names]}))
  
  # number of genes selected in each group
  num.list <- lapply(feat.names.groups, length)
  
  # create tibble
  num.list.tibble <- tibble::enframe(num.list, name = 'module')
  
  # convert value column from list to numeric
  num.list.tibble$value <- as.numeric(num.list.tibble$value)
  
  return(num.list.tibble)
  
}

# helper function

headerTrue <- function(df) {
  # function to replace column names in a dataframe with the first row
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

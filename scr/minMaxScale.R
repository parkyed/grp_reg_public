# min - max scaling of counts by gene
# note: the apply function transposes the output

minMaxScale <- function(x){(x-min(x))/(max(x)-min(x))}
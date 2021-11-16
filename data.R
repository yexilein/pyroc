
library(tidyverse)


mouse_go = function() {
    result = load_go("data/go_mouse")
    return(result)
}

load_go = function(file_prefix) {
    result = as.matrix(Matrix::readMM(paste0(file_prefix, ".txt")))
    rownames(result) = readLines(paste0(file_prefix, "_row_labels.txt"))
    colnames(result) = readLines(paste0(file_prefix, "_col_labels.txt"))
    return(result)
}
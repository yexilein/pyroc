
library(tidyverse)
library(Matrix)
library(ontologyIndex)
library(rjson)

this_directory__ = getwd()


main = function() {
    export_go_matrix(make_go_mouse(), "mouse", ".")
}

make_go_mouse = function(go_slim = FALSE) {
    result = t(make_go_matrix(file.path(this_directory__, "gene_association.mgi.gz"), go_slim))
    return(result)
}
                    
make_go_matrix = function(filename, go_slim = FALSE) {
    go_table = read_go_annotation(filename)
    go_terms = go_table$V5
    gene_symbols = go_table$V3
    
    ontology = get_ontology(file.path(this_directory__, "go-basic_181219.obo"), extract_tags = "everything")
    ancestor_matrix = get_ancestor_matrix_from_obo(ontology)
    alt_id_to_id = alt_id_to_id_from_obo(ontology)
    term_names = term_names_from_obo(ontology)
    
    is_alt_id = levels(go_terms) %in% names(alt_id_to_id)
    levels(go_terms)[is_alt_id] = alt_id_to_id[levels(go_terms)[is_alt_id]]
        
    result = sparseMatrix(as.numeric(gene_symbols), as.numeric(go_terms), dimnames = list(levels(gene_symbols), levels(go_terms)))
    is_known = colnames(result) %in% rownames(ancestor_matrix)
    warning(paste0(sum(!is_known), " unknown terms will be removed."))
    result = result[, is_known]
    result = result %*% ancestor_matrix[colnames(result),]
    if (go_slim) {
        is_go_slim = map_lgl(ontology$subset, ~ "goslim_generic" %in% .)
        result = result[, is_go_slim[colnames(result)]]
    }
    result = result[, colSums(result) > 0]
    colnames(result) = term_names[colnames(result)]
    return(result)
}
                   
read_go_annotation = function(filename) {
  return(read.csv(filename, header = FALSE, sep = "\t", comment.char = "!"))
}
                   
alt_id_to_id_from_obo = function(ontology) {
    result = ontology$alt_id %>%
        enframe("id", "alt_id") %>%
        unnest(alt_id) %>%
        dplyr::select(alt_id, id) %>%
        deframe
    return(result)
}
                   
term_names_from_obo = function(ontology) {
  result = paste(ontology$id, ontology$name, ontology$namespace, sep = "|")
  names(result) = ontology$id
  return(result)
}
                   
get_ancestor_matrix_from_obo = function(ontology) {
    ancestors = ontology$ancestors %>%
        enframe("term", "ancestor") %>%
        unnest(ancestor)
    terms = factor(ancestors$term)
    ancestors = factor(ancestors$ancestor, levels = levels(terms))
    result = sparseMatrix(as.numeric(terms), as.numeric(ancestors),
                          dimnames = list(levels(terms), levels(ancestors)))
}

export_go_matrix = function(go_matrix, file_prefix, subdir = "go") {
    writeMM(go_matrix, file.path(this_directory__, "data", subdir, paste0(file_prefix, ".txt")))
    write(rownames(go_matrix), file.path(this_directory__, "data", subdir, paste0(file_prefix, "_row_labels.txt")))
    write(colnames(go_matrix), file.path(this_directory__, "data", subdir, paste0(file_prefix, "_col_labels.txt")))
}

if (sys.nframe() == 0) {
    main()
}
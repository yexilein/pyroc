
library(tidyverse)
library(Matrix)


this_directory__ = getwd()
BIOGRID_NAME = "BIOGRID-ALL-4.4.197.tab2.txt"

export_mouse_ppi = function() {
    biogrid = load_biogrid_data()
    export_ppi(biogrid, mouse_tax_id(), "biogrid_mouse")
}

load_biogrid_data = function() {
    read.csv(file.path(this_directory__, BIOGRID_NAME), header = TRUE, comment.char = "", sep = "\t")
}

export_ppi = function(biogrid_data, tax_id, file_prefix) {
    result = build_ppi_from_biogrid(biogrid_data, tax_id)
    writeMM(result, paste0(file_prefix, "_ppi.txt"))
    write(colnames(result), paste0(file_prefix, "_genes.txt"))
}

build_ppi_from_biogrid = function(biogrid_data, tax_id) {
    ppi = biogrid_data %>%
        filter_physical() %>%
        filter_species(tax_id) %>%
        select(A = Official.Symbol.Interactor.A, B = Official.Symbol.Interactor.B) %>%
        droplevels()
    
    all_genes = unique(c(unique(ppi$A), unique(ppi$B)))
    ppi$A = factor(ppi$A, levels = all_genes)
    ppi$B = factor(ppi$B, levels = all_genes)
    result = sparseMatrix(as.numeric(ppi$A), as.numeric(ppi$B), dims = c(length(all_genes), length(all_genes)), dimnames = list(all_genes, all_genes))
    return(result)
}

filter_species = function(data, tax_id) {
    data %>%
        filter(Organism.Interactor.A == tax_id) %>%
        filter(as.character(Organism.Interactor.A) == as.character(Organism.Interactor.B))
}

mouse_tax_id = function() {
    return(10090)
}

filter_physical = function(data) {
    data %>%
        filter(Experimental.System.Type == "physical")
}

filter_known_techniques = function(data) {
    data %>%
        filter(Experimental.System %in% c("Two-hybrid", "Affinity Capture-MS"))
}

if (sys.nframe() == 0) {
    export_mouse_ppi()
}
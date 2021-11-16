
library(tidyverse)
library(SingleCellExperiment)
library(MetaMarkers)


main = function() {
    compute_tm_markers()
}

compute_tm_markers = function() {
    dataset = readRDS("~/data/tabula_muris/10x.rds")    
    cpm(dataset) = convert_to_cpm(assay(dataset), scale_factor = 1e4)
    labels = as.character(dataset$cell_ontology_class)
    labels[labels == ""] = "unknown"
    mouse_id = paste0(dataset$mouse.id, "_10x")
    markers_10x = lapply(set_names(unique(mouse_id)), function(id) {
        keep_cell = mouse_id == id
        compute_markers(cpm(dataset)[, keep_cell], labels[keep_cell], dataset$tissue[keep_cell])
    })
    markers_10x = bind_rows(markers_10x, .id = "mouse_id")

    dataset = readRDS("~/data/tabula_muris/smart_seq.rds")
    cpm(dataset) = convert_to_cpm(assay(dataset), scale_factor = 1e6)
    labels = as.character(dataset$cell_ontology_class)
    labels[labels == ""] = "unknown"
    mouse_id = paste0(dataset$mouse.id, "_ss")
    markers_ss = lapply(set_names(unique(mouse_id)), function(id) {
        keep_cell = mouse_id == id
        compute_markers(cpm(dataset)[, keep_cell], labels[keep_cell], dataset$tissue[keep_cell])
    })
    markers_ss = bind_rows(markers_ss, .id = "mouse_id")
    
    markers = bind_rows(markers_10x, markers_ss)
    export_markers(markers, "de_lists/markers_tm.csv")
}

if (!interactive()) {
    main()
}

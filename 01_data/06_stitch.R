
library(tidyverse)
library(Matrix)


this_directory__ = getwd()
STITCH_NAME = "10090.protein_chemical.links.detailed.v5.0.tsv.gz"
ALIAS_NAME = "10090.protein.aliases.v11.5.txt.gz"

export_mouse_dti = function() {
    stitch = load_stitch_data()
    export_dti(stitch, "stitch_mouse")
}

load_stitch_data = function() {
    result = read.csv(file.path(this_directory__, STITCH_NAME), header = TRUE, comment.char = "", sep = "\t")
    gene_names = read.csv(file.path(this_directory__, ALIAS_NAME), header=TRUE, sep="\t") %>%
        filter(source=="Ensembl_MGI") %>%
        select(protein = `X.string_protein_id`, gene = alias) %>%
        group_by(protein) %>%
        filter(n()==1) %>%
        ungroup()
    result = inner_join(result, gene_names)
    return(result)
}

export_dti = function(stitch_data, file_prefix) {
    result = t(build_dti_from_stitch(stitch_data))
    writeMM(result, paste0(file_prefix, ".txt"))
    write(rownames(result), paste0(file_prefix, "_row_labels.txt"))
    write(colnames(result), paste0(file_prefix, "_col_labels.txt"))
}

build_dti_from_stitch = function(stitch_data, score_name="combined_score", cut_off=0.9) {
    score_cut_off = quantile(stitch_data[[score_name]], probs = cut_off)
    stitch_data = filter(stitch_data, .data[[score_name]] >= score_cut_off)
    chemicals = as.factor(stitch_data$chemical)
    genes = as.factor(stitch_data$gene)
    result = sparseMatrix(as.numeric(genes), as.numeric(chemicals),
                          dims = c(nlevels(genes), nlevels(chemicals)),
                          dimnames = list(levels(genes), levels(chemicals)))
    return(result)
}

if (sys.nframe() == 0) {
    export_mouse_dti()
}

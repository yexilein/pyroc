
library(tidyverse)
library(Matrix)


this_directory__ = getwd()
UNIPROT_NAME = "uniprot_220812.tsv.gz"

export_mouse_uniprot = function() {
    uniprot = load_uniprot()
    export_uniprot(uniprot, "codomain_mouse")
}

load_uniprot = function() {
    result = read.csv(file.path(this_directory__, UNIPROT_NAME), header = TRUE, comment.char = "", sep = "\t")
    result = result %>%
        summarize(gene_name = parse_gene_name(Gene.Names),
                  domain = parse_domain(Domain..FT.))
    result = unnest(result, domain)
    return(drop_na(result))
}

parse_gene_name = function(x) {
    sapply(strsplit(x, split=" ", fixed=TRUE), "[", 1)
}

parse_domain = function(x) {
    lapply(str_match_all(x, "/note=([^;]+)"), function(x) {x[,2]})
}

export_uniprot = function(uniprot, file_prefix) {
    result = build_codomain_from_uniprot(uniprot)
    writeMM(result, paste0(file_prefix, ".txt"))
    write(rownames(result), paste0(file_prefix, "_genes.txt"))
}

build_codomain_from_uniprot = function(uniprot) {
    domain = as.factor(uniprot$domain)
    genes = as.factor(uniprot$gene_name)
    result = sparseMatrix(as.numeric(genes), as.numeric(domain),
                          dims = c(nlevels(genes), nlevels(domain)),
                          dimnames = list(levels(genes), levels(domain)))
    # TF-IDF + cosine similarity or Jaccard?
    #ji = jaccard_index(result)
    tfidf = tf_idf(result)
    return(tfidf)
}

jaccard_index = function(x) {
    intersect = tcrossprod(1*x)
    nnz = rowSums(x)
    union = outer(nnz, nnz, FUN = "+")
    result = intersect/union
    diag(result)=1
    return(result)
}

tf_idf = function(x) {
    tf = x / rowSums(x)
    idf = log(nrow(x)) - log(colSums(x))
    tfidf = t(tf)/idf
    result = crossprod(tfidf)
    l2_norm = sqrt(colSums(tfidf**2))
    result = result / outer(l2_norm, l2_norm)
    return(result)
}

if (sys.nframe() == 0) {
    export_mouse_uniprot()
}

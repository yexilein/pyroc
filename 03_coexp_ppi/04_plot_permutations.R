
library(tidyverse)
source("../data.R", chdir=TRUE)


main = function() {
    go_terms = read_delim("curves/go_terms.txt", delim="\t")
    stats = read_csv("stats/coexp/stats.csv")
    colnames(stats)[1] = "Curve index"
    stats = inner_join(go_terms, stats)
    genes = read.table("curves/genes.txt")[,1]

    go_term = "GO:0008083|growth factor activity|molecular_function" # 0.135
    #go_term = "GO:0043413|macromolecule glycosylation|biological_process" # 0.4
    i = stats$`Curve index`[stats$`Gene set` == go_term]
    predictor = t(read.table("curves/coexp/predictors.txt", skip = i, nrows = 1))[,1]
    labels = mouse_go()[go_term,genes]
    
    original = compute_roc(predictor, labels)
    
    fec_min=0
    fec_max = 0.135
    n_permutations = 10
    set.seed(17)
    fec_permutations = replicate(n = n_permutations, permute_curve(predictor, labels, fec_max, fec_min), simplify = FALSE)
    fec_permutations = transpose(fec_permutations)
    perm_auroc = enframe(fec_permutations$auroc, name = "permutation_id", value="auroc")
    perm_roc = enframe(fec_permutations$roc, name = "permutation_id") %>%
        unnest(value)
    
    global_permutations = replicate(n = n_permutations, permute_curve(predictor, labels, 1,0), simplify = FALSE)
    global_permutations = transpose(global_permutations)
    global_auroc = enframe(global_permutations$auroc, name = "permutation_id", value="auroc")
    global_roc = enframe(global_permutations$roc, name = "permutation_id") %>%
        unnest(value)

    original$roc %>%
        ggplot(aes(x=fpr,y=tpr)) +
        geom_line() +
        theme_classic(base_size = 20)
    ggsave("figs/permutations_original_roc.pdf")

    original$roc %>%
        filter(fpr > fec_max | fpr < fec_min) %>%
        ggplot(aes(x=fpr,y=tpr)) +
        geom_line() +
        geom_line(data=filter(original$roc, fpr<=fec_max & fpr>=fec_min)) +
        geom_line(data=filter(perm_roc, fpr <= fec_max & fpr>=fec_min), aes(group=permutation_id), alpha=0.2) +
        theme_classic(base_size = 20)
    ggsave("figs/permutations_full_roc.pdf")

    original$roc %>%
        ggplot(aes(x=fpr,y=tpr)) +
        geom_line() +
        geom_line(data=global_roc, aes(group=permutation_id), alpha=0.2) +
        theme_classic(base_size = 20)
    ggsave("figs/permutations_global_perm.pdf")

    original$roc %>%
        filter(fpr > fec_max) %>%
        ggplot(aes(x=fpr,y=tpr)) +
        geom_line(data=filter(original$roc, fpr<=fec_max)) +
        geom_line(data=filter(perm_roc, fpr <= fec_max), aes(group=permutation_id), alpha=0.2) +
        theme_classic(base_size = 20)
    ggsave("figs/permutations_zoom_roc.pdf")
    

    fec_max = 0.135
    n_permutations = 100
    fec_permutations = replicate(n = n_permutations, permute_curve(predictor, labels, fec_max), simplify = FALSE)
    fec_permutations = transpose(fec_permutations)
    perm_auroc = enframe(fec_permutations$auroc, name = "permutation_id", value="auroc") %>%
        unnest(auroc)
    
    global_permutations = replicate(n = n_permutations, permute_curve(predictor, labels, 1), simplify = FALSE)
    global_permutations = transpose(global_permutations)
    global_auroc = enframe(global_permutations$auroc, name = "permutation_id", value="auroc") %>%
        unnest(auroc)

    all_aurocs = bind_rows(list(fec = perm_auroc, global = global_auroc), .id = "permutation_type")
    all_aurocs %>%
        ggplot(aes(x = auroc, y=..ndensity.., fill = permutation_type)) +
        geom_density() +
        theme_classic(base_size=20) +
        theme(legend.position = c(0.6, 0.8)) +
        geom_vline(xintercept = original$auroc, linetype = "dashed") +
        expand_limits(x=0.8) +
        scale_fill_brewer(palette = "Set1")
    ggsave("figs/permutations_auroc.pdf")
}

compute_roc = function(predictor, labels) {
    pred = ROCR::prediction(predictor, labels)
    roc = ROCR::performance(pred, "tpr", "fpr")
    roc = tibble(fpr = roc@x.values[[1]], tpr = roc@y.values[[1]])
    auc = ROCR::performance(pred, measure = "auc")@y.values[[1]]
    return(list(roc = roc, auroc=auc))
}

permute_curve = function(predictor, labels, fec_max, fec_min=0) {
    n_max = floor(fec_max*length(labels))
    n_min = floor(fec_min*length(labels))
    r = rank(-predictor)
    randomize = r <= n_max & r >=n_min
    new_labels = labels
    new_labels[randomize] = sample(new_labels[randomize])
    return(compute_roc(predictor, new_labels))
}

plot_curve = function(curves, index) {
    to_plot = as_tibble(curves[,index+1,drop=FALSE]) %>%
        mutate(fpr = seq(0,1,l=nrow(curves))) %>%
        pivot_longer(-fpr, names_to="curve_id", values_to="tpr") %>%
        separate(curve_id, c("go_id", "go_name", "go_ontology"), sep = "\\|")
    to_plot %>%
        ggplot(aes(x=fpr, y=tpr, col=go_name)) +
        geom_line() +
        theme_classic(base_size = 20) +
        theme(legend.position=c(0.7, 0.2)) +
        labs(col=NULL)
}

if (sys.nframe() == 0) {
    main()
}
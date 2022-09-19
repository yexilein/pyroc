
library(tidyverse)


main = function() {
    set.seed(17)
    all_params = all_sim_params()
    for (i in seq_along(all_params)) {
        net = make_network(all_params[[i]])
        write.table(net, file=paste0("simulated_networks/net", i, ".txt"), row.names=FALSE, col.names=FALSE)
        write(rownames(net), file=paste0("simulated_networks/net", i, "_genes.txt"))        
    }
}

make_network = function(sim_params) {
    total_genes = sum(sim_params$n_genes)
    result = tibble(gene = 1:total_genes, block = rep(1:length(sim_params$n_genes), sim_params$n_genes), dummy=1)
    result = inner_join(result, result, by="dummy", suffix = c("_a", "_b")) %>%
        mutate(mu = ifelse(block_a==block_b, sim_params$mu_within, sim_params$mu_across)) %>%
        mutate(edge = pmax(rnorm(n(), mu, sim_params$sigma),0))
    result = result %>%
        mutate(gene_a = paste(gene_a, block_a, sep="|"), gene_b = paste(gene_b, block_b, sep="|")) %>%
        pivot_wider(gene_a, names_from = gene_b, values_from = edge) %>%
        column_to_rownames("gene_a") %>%
        as.matrix()
    result = result[,rownames(result)]
    return(result)
}

all_sim_params = function() {
    list(
        list(n_genes=c(2500,2500,2500,2500), sigma=1, mu_within=2, mu_across=2),
        list(n_genes=c(2500,2500,2500,2500), sigma=1, mu_within=2.1, mu_across=2),
        list(n_genes=c(2500,2500,2500,2500), sigma=1, mu_within=2.2, mu_across=2),
        list(n_genes=c(2500,2500,2500,2500), sigma=1, mu_within=2.3, mu_across=2)
    )
}

if (sys.nframe() == 0) {
    main()
}
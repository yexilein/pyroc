
library(tidyverse)
library(ROCR)


main = function() {
    all_params = sim_params()
    sigma_test = c(0.1, 0.2, 0.3, 0.4)
    pdf("figs/gaussian_simulation.pdf")
    for (p in all_params) {
        my_genes = make_gene_params(p)
        ref_study = simulate_functionality(my_genes) %>%
            mutate(is_positive = (rank(-functionality) / n()) < my_sim_params$threshold)
        t = quantile(ref_study$functionality, probs=1-my_sim_params$threshold)
        all_plots = list(plot_functionality(ref_study, "mu") +
                         geom_vline(xintercept = t, linetype="dashed") +
                         scale_fill_brewer(palette = "Purples"))
        for (st in sigma_test) {
            my_genes$sigma = st
            test_study = simulate_functionality(my_genes) %>%
                mutate(is_positive = ref_study$is_positive)
            all_plots[[length(all_plots)+1]] = plot_functionality(test_study)
            roc = make_roc(test_study$functionality, test_study$is_positive)
            all_plots[[length(all_plots)+1]] = plot_roc(roc)
        }
        do.call(gridExtra::grid.arrange, all_plots)
    }
    dev.off()
}

simulate_functionality = function(gene_params) {
    result = gene_params %>%
        mutate(functionality = rnorm(n(), mu, sigma))
}

make_gene_params = function(sim_params) {
    gaussian_params = tibble(gaussian_id=1:length(sim_params$w_gaussians),
                             mu = (gaussian_id-1)/(max(gaussian_id)-1),
                             sigma = sim_params$sigma_noise)
    gaussian_params$mu = ifelse(is.na(gaussian_params$mu), 0, gaussian_params$mu)
    result = tibble(gene_id = 1:sim_params$n_genes,
                    gaussian_id = rep(1:length(sim_params$w_gaussians),
                                      rmultinom(1, sim_params$n_genes, sim_params$w_gaussians))) %>%
        inner_join(gaussian_params) %>%
        select(-gaussian_id)
    return(result)
}

default_sim_params = function() {
    list(n_genes=10000, w_gaussians = c(0.9,0.1), sigma_noise = 0.5, threshold = 0.1)
}

sim_params = function() {
    list(
        list(n_genes=10000, w_gaussians=c(1), sigma_noise = 1, threshold = 0.2),
        list(n_genes=10000, w_gaussians=c(0.8, 0.2), sigma_noise = 1/2, threshold = 0.2),
        list(n_genes=10000, w_gaussians=c(0.5, 0.3, 0.2), sigma_noise = 1/3, threshold = 0.2),
        list(n_genes=10000, w_gaussians=c(0.4,0.3,0.2,0.1), sigma_noise = 1/4, threshold = 0.2)
    )
}

plot_functionality = function(fn_table, fill_col = "is_positive") {
    fn_table %>%
        ggplot(aes(x=functionality, y=after_stat(count), fill = as.factor(.data[[fill_col]]))) +
        geom_density(show.legend = FALSE, alpha=0.5) +
        theme_classic() + 
        theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line.y = element_blank()) +
        labs(x="Functionality",y=NULL)
}

make_roc = function(predictor, labels) {
    pred = ROCR::prediction(predictor, labels)
    result = ROCR::performance(pred, "tpr", "fpr")
    result = tibble(tpr = unlist(result@x.values), fpr = unlist(result@y.values))
    return(result)
}

plot_roc = function(roc_curve) {
    ggplot(roc_curve) +
        geom_line(aes(x=tpr,y=fpr)) +
        theme_void() +
        theme(panel.border = element_rect(fill=NA,size=1))
}
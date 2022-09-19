
library(tidyverse)
theme_set(theme_classic(base_size = 20))


main = function() {
    plot_curves()
}

plot_curves = function() {
    all_curves = lapply(set_names(c("net1", "net2", "net3", "net4")), read_curves)
    all_curves = bind_rows(all_curves, .id="modularity")
    all_curves_tidy = all_curves %>%
        pivot_longer(c(-modularity, -fpr), names_to = "state", values_to = "tpr")

    p1 = ggplot(all_curves_tidy, aes(x=fpr, y=tpr)) +
        geom_line() +
        labs(col=NULL, x="FPR",y="TPR") +
        facet_grid(modularity ~ state) +
        theme_void() +
        theme(panel.border = element_rect(fill=NA,size=1))
    ggsave("figs/network simulation.pdf", p1)
}

read_curves = function(prefix) {
    as.data.frame(t(read.table(paste0("curves/", prefix, "_curves.txt")))) %>%
        set_names(c("state_1", "state_2","state_3","state_4")) %>%
        mutate(fpr = seq(0,1,length.out = nrow(.)))
}

plot_curve = function(curves) {
    to_plot = as_tibble(curves[,index+1,drop=FALSE]) %>%
        mutate(fpr = seq(0,1,l=nrow(curves))) %>%
        pivot_longer(-fpr, names_to="curve_id", values_to="tpr") %>%
        separate(curve_id, c("go_id", "go_name", "go_ontology"), sep = "\\|")
    result = ggplot(to_plot, aes(x=fpr, y=tpr)) +
        labs(col=NULL, x="FPR",y="TPR")
    return(result)
}

get_tpr = function(curves, index, fpr) {
    curves[200*fpr+1,index+1]
}

get_name = function(curves, index) {
    colnames(curves)[index+1]
}

if (sys.nframe() == 0) {
    main()
}

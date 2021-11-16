
library(tidyverse)
library(ggupset)
library(MetaMarkers)
library(SingleCellExperiment)
theme_set(theme_classic(base_size = 20))


main = function() {
    make_fecs()
    plot_panels()
}

make_fecs = function() {
    markers = read_markers("de_lists/markers_tm.csv.gz") %>%
        filter(population_size >= 20) %>%
        group_by(group, cell_type, gene) %>%
        filter(min(detection_rate) > 0.1) %>%
        ungroup()

    my_cell_type = "B cell"
    my_train_mouse = "3_10_M_ss"
    
    my_train_tissue = "Fat"
    fat20_panel = select_initial_panel(markers, my_cell_type, my_train_tissue, my_train_mouse, 20)
    length(fat20_panel)
    fat20_rocs = compute_all_rocs(markers, my_cell_type, fat20_panel)
    write(fat20_panel, "fec/fat20_panel.txt")
    write_csv(fat20_rocs$roc, "fec/fat20_roc.csv")
    write_csv(fat20_rocs$auroc, "fec/fat20_auroc.csv")

    my_tissue = "Lung"
    lung_fec = select_fec(markers, my_cell_type, my_tissue, 0.05)
    length(lung_fec$gene)
    sum(fat20_panel %in% lung_fec$gene)
    lung_rocs = compute_all_rocs(markers, my_cell_type, unique(lung_fec$gene))
    write_csv(lung_fec, "fec/lung_fec.txt")
    write_csv(lung_rocs$roc, "fec/lung_roc.csv")
    write_csv(lung_rocs$auroc, "fec/lung_auroc.csv")
    
    my_tissue = "Spleen"
    spleen_fec = select_fec(markers, my_cell_type, my_tissue, 0.01)
    length(spleen_fec$gene)
    spleen_rocs = compute_all_rocs(markers, my_cell_type, unique(spleen_fec$gene))
    write_csv(spleen_fec, "fec/spleen_fec.txt")
    write_csv(spleen_rocs$roc, "fec/spleen_roc.csv")
    write_csv(spleen_rocs$auroc, "fec/spleen_auroc.csv")

    my_tissue = "Fat"
    fat_fec = select_fec(markers, my_cell_type, my_tissue, 0.05)
    fat_rocs = compute_all_rocs(markers, my_cell_type, unique(fat_fec$gene))
    write_csv(fat_fec, "fec/fat_fec.txt")
    write_csv(fat_rocs$roc, "fec/fat_roc.csv")
    write_csv(fat_rocs$auroc, "fec/fat_auroc.csv")
}

select_initial_panel = function(markers, my_cell_type, my_train_tissue, my_train_mouse, top_n) {
    top_markers = markers %>%
        filter(cell_type == my_cell_type & group == my_train_tissue & mouse_id == my_train_mouse) %>%
        filter(rank(-auroc) <= top_n) %>%
        pull(gene)
}

compute_all_rocs = function(markers, my_cell_type, top_markers) {
    predictors = markers %>%
        filter(cell_type == my_cell_type) %>%
        mutate(mouse_id = paste0(mouse_id, "|", group)) %>%
        select(test_mouse = mouse_id, gene, auroc) %>%
        pivot_wider(everything(), names_from = test_mouse, values_from = auroc) %>%
        drop_na()
    labels = predictors$gene %in% top_markers
    rocs = apply(as.matrix(predictors[,-1]), 2, compute_roc, labels) %>%
        bind_rows(.id = "test_mouse") %>%
        separate(test_mouse, c("mouse_id", "tissue"), "\\|") %>%
        add_column(test_ct = my_cell_type)
    aurocs = apply(as.matrix(predictors[,-1]), 2, compute_auroc, labels) %>%
        bind_rows(.id = "test_mouse") %>%
        separate(test_mouse, c("mouse_id", "tissue"), "\\|") %>%
        add_column(test_ct = my_cell_type)
    return(list(roc = rocs, auroc = aurocs))
}

compute_roc = function(predictor, labels) {
    roc = ROCR::prediction(predictor, labels)
    roc = ROCR::performance(roc, "tpr", "fpr")
    result = data.frame(
        fpr = unlist(roc@x.values),
        tpr = unlist(roc@y.values),
        score = unlist(roc@alpha.values)
    )
    return(result)
}

compute_auroc = function(predictor, labels) {
    roc = ROCR::prediction(predictor, labels)
    result = ROCR::performance(roc, "auc")
    return(data.frame(auroc = unlist(result@y.values)))
}

select_fec = function(markers, my_cell_type, my_tissue, fec_size) {
    markers %>%
        filter((group == my_tissue) & cell_type == my_cell_type) %>%
        select(mouse_id, gene, auroc) %>%
        group_by(mouse_id) %>%
        filter(rank(-auroc) / n() < fec_size) %>%
        ungroup()
}

plot_panels = function() {
    rocs = read_csv("fec/fat20_roc.csv")
    tissue_cols = make_tissue_cols(rocs$tissue)
    rocs %>%
        filter(tissue %in% c("Fat", "Lung")) %>%
        plot_rocs() +
        scale_color_manual(values = tissue_cols) +
        annotate("rect", xmin = 0, xmax=0.05, ymin=0, ymax=1, linetype="dashed", fill=NA,color="black")
    ggsave("figs/fat20.pdf")
    
    rocs %>%
        filter(tissue %in% c("Lung")) %>%
        plot_rocs() +
        scale_color_manual(values = tissue_cols) +
        lims(x=c(0,0.05),y=c(0,1))
    ggsave("figs/fat20_zoom.pdf")

    rocs %>%
        filter(tissue %in% c("Lung") & mouse_id =="3-F-57_10x") %>%
        plot_rocs() +
        scale_color_manual(values = tissue_cols) +
        lims(x=c(0,0.05),y=c(0,1))
    ggsave("figs/fat20_example.pdf")
    
    rocs %>%
        filter(tissue %in% c("Lung") & fpr <0.05) %>%
        group_by(mouse_id) %>%
        filter(fpr > 0 | tpr == max(tpr[fpr==0])) %>%
        summarize(compute_ks(fpr, tpr))

    rocs = read_csv("fec/lung_roc.csv")
    fec = read_csv("fec/lung_fec.txt")
    pdf("figs/lung_upset.pdf")
    fec %>%
        group_by(gene) %>%
        summarize(mouse_id = list(sort(mouse_id))) %>%
        ggplot(aes(x = mouse_id)) +
        geom_bar() +
        scale_x_upset(order_by = "degree", reverse = TRUE, n_intersections = 20)
    dev.off()
    rocs %>%
        filter(tissue %in% c("Fat", "Lung", "Spleen")) %>%
        plot_rocs() +
        scale_color_manual(values = tissue_cols) +
        annotate("rect", xmin = 0, xmax=0.01, ymin=0, ymax=0.3, linetype="dashed", fill=NA,color="black")
    ggsave("figs/lung.pdf")
    
    rocs %>%
        filter(tissue %in% c("Spleen")) %>%
        plot_rocs() +
        scale_color_manual(values = tissue_cols) +
        lims(x=c(0,0.01),y=c(0,0.3))
    ggsave("figs/lung_zoom.pdf") 
    
    rocs %>%
        filter(tissue %in% c("Spleen") & fpr<0.01) %>%
        group_by(mouse_id) %>%
        filter(fpr > 0 | tpr == max(tpr[fpr==0])) %>%
        summarize(compute_ks(fpr, tpr)) %>%
        mutate(fdr = p.adjust(pval, "fdr"))

    rocs %>%
        filter(tissue %in% c("Spleen")) %>%
        compute_average_roc() %>%
        filter(fpr<0.01) %>%
        filter(fpr > 0 | tpr == max(tpr[fpr==0])) %>%
        summarize(compute_ks(fpr, tpr))

    rocs = read_csv("fec/spleen_roc.csv")
    rocs %>%
#        filter(tissue %in% c("Fat", "Lung")) %>%
        plot_rocs(summary_only = TRUE) +
        scale_color_manual(values = tissue_cols)
    ggsave("figs/spleen.pdf")

    aurocs = list(
        fat = read_csv("fec/fat20_auroc.csv"),
        lung = read_csv("fec/lung_auroc.csv"),
        spleen = read_csv("fec/spleen_auroc.csv")
    ) %>%
        bind_rows(.id = "fec") %>%
        mutate(fec = factor(fec, levels = c("fat", "lung", "spleen")))
    aurocs %>%
        group_by(tissue, fec) %>%
        summarize(med = median(auroc), q1 = quantile(auroc, p=0.25), q3 = quantile(auroc, p=0.75)) %>%
        ggplot(aes(x = tissue)) +
        geom_line(aes(y=med, col = fec, group=fec), size=1) +
        geom_ribbon(aes(ymin=q1,ymax=q3,fill=fec, group=fec), alpha=0.3, show.legend = FALSE) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x = "Tissue", y = "Marker agreement (AUROC)", col = "Marker set")        
    ggsave("figs/aurocs.pdf")
    
    aurocs %>%
        group_by(tissue, fec) %>%
        summarize(med = median(auroc)) %>%
        pivot_wider("tissue", names_from = "fec", values_from = "med") %>%
        write_csv("figs/median_auroc.csv")
    
    fat = filter(aurocs, fec == "fat") %>%
        mutate(is_10x=1*endsWith(mouse_id, "_10x"))
    total_var = with(fat, var(auroc))
    group_by(fat, tissue) %>% summarize(auroc=mean(auroc)) %>% with(., var(auroc) / total_var)
    group_by(fat, is_10x) %>% summarize(auroc=mean(auroc)) %>% with(., var(auroc) / total_var)
    group_by(fat, mouse_id) %>% summarize(auroc=mean(auroc)) %>% with(., var(auroc) / total_var)
    group_by(fat, tissue) %>% summarize(auroc=mean(auroc))
}

make_tissue_cols = function(tissues) {
    unique_tissues = sort(unique(tissues))
    result = RColorBrewer::brewer.pal(n=length(unique_tissues)+1, name = "Set1")
    result = result[-6]
    names(result) = unique_tissues
    return(result)
}

plot_rocs = function(rocs, summary_only=FALSE) {
    summarized_rocs = rocs %>%
        group_by(tissue) %>%
        group_split() %>%
        map(compute_average_roc) %>%
        bind_rows() %>%
        mutate(mouse_id = paste0("average_", tissue))
    result = ggplot(rocs, aes(x = fpr, y = tpr, col = tissue, group=interaction(mouse_id, tissue)))
    if (!summary_only) {
        result = result + geom_line(size=0.5)
    }
    result +
        geom_line(data = summarized_rocs, size=1.5) +
        theme(legend.position = c(0.7, 0.2)) +
        labs(x="False Positive Rate", y="True Positive Rate", col=NULL)
}

compute_average_roc = function(roc_data) {
    n_mice = length(unique(roc_data$mouse_id))
    result = roc_data %>%
        select(mouse_id, tpr, fpr) %>%
        group_by(mouse_id, fpr) %>%
        summarize(min_tpr = min(tpr), max_tpr = max(tpr)) %>%
        group_by(fpr) %>%
        filter(n() == n_mice) %>%
        summarize(min_tpr = mean(min_tpr), max_tpr = mean(max_tpr)) %>%
        pivot_longer(-fpr, names_to = NULL, values_to = "tpr") %>%
        distinct()
    result$tissue = roc_data$tissue[1]
    return(result)
}

plot_detailed_rocs = function(rocs) {
    rocs %>%
        mutate(technology = ifelse(grepl("10x", mouse_id), "10x", "SmartSeq")) %>%
        mutate(sex = ifelse(grepl("M", mouse_id), "M", "F")) %>%
        ggplot(aes(x = fpr, y = tpr, col = mouse_id, linetype = technology)) +
        geom_line() +
        facet_wrap(~ tissue)
}

compute_ks = function(fpr, tpr) {
    n_positives = length(unique(tpr))-1
    delta_x = (fpr - min(fpr)) / (max(fpr)-min(fpr))
    delta_y = (tpr - min(tpr)) / (max(tpr)-min(tpr))
    local_deviation = abs(delta_y - delta_x)
    ks = max(local_deviation, na.rm = TRUE)
    normalized_ks = ks * sqrt(n_positives)
    pval = ks_pval_approx(ks, n_positives)
    return(data.frame(ks=normalized_ks,pval=pval, n=n_positives))
}

ks_pval_exact = function(ks, n) {
    sapply(ks, function(x)
        1-.Call(stats:::C_pKolmogorov2x,x,n)
    )
}

ks_pval_approx = function(ks, n) {
    sapply(ks, function(x)
        1-pkstwo(sqrt(n)*x)
    )
}

pkstwo <- function(x, tol = 1e-06) {
            if (is.numeric(x)) 
                x <- as.double(x)
            else stop("argument 'x' must be numeric")
            p <- rep(0, length(x))
            p[is.na(x)] <- NA
            IND <- which(!is.na(x) & (x > 0))
            if (length(IND)) 
                p[IND] <- .Call(stats:::C_pKS2, p = x[IND], tol)
            p
        }

if (sys.nframe() == 0) {
    main()
}

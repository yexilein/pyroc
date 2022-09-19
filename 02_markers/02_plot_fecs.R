
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
    length(unique(lung_fec$gene))
    sum(fat20_panel %in% lung_fec$gene)
    lung_rocs = compute_all_rocs(markers, my_cell_type, unique(lung_fec$gene))
    write_csv(lung_fec, "fec/lung_fec.txt")
    write_csv(lung_rocs$roc, "fec/lung_roc.csv")
    write_csv(lung_rocs$auroc, "fec/lung_auroc.csv")
    
    my_tissue = "Spleen"
    spleen_fec = select_fec(markers, my_cell_type, my_tissue, 0.01)
    length(unique(spleen_fec$gene))
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
    
    # omnigenic vs polygenic?
    # what happens when we pool FECs?
    pooled_markers = unique(c(fat20_panel, spleen_fec$gene, lung_fec$gene))
    length(pooled_markers)
    pooled_fec_rocs = compute_all_rocs(markers, my_cell_type, pooled_markers)
    write_csv(pooled_fec_rocs$roc, "fec/pooled_markers_roc.csv")
    write_csv(pooled_fec_rocs$auroc, "fec/pooled_markers_auroc.csv")
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
    
    # omnigenic vs polygenic
    rocs = read_csv("fec/spleen_roc.csv")
    n_markers = 50
    # what happens when we pool individuals?
    average_rocs = rocs %>%
        group_by(tissue) %>%
        group_split() %>%
        map(compute_average_roc) %>%
        bind_rows() %>%
        mutate(mouse_id = paste0("average_", tissue)) %>%
        filter(n_curves > 2)
    fecs = average_rocs %>%
        mutate(n_positives = n_markers * n_curves) %>%
        group_by(tissue, mouse_id) %>%
        summarize(find_fecs(fpr, tpr, n_positives[1], 1), .groups="drop")
    p1 = average_rocs %>%
        plot_rocs(summary_only = TRUE) +
        geom_segment(data = filter(fecs, fpr_end-fpr_start>0.05),
                     aes(x=fpr_start,xend=fpr_end,y=tpr_start,yend=tpr_end), col="gray0") +
        scale_color_manual(values = tissue_cols) +
        facet_wrap(~ tissue)
    ggsave("figs/pooled_individuals.pdf", p1)
    p2 = p1 + coord_cartesian(xlim = c(0, 0.1)) + guides(col = "none")
    ggsave("figs/pooled_individuals_zoom.pdf", p2)

    # what happens when we pool tissues?
    rocs %>%
        select(mouse_id, tissue) %>%
        distinct() %>%
        group_by(mouse_id) %>%
        filter(n()>2) %>%
        nest()
    average_rocs = rocs %>%
        group_by(mouse_id) %>%
        group_split() %>%
        map(compute_average_roc, average_over = "tissue", group_var="mouse_id") %>%
        bind_rows() %>%
        mutate(tissue = mouse_id) %>%
        filter(n_curves == 3)
    fecs = average_rocs %>%
        mutate(n_positives = n_markers * n_curves) %>%
        group_by(tissue, mouse_id) %>%
        summarize(find_fecs(fpr, tpr, n_positives[1], 1), .groups="drop")
    p1 = average_rocs %>%
        plot_rocs(summary_only = TRUE) +
        geom_segment(data = mutate(filter(fecs, fpr_end-fpr_start>0.05), tissue = mouse_id),
                     aes(x=fpr_start,xend=fpr_end,y=tpr_start,yend=tpr_end), col="gray0") +
        facet_wrap(~ tissue)
    ggsave("figs/pooled_tissues.pdf", p1)
    p2 = p1 + coord_cartesian(xlim = c(0, 0.1)) + guides(col = "none")
    ggsave("figs/pooled_tissues_zoom.pdf", p2)

    global_average = mutate(rocs, mouse_tissue=paste(mouse_id, tissue), tissue="all") %>%
        compute_average_roc(average_over = "mouse_tissue", group_var = "tissue") %>%
        mutate(mouse_id = "all")
    fecs = global_average %>%
        mutate(n_positives = n_markers * n_curves) %>%
        group_by(tissue, mouse_id) %>%
        summarize(find_fecs(fpr, tpr, n_positives[1], 1), .groups="drop")
    global_average %>%
        plot_rocs(summary_only = TRUE) +
        geom_segment(data = mutate(filter(fecs, fpr_end-fpr_start>0.05), tissue = mouse_id),
                     aes(x=fpr_start,xend=fpr_end,y=tpr_start,yend=tpr_end), col="gray0")
    ggsave("figs/global_average.pdf")
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

compute_average_roc = function(roc_data, average_over = "mouse_id", group_var = "tissue") {
    n_curves = length(unique(roc_data[[average_over]]))
    result = roc_data %>%
        select(.data[[average_over]], tpr, fpr) %>%
        group_by(.data[[average_over]], fpr) %>%
        summarize(min_tpr = min(tpr), max_tpr = max(tpr)) %>%
        group_by(fpr) %>%
        filter(n() == n_curves) %>%
        summarize(min_tpr = mean(min_tpr), max_tpr = mean(max_tpr)) %>%
        pivot_longer(-fpr, names_to = NULL, values_to = "tpr") %>%
        distinct()
    result[[group_var]] = roc_data[[group_var]][1]
    result$n_curves = n_curves
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

compute_ks = function(fpr, tpr, n_positives = length(unique(tpr))-1) {
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

find_fecs = function(fpr, tpr, n_positives, max_ks = 1) {
    ks = compute_all_ks(fpr, tpr, n_positives)
    find_linear_segments(fpr, tpr, ks, max_ks)
}

compute_all_ks = function(fpr, tpr, n_positives, min_genes=20) {
    N = length(fpr)
    result = matrix(NA, N, N)
    delta_x = t(outer(fpr, fpr, "-"))
    delta_y = t(outer(tpr, tpr, "-"))
    for (start in seq(1, N-min_genes+1)) {
        for (end in seq(start+min_genes-1, N)) {
            scaled_x = delta_x[start, start:end] / delta_x[start, end]
            scaled_y = delta_y[start, start:end] / delta_y[start, end]
            local_deviation = abs(scaled_y - scaled_x)
            local_deviation[is.na(local_deviation)] = 0
            ks = max(local_deviation, na.rm = TRUE)
            result[start, end] = ks * sqrt(n_positives*delta_y[start, end])
        }
    }
    return(result)
}

find_linear_segments = function(fpr, tpr, ks, max_ks) {
    valid_segments = ks <= max_ks
    valid_segments[is.na(valid_segments)] = FALSE
    segment_length = t(outer(fpr, fpr, "-"))
    segment_length[!valid_segments] = 0
    segments = c()
    while (max(segment_length) > 0) {
        start_end = arrayInd(which.max(segment_length), dim(segment_length))
        segments = rbind(segments, start_end)
        segment_length[,start_end[1]:start_end[2]] = 0
        segment_length[start_end[1]:start_end[2],] = 0
    }
    colnames(segments) = c("start", "end")
    result = as_tibble(segments) %>%
        rowwise() %>%
        mutate(ks = ks[start, end], fpr_start = fpr[start], fpr_end = fpr[end], tpr_start = tpr[start], tpr_end = tpr[end])
    return(result)
}

if (sys.nframe() == 0) {
    main()
}

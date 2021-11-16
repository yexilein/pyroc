
library(tidyverse)
theme_set(theme_classic(base_size = 20))


main = function() {
    plot_global_panels()
    plot_segment_panels()
}

plot_global_panels = function() {
    coexp = read_stats("coexp")
    ppi = read_stats("ppi")
    all_stats = bind_rows(list(COEXP=coexp,PPI=ppi), .id="modality")

    all_stats %>% group_by(modality) %>% summarize(auroc=median(auroc))
    all_stats %>%
        ggplot(aes(x=auroc, fill=modality)) +
        geom_histogram(position = "identity", alpha=0.5) +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Functional enrichment (AUROC)", y="# GO terms",fill="Modality") +
        theme(legend.position=c(0.2,0.8))
    ggsave("panels/auroc.pdf")

    all_stats %>% group_by(modality) %>% summarize(f=mean(fraction_linear))
    all_stats %>% group_by(modality) %>% summarize(f=sum(fraction_linear==0),n=n())
    all_stats %>% group_by(modality) %>% summarize(f=sum(fraction_linear==1),n=n(), ff=mean(fraction_linear==1))
    all_stats %>%
        ggplot(aes(x=fraction_linear, fill=modality)) +
        geom_histogram(position = "identity", alpha=0.5) +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Fraction of ROC curve detected as linear", y="# GO terms",fill="Modality") +
        theme(legend.position=c(0.2,0.8))
    ggsave("panels/f_linear.pdf")

    m = round(min(all_stats$auroc)-0.005,2)
    wide_stats = all_stats %>%
        pivot_wider(id_cols = "GO term", names_from = "modality", values_from = "auroc")
    rho = with(wide_stats, round(cor(COEXP, PPI, method = "spearman"),2))
    wide_stats %>%
        ggplot(aes(x=PPI, y=COEXP)) +
        geom_point(size=0.2, alpha=0.5) +
        annotate("text", label=paste0("rho==", rho),x=0.5, y=0.95, size=8, parse = T) +
        labs(x = "Enrichment PPI (AUROC)", y="Enrichment co-expression (AUROC)") +
        lims(x=c(m,1), y=c(m,1)) +
        geom_abline(slope = 1, linetype="dashed")
    ggsave("panels/ppi_vs_coexp.pdf")
    
    go = read_go()
    curves_ppi = read_curves("ppi", go)
    curves_coexp = read_curves("coexp", go)
    top_ppi = ppi %>%
        filter(`Number of genes`>50 & `Number of genes`<100 & auroc>0.7) %>%
        filter(number_null_segments==2) %>%
        arrange(desc(fraction_linear)) %>%
        head(4)
    top_coexp = coexp %>%
        filter(`Number of genes`>50 & `Number of genes`<100 & auroc>0.7) %>%
        filter(number_null_segments==2) %>%
        arrange(desc(fraction_linear)) %>%
        head(4)
    segments_ppi = read_segments("ppi/")
    segments_coexp = read_segments("coexp/")
    plot_curve(curves_ppi, top_ppi$`Curve index`, segments=segments_ppi) +
        facet_wrap(~ go_name) +
        theme(strip.background = element_blank(), strip.text.x = element_blank()) +
        coord_fixed() +
        ggtitle("PPI")
    ggsave("panels/ppi_example.pdf")
    plot_curve(curves_coexp, top_coexp$`Curve index`, segments=segments_coexp) +
        facet_wrap(~ go_name) +
        theme(strip.background = element_blank(), strip.text.x = element_blank()) +
        coord_fixed() +
        ggtitle("COEXP")
    ggsave("panels/coexp_example.pdf")
}

read_stats = function(prefix) {
    go_terms = read_go()
    stats = read_csv(paste0("stats/", prefix, "/stats.csv"))
    colnames(stats)[1] = "Curve index"
    stats = inner_join(go_terms, stats)
    return(stats)
}

read_go = function() {
    read_delim("curves/go_terms.txt", delim="\t")
}

read_curves = function(prefix, go_terms) {
    as.data.frame(t(read.table(paste0("curves/", prefix, "/curves.txt")))) %>%
        set_names(go_terms$`GO term`) %>%
        mutate(fpr = seq(0,1,length.out = nrow(.)))
}

read_segments = function(prefix) {
    read_delim(paste0("stats/", prefix, "all_linear_segments.txt"), delim=" ", col_names = c("curve_id", "start", "end")) %>%
        mutate(start = start/200, end = end/200) %>%
        mutate(length = end-start)
}

modality_cols = function() {
    c(COEXP="darkblue",PPI="darkorange")
}

plot_curve = function(curves, index, flip=FALSE, segments=NULL) {
    to_plot = as_tibble(curves[,index+1,drop=FALSE]) %>%
        mutate(fpr = seq(0,1,l=nrow(curves))) %>%
        pivot_longer(-fpr, names_to="curve_id", values_to="tpr") %>%
        separate(curve_id, c("go_id", "go_name", "go_ontology"), sep = "\\|")
    if (flip) {
        to_plot = mutate(to_plot, fpr=1-fpr, tpr=1-tpr)
    }
    result = ggplot(to_plot, aes(x=fpr, y=tpr))
    if (!is.null(segments)) {
        segments = segments %>%
            filter(curve_id %in% index) %>%
            rowwise() %>%
            mutate(tpr_start = get_tpr(curves, curve_id, start),
                   tpr_end = get_tpr(curves, curve_id, end),
                   go_id = get_name(curves, curve_id)) %>%
            separate(go_id, c("go_id", "go_name", "go_ontology"), sep = "\\|")
        result = result +
            geom_segment(data=segments, aes(x=start,xend=end,y=tpr_start,yend=tpr_end), col="gray40")
    }
    result = result +
        geom_line(aes(col=go_name)) +
        theme(legend.position=c(0.6, 0.2)) +
        labs(col=NULL, x="False Positive Rate",y="True Positive Rate")
    return(result)
}

get_tpr = function(curves, index, fpr) {
    curves[200*fpr+1,index+1]
}

get_name = function(curves, index) {
    colnames(curves)[index+1]
}

plot_segment_panels = function() {
    coexp = read_stats("coexp")
    ppi = read_stats("ppi")
    all_stats = bind_rows(list(COEXP=coexp,PPI=ppi), .id="modality")

    all_stats %>%
        group_by(modality) %>%
        summarize(f=sum(number_linear_segments==2),n=n(), ff=mean(number_linear_segments==2))
    all_stats %>%
        ggplot(aes(x=as.factor(number_linear_segments), fill=modality)) +
        geom_bar(position="dodge") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "# Linear segments", y="# GO terms",fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("panels/n_segments.pdf")
    
    segments = bind_rows(list(
        PPI = read_segments("ppi/"),
        COEXP = read_segments("coexp/")
    ), .id="modality")

    segments %>%
        ggplot(aes(x=length, fill=modality)) +
        geom_histogram(position="identity", alpha=0.5) +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Linear segment length", y="# Segments", fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("panels/segment_length.pdf")
    
    segments %>% group_by(modality) %>% summarize(short=mean(length<=0.4), long=mean(length>0.4))
    segments %>%
        group_by(modality, curve_id) %>%
        filter(start == min(start) & end<=0.4) %>%
        pull(length) %>%
        quantile(probs=seq(0,1,by=0.1))
    segments %>%
        group_by(modality, curve_id) %>%
        filter(start == min(start)) %>%
        ggplot(aes(x=length, fill=modality)) +
        geom_histogram(position="identity", alpha=0.5) +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Length of initial linear segment", y="# GO Terms", fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("panels/initial_segment_length.pdf")
    
    go = read_go()
    curves_coexp = read_curves("coexp", go)
    top_coexp = coexp %>%
        filter(`Number of genes`>50 & `Number of genes`<100 & auroc>0.7) %>%
        arrange(desc(fraction_linear)) %>%
        head(1)
    top_coexp
    fpr_0 = 200-173+1
    tpr_0 = curves_coexp[fpr_0,top_coexp$`Curve index`+1]
    plot_curve(curves_coexp, top_coexp$`Curve index`) +
        annotate("segment", x=0,xend=fpr_0/200,y=0,yend=tpr_0) +
        annotate("segment", x=fpr_0/200,xend=1,y=tpr_0,yend=1)
    ggsave("panels/two_segment_example.pdf")    

    all_stats %>%
        mutate(x=(auroc-reduced_curve_auroc)/auroc) %>%
        summarize(f = mean(abs(x)<=0.05, na.rm = TRUE))
    all_stats %>%
        ggplot(aes(x=(auroc-reduced_curve_auroc)/auroc, fill=modality)) +
        geom_histogram(position="identity", alpha=0.5) +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "2-segment approximation error (% original AUROC)", y="# GO terms",fill="Modality") +
        theme(legend.position=c(0.8,0.8)) +
        scale_x_continuous(labels = scales::percent)
    ggsave("panels/two_segment_approximation.pdf")
    
    all_stats %>%
        mutate(x=(auroc-reduced_curve_auroc)/auroc) %>%
        mutate(auroc = cut(auroc, breaks = seq(0,1,by = 0.1))) %>%
        group_by(auroc) %>%
        summarize(f = mean(abs(x)<=0.05, na.rm = TRUE))
    all_stats %>%
        mutate(error=(auroc-reduced_curve_auroc)/auroc) %>%
        mutate(auroc = cut(auroc, breaks = seq(0,1,by = 0.1))) %>%
        ggplot(aes(x=auroc, y=error, fill=modality)) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_manual(values = modality_cols()) +
        labs(y = "Approximation error", x="Original AUROC") +
        scale_y_continuous(labels = scales::percent) +
        geom_hline(yintercept = 0.05, linetype="dashed")
    ggsave("panels/two_segment_approximation_inset.pdf")

    all_stats %>%
        filter(modality=="COEXP") %>%
        ggplot(aes(x=auroc, y=reduced_curve_auroc)) +
        geom_point(size=0.2, alpha=0.2, color=modality_cols()["COEXP"]) +
        labs(x = "Original AUROC", y="2-segment AUROC") +
        geom_abline(slope=1, linetype="dashed")
    ggsave("panels/two_segment_coexp.pdf")

    all_stats %>%
        filter(modality=="PPI") %>%
        ggplot(aes(x=auroc, y=reduced_curve_auroc)) +
        geom_point(size=0.2, alpha=0.2, color=modality_cols()["PPI"]) +
        labs(x = "Original AUROC", y="2-segment AUROC") +
        geom_abline(slope=1, linetype="dashed")
    ggsave("panels/two_segment_ppi.pdf")
}

if (sys.nframe() == 0) {
    main()
}
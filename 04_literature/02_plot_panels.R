
library(tidyverse)
theme_set(theme_classic(base_size = 20))


main = function() {
    curves = read_curves()
    stats = read_stats() %>%
        mutate(modality="ARTICLE")
    all_segments = read_segments() %>%
        mutate(modality="ARTICLE")
    
    stats %>% summarize(auroc=median(auroc))
    stats %>%
        ggplot(aes(x=auroc, fill=modality)) +
        geom_histogram(position = "identity") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Functional enrichment (AUROC)", y="# Curves",fill="Modality") +
        theme(legend.position=c(0.2,0.8))
    ggsave("figs/auroc.pdf")

    stats %>% summarize(f=sum(fraction_linear==0),n=n())
    stats %>% summarize(f=mean(fraction_linear))
    stats %>% summarize(f=sum(fraction_linear>0.9),n=n(), f/n)
    stats %>%
        ggplot(aes(x=fraction_linear, fill=modality)) +
        geom_histogram(position = "identity") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Fraction of ROC curve detected as linear", y="# Curves",fill="Modality") +
        theme(legend.position=c(0.2,0.8))
    ggsave("figs/f_linear.pdf")

    table(stats$number_linear_segments)
    stats %>%
        ggplot(aes(x=as.factor(number_linear_segments), fill=modality)) +
        geom_bar(position="dodge") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "# Linear segments", y="# Curves",fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("figs/n_segments.pdf")
    
    all_segments %>%
        ggplot(aes(x=length, fill=modality)) +
        geom_histogram(position="identity") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Linear segment length", y="# Segments", fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("figs/segment_length.pdf")
    
    all_segments %>%
        group_by(modality, curve_id) %>%
        filter(start == min(start)) %>%
        ggplot(aes(x=length, fill=modality)) +
        geom_histogram(position="identity") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "Length of initial linear segment", y="# Curves", fill="Modality") +
        theme(legend.position=c(0.8,0.8))
    ggsave("figs/initial_segment_length.pdf")

    initial_segments = all_segments %>%
        group_by(modality, curve_id) %>%
        filter(start == min(start))
    my_curves = stats %>%
        inner_join(initial_segments, by = c("Curve index"="curve_id")) %>%
        filter(length > 0.2 & length < 0.8) %>%
        slice_max(fraction_linear, n=5) %>%
        pull("Curve index")
    curves %>%
        filter(curve_id %in% my_curves) %>%
        mutate(curve_id = as.factor(curve_id)) %>%
        ggplot(aes(x=fpr, y=tpr, col=curve_id)) +
        geom_line() +
        theme(legend.position=c(0.7, 0.2)) +
        labs(col=NULL,x="False Positive Rate",y="True Positive Rate")
    ggsave("figs/linear_length_example.pdf")
    
    stats %>% mutate(x=(auroc-reduced_curve_auroc)/auroc) %>% summarize(mean(abs(x)<=0.05,na.rm=TRUE))
    stats %>%
        mutate(x=(auroc-reduced_curve_auroc)/auroc) %>%
        mutate(auroc = cut(auroc, breaks = seq(0,1,by = 0.1))) %>%
        group_by(auroc) %>%
        summarize(f = mean(abs(x)<=0.05, na.rm = TRUE))
    stats %>%
        ggplot(aes(x=(auroc-reduced_curve_auroc)/auroc, fill=modality)) +
        geom_histogram(position="identity") +
        scale_fill_manual(values = modality_cols()) +
        labs(x = "2-segment approximation error (% original AUROC)", y="# Curves",fill="Modality") +
        theme(legend.position=c(0.8,0.8)) +
        scale_x_continuous(labels = scales::percent)
    ggsave("figs/two_segment_approximation.pdf")

    stats %>%
        ggplot(aes(x=auroc, y=reduced_curve_auroc)) +
        geom_point(color=modality_cols()["ARTICLE"]) +
        labs(x = "Original AUROC", y="2-segment AUROC") +
        geom_abline(slope=1, linetype="dashed")
    ggsave("figs/two_segment_article.pdf")
    
    my_curves = stats %>%
        filter(number_linear_segments == 2) %>%
        slice_max(fraction_linear, n=5) %>%
        pull("Curve index")
    curves %>%
        filter(curve_id %in% my_curves) %>%
        mutate(curve_id = as.factor(curve_id)) %>%
        ggplot(aes(x=fpr, y=tpr, col=curve_id)) +
        geom_line() +
        theme_classic(base_size = 20) +
        theme(legend.position=c(0.7, 0.2)) +
        labs(col=NULL,x="False Positive Rate",y="True Positive Rate")
    ggsave("figs/two_segment_example.pdf")

    min_auroc = min(stats$auroc)
    stats %>%
        ggplot(aes(x = auroc, y=auroc+flip_performance_increase)) +
        geom_point(color=modality_cols()["ARTICLE"]) +
        geom_abline(slope = 1) +
        geom_abline(slope = 1, intercept = 0.1, linetype="dashed") +
        geom_text(data=filter(stats, flip_performance_increase>0.1), aes(label=`Curve index`)) +
        theme_classic(base_size = 20) +
        lims(x=c(min_auroc,1), y=c(min_auroc,1)) +
        labs(x="Original AUROC",y="AUROC after optimal flip",col="Modality")
    ggsave("figs/flip.pdf")
    
    my_curves = stats %>%
        slice_max(flip_performance_increase, n=4) %>%
        pull("Curve index")
    normal_curves = filter(curves, curve_id %in% my_curves) %>%
        mutate(curve_id = as.factor(curve_id))
    flipped_curves = lapply(my_curves, function(my_curve) {
        my_stats = filter(stats, `Curve index` == my_curve)
        flip_curve(filter(curves, curve_id == my_curve), my_stats$flip_start, my_stats$flip_end)
    })
    flipped_curves = bind_rows(flipped_curves) %>%
        mutate(curve_id = as.factor(curve_id))
    normal_curves %>%
        ggplot(aes(x=fpr, y=tpr, col=curve_id)) +
        geom_line() +
        geom_line(data=flipped_curves, linetype="dashed") +
        theme_classic(base_size = 20) +
        theme(legend.position=c(0.7, 0.2)) +
        labs(col=NULL,x="False Positive Rate",y="True Positive Rate")
    ggsave("figs/flip_example.pdf")
}

read_stats = function() {
    stats = read_csv(paste0("stats/stats.csv"))
    colnames(stats)[1] = "Curve index"
    return(stats)
}

read_curves = function() {
    as.data.frame(t(read.table(paste0("curves/curves.txt")))) %>%
        set_names(as.character(1:ncol(.))) %>%
        mutate(fpr = seq(0,1,length.out = nrow(.))) %>%
        pivot_longer(-fpr, names_to = "curve_id", values_to = "tpr") %>%
        mutate(curve_id = as.numeric(curve_id)-1)
}

read_segments = function() {
    read_delim(paste0("stats/all_linear_segments.txt"), delim=" ", col_names = c("curve_id", "start", "end")) %>%
        mutate(start = start/200, end = end/200) %>%
        mutate(length = end-start)
}

modality_cols = function() {
    c(COEXP="darkblue",PPI="darkorange",ARTICLE="darkgreen")
}

flip_curve = function(curve, flip_start, flip_end) {
    result = curve$tpr
    m = result[flip_start+1]
    M = result[flip_end+1]
    flipped_segment = rev(result[(flip_start+1):(flip_end+1)])
    flipped_segment = m + M - flipped_segment
    curve$tpr[(flip_start+1):(flip_end+1)] = flipped_segment
    return(curve)
}

if (sys.nframe()==0) {
    main()
}
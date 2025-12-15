# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(qgraph)
library(bootnet)
library(flextable)
library(officer)
library(ggthemes)

sect_properties <- prop_section(
    page_size = page_size(
        orient = "landscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0)
)


# Data I/O ---------------------------------------------------------------------
master_srs_df <- read_rds(
    file = here("outputs", "caches", "master_srs_df.rds")
) %>%
    group_by(src_subject_id) %>%
    slice(1) %>%
    ungroup()

discovery_sig_labels <- read_rds(
    file = here("outputs", "caches", "discovery_sig_labels.rds")
)

discovery_df <- filter(master_srs_df, dataset == "Discovery")
holdout_df <- filter(master_srs_df, dataset == "Holdout")

n_label <- length(discovery_sig_labels)

network_data <- list()
network_data$data <- holdout_df %>%
    select(all_of(discovery_sig_labels), avg_srs) %>% 
    rename(`average stressor reactivity score` = avg_srs)
network_data$node <- c(
    paste0("ENV", str_pad(1:n_label, width = 2, pad = 0)), "aSRS"
)
network_data$node_col <- c("#DAC9A6", "#B5CAA0")
network_data$group <- c(
    rep("1. Environment at Baseline", n_label),
    "2. Resilience at Follow-ups"
)
network_data$node_detail <- c(
    discovery_sig_labels,
    "average stressor reactivity score"
)
set.seed(1234)
g <- estimateNetwork(
    network_data$data,
    default = "ggmModSelect",
    stepwise = TRUE,
    corMethod = "cor_auto"
)

network_graph <- qgraph(
    g$graph,
    # layout
    layout = "spring",
    repulsion = 1.3,
    color = network_data$node_col,
    labels = network_data$node,
    groups = network_data$group,
    nodeNames = network_data$node_detail,
    # node
    vsize = 4,
    label.fill.vertical = 0.2,
    borders = FALSE,
    label.cex = 1.25,
    # edge
    esize = 10,
    fade = TRUE,
    # edge curvature
    curve = 0.5,
    curveAll = TRUE,
    # legend
    legend = TRUE,
    legend.cex = 0.4,
    GLratio = 2,
    cut = 0,
    minimum = 0.01,
    theme = "TeamFortress",
    mar = c(2, 2, 2, 2),
    filetype = "pdf",
    filename = here("outputs", "figures", "srs_network"),
    width = 8,
    height = 6
)

# plot centrality
centralityplot <- centralityPlot(
    g, 
    include = c("ExpectedInfluence"),
    orderBy = "ExpectedInfluence"
) +
    scale_x_continuous(limits = c(-1.8, 1.8)) +
    theme_pander()

ggsave(
    plot = centralityplot,
    filename = here("outputs", "figures", "centralityplot.pdf"),
    device = cairo_pdf,
    height = 6,
    width = 6
)


# Check Stability --------------------------------------------------------------
g_bs_np <- bootnet(
    g,
    statistics = c("Edge"),
    nBoots = 1000,
    nCores = 10,
    type = "nonparametric"
)

write_rds(g_bs_np, here("outputs", "caches", "g_bs_np.rds"))
stability_edge_plot <- plot(g_bs_np, statistics = "Edge", order = "sample")

stability_edge_plot$data %>%
    filter(node2 == "average stressor reactivity score") %>%
    mutate(sig = if_else(CIlower * CIupper >= 0, 1, 0)) %>%
    filter(sig == 1) %>%
    slice(1:4) %>%
    select(3, value, 11, 6:7) %>%
    mutate_if(is.numeric, ~ round(., 3)) %>%
    mutate(
        node1 = str_to_title(node1)
    ) %>%
    rename_with(
        ~ c("Connecting Node with SRS", "Empirical", "Bootstrapped", "Lower 95% CI",
            "Upper 95% CI")
    ) %>%
    flextable() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(j = 2:5, align = "center", part = "all") %>%
    save_as_docx(
        path = here("outputs", "tables", "table_edge_stability.docx"),
        pr_section = sect_properties
    )

g_bs_case <- bootnet(
    g,
    statistics = c("expectedInfluence"),
    nBoots = 1000,
    nCores = 6,
    type = "case"
)
write_rds(
    g_bs_case,
    here("outputs", "caches", "g_bs_case.rds")
)
ei_drop_plot <- plot(g_bs_case, statistics = "expectedInfluence") +
    scale_y_continuous(
        breaks = c(0.9, 0.95, 1),
        limits = c(0.9, 1)
    ) +
    labs(
        x = "Sampled Case",
        y = "Averaged Corelation with Original Sample",
        color = "",
        fill = ""
    ) +
    scale_color_manual(
        values = c("expectedInfluence" = "gray30"),
        labels = c("expectedInfluence" = "Expected Influence")
    ) +
    scale_fill_manual(
        values = c("expectedInfluence" = "gray30"),
        labels = c("expectedInfluence" = "Expected Influence")
    ) +
    theme_pander() +
    theme(
        legend.position = "top",
        plot.margin = margin(5, 5, 5, 5, "mm"),
        text = element_text(family = "Helvetica Neue")
    )
ei_drop_plot
ggsave(
    plot = ei_drop_plot,
    filename = here("outputs", "figures", "ei_drop_plot.pdf"),
    width = 6,
    height = 4,
    device = cairo_pdf
)
corStability(g_bs_case)

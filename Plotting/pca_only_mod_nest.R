library(xtable)
library(igraph)
library(scales)
library(stringr)
library(tidyverse)
library(broom)

source("ggbiplot.R")
theme_set(theme_bw())

source("ancillary_functions.R")

BIPLOT_HEIGHT <- 4
BIPLOT_WIDTH <- 6

set.seed(0) # to replicate sampling

#### KEY PARAMETERS
MIN_SIZE <- 5 # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)
results_file <- "../Results/full_results.csv"

#### scales for ggplots
cols <- c("antagonism"              = "#7570b3",
          "mutualism"               = "#66a61e")
shps <- c("antagonism"              = 4,
          "mutualism"               = 3)

#### read in the results file
results_df <- read_csv(results_file,
                       col_types=cols(nrows = col_double(),
                                      feature1 = col_character(),
                                      feature2 = col_character(),
                                      feature3 = col_character(),
                                      feature4 = col_character())) %>%
    mutate(type = tolower(type)) %>% 
    filter(type %in% c("antagonism", "mutualism", "ecologicalinteractions")) %>% 
    # bind_rows(read_csv("../Results/SimpleDemo_old_results.csv")) %>%
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>%
    mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature1), tolower(type)))
    

Metrics <- list(
    "Q"                  = ~Q,
    "N.olap"             = ~N.olap,
    "N.temp"             = ~N.temp,
    "N.nodf"             = ~N.nodf
)

#### run the PCA
to_fit_df <- results_df %>% filter(randomization == "None")
pca_df <- to_fit_df %>% transmute_(.dots=Metrics)
pca_results <- prcomp(pca_df, center=TRUE, scale=TRUE)

#### TABLE OF PCA LOADINGS
pca_results$rotation %>% rbind(pca_results$sdev^2 / sum(pca_results$sdev^2) * 100) %>% xtable(digits=3) %>% print()

#### plot the fit data
g_fit <- ggbiplot(pca_results, groups=to_fit_df$type, subgroups=to_fit_df$type, ellipse=TRUE, var.axes=FALSE, alpha=0.35) +
    scale_colour_manual(name="Network Type", values=cols) + scale_shape_manual(name="Network Type", values=shps)

emp_df <- results_df %>%
    filter(randomization == "None") %>% 
    mutate_(.dots=Metrics) %>%
    ## center and scale according to the fitting data
    mutate_at(names(Metrics),
              funs((. - (pca_df %>% summarise_all(mean))$.) /
                       (pca_df %>% summarise_all(sd))$.))
## transform into PCA space
emp_df <- emp_df %>%
    rowwise() %>%
    ## the ggbiplot code divides the coordinates by the standard deviation
    ## in the pca results. I don't know why, but to make the points line up
    ## with the ellipses, we do so here as well.
    do(as.numeric(.[names(Metrics)]) %*% pca_results$rotation %>%
           sweep(2, pca_results$sdev, '/') %>% 
           tbl_df()) %>%
    ungroup() %>% 
    bind_cols(emp_df)
pca_table_plot <- function(pcs_with_labels, n_pcs=5) {
    ellipses_df <- pcs_with_labels %>%
        select(type, matches(str_c("PC", 1:n_pcs, "$", collapse="|"))) %>%
        bind_cols(rename_at(., vars(starts_with("PC")), funs(str_c("x", .))) %>% select(-type)) %>%
        gather("Axis", "Value", matches("^PC\\d+$")) %>% 
        gather("Axis1", "Value1", matches("^xPC\\d+$")) %>%
        filter(as.integer(str_extract(Axis, "\\d+")) > as.integer(str_extract(Axis1, "\\d+"))) %>%
        mutate(Axis1 = str_replace_all(Axis1, "x", "")) %>% 
        group_by_at(vars(-Value, -Value1)) %>%
        do(ellipse_path(., "Value1", "Value")) %>% 
        ungroup()
    pcs_with_labels %>%
        group_by(type) %>%
        do(pca_table_plot_data(., n_pcs)) %>%
        filter(as.integer(str_extract(Axis, "\\d+")) > as.integer(str_extract(Axis1, "\\d+"))) %>%
        ggplot() +
        aes_string(x="Value1", y="Value", colour=names(ellipses_df)[1], shape=names(ellipses_df)[1]) +
        geom_point(alpha=0.5, size=2) +
        geom_path(data=ellipses_df) +
        facet_grid(Axis~Axis1) +
        coord_equal() +
        scale_colour_manual(values=c("antagonism"   = "#e07699",
                                     "mutualism"    = "#95a7d5",
                                     "biogeography" = "#79c8a5")) +
        scale_shape_manual(values=c("antagonism"   = 4,
                                    "mutualism"    = 3,
                                    "biogeography" = 1)) +
        theme(legend.position="bottom",
              legend.title=element_blank(),
              axis.title = element_blank())
}
pca_table_plot(emp_df %>% select(type, matches("PC\\d+")), 4)
ggsave(filename="../Figures/NestMod_PCA.pdf", width=10, height=10.5)

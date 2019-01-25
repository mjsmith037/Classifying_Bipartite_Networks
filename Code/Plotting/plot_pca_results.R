## afterwards, run "mogrify -format png *.pdf"

library(xtable)
library(igraph)
library(scales)
library(stringr)
library(tidyverse)
library(broom)

source("ggbiplot.R")
theme_set(theme_bw())

source("ancillary_functions.R")

BIPLOT_HEIGHT <- 6
BIPLOT_WIDTH <- 6

set.seed(0) # to replicate sampling

#### KEY PARAMETERS
N_FIT <- 100 # size of subsets used for fitting 
MIN_SIZE <- 5 # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)


# wemp <- ""
# results_file <- "../Results/only_new-data_results.csv"
# MAX_SIG <- Inf

# wemp <- ""
# results_file <- "../Results/full_results.csv"
# MAX_SIG <- Inf

# wemp <- "wEmp_"
# results_file <- "../Results/full_results.csv"
# MAX_SIG <- Inf

wemp <- ""
results_file <- "../Results/SimpleDemo_results.csv"
MAX_SIG <- NA

# wemp <- "wEmp_"
# results_file <- "../Results/SimpleDemo_results.csv"
# MAX_SIG <- NA

# wemp <- "wNewEmp_"
# results_file <- "../Results/SimpleDemo_results.csv"
# MAX_SIG <- NA
# N_FIT <- 48




#### scales for ggplots
cols <- c("crimes"                  = "#1b9e77",
          "antagonism"              = "#e07699",
          "mutualism"               = "#95a7d5",
          "microbiome"              = "#a6761d",
          "actor collaboration"     = "#d95f02",
          "biogeography"            = "#79c8a5",
          "legislature"             = "#e6ab02",
          "authorship"              = "#666666")
shps <- c("crimes"                  = 0,
          "antagonism"              = 4,
          "mutualism"               = 3,
          "microbiome"              = 2,
          "actor collaboration"     = 7,
          "biogeography"            = 1,
          "legislature"             = 5,
          "authorship"              = 6)


#### read in the results file
results_df <- read_csv(results_file,
                       ## specify type for these columns because they vary between networks or cause errors in parsing
                       col_types=cols(nrows = col_double(),
                                      # H2 = col_double(),
                                      # H3 = col_double(),
                                      # H4 = col_double(),
                                      feature1 = col_character(),
                                      feature2 = col_character(),
                                      feature3 = col_character(),
                                      feature4 = col_character()))

if (!str_detect(wemp, "New") & str_detect(results_file, "SimpleDemo")) {
    results_df <- results_df %>%
        bind_rows(read_csv("../Results/SimpleDemo_old_results.csv"))
}

results_df <- results_df %>%
    ## clean up some names
    mutate(type = tolower(type),
           type = str_replace(type, "movies", "actor collaboration"),
           type = str_replace(type, "ecologicalinteractions", "ecological interactions")) %>%
    mutate(type = ifelse(type == "ecological interactions", tolower(feature1), type)) %>%
    mutate(type = str_replace(type, "fungi|islands|mountains", "biogeography"),
           type = str_replace(type, "anemonefish|ant-plant|pollination|seeddispersal", "mutualism"),
           type = str_replace(type, "bacteria-phage|host-parasitoid|parasitism|plant-herbivore", "antagonism")) %>%
    ## remove very small networks
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>% 
    ## remove outliers
    # filter(name != "Washington_DC_crimes_2016-12-25") %>%
    filter(name != "EP_2003")


#### specify/calculate the variables to consider in PCA
available_metrics <- list(
    "CM_Ratio"           = ~l1_cm / l1,
    "ER_Ratio"           = ~l1_er / l1,
    "MP_Ratio"           = ~((1 + sqrt(ncols / nrows)) *
                                 sqrt(nlinks / ncols *
                                          (1 - nlinks / nrows / ncols))) / l2,
    "Reg_Ratio"          = ~l1 / l1_reg,
    "Lap_Ratio"          = ~l1 / l1_lap,
    "Gap_Ratio_21"       = ~(l1 - l2) / l1,
    "Gap_Ratio_32"       = ~(l2 - l3) / l2,
    "Lap_Gap_Ratio_21"   = ~(l1_lap - l2_lap) / l1_lap,
    "Lap_Gap_Ratio_32"   = ~(l2_lap - l3_lap) / l2_lap,
    "Gap_Lap_Gap_Ratio"  = ~(l1_lap - l2_lap) / (l1 - l2),
    "l2l3_Ratio"         = ~l3 / l2,
    "ipr1"               = ~ipr1,
    "ipr2"               = ~ipr2,
    "ipr3"               = ~ipr3,
    "Alg_Connectance"    = ~alg_conn,
    "cent_between"       = ~cent_between,
    "cent_close"         = ~cent_close,
    "cent_eigen"         = ~cent_eigen,
    "av_short_path"      = ~mean_path_length,
    "diam"               = ~diam,
    "log_H2"             = ~log(H2),
    "log_H3"             = ~log(H3),
    "log_H4"             = ~log(H4),
    "log_H17"            = ~log(H17 + 1), # +1 because not all webs have an H17 subgraph
    "H4H2"               = ~H4 / H2,
    "H17H4"              = ~H17 / H4,
    ## "Q"                 = ~Q,
    ## "N.olap"            = ~N.olap,
    ## "N.temp"            = ~N.temp,
    ## "N.nodf"            = ~N.nodf,
    "deg_het_row"        = ~deg_het_row,
    "deg_het_col"        = ~deg_het_col
)

if (is.na(MAX_SIG)) {
    acceptable_metrics <- c("CM_Ratio", "ER_Ratio", "MP_Ratio")
} else if (is.infinite(MAX_SIG)) {
    acceptable_metrics <- names(available_metrics)
} else {
    ## Check that the metrics are not dependent on matrix size/connectance
    zzz <- results_df %>%
        filter(randomization == "None") %>%
        mutate_(.dots=available_metrics) %>%
        select(type, nrows, ncols, nlinks, names(available_metrics)) %>%
        ## dependent variables
        gather("Variable1", "Value1", -type, -nrows, -ncols, -nlinks) %>%
        ## independent variables
        gather("Variable2", "Value2", -type, -Variable1, -Value1) %>%
        group_by(type, Variable1, Variable2) %>%
        filter(is.finite(Value1)) %>%
        do(lm(Value2~Value1, data=.) %>%
               tidy() %>%
               filter(term == "Value1") %>%
               mutate(direction = sign(estimate)) %>%
               select(direction, p.value)) %>%
        ungroup() %>%
        ## extract the significant linear relations (accounting for multiple testing)
        mutate(n_tests = nrow(.),
               adjusted_p.value = p.adjust(p.value, method="bonferroni")) %>%
        filter(adjusted_p.value < 0.05)
    
    # zzz %>% count(Variable1, wt=direction) %>% as.data.frame() %>% arrange(n)
    zzz %>% group_by(Variable2) %>% count(Variable1, wt=direction) %>%
        spread(Variable2, n, fill=0) %>% mutate(sum=ncols + nlinks + nrows,
                                                total = abs(ncols) + abs(nlinks) + abs(nrows)) %>%
        as.data.frame() %>% arrange(total) %>% print()
    
    # results_df %>%
    #     # mutate_(.dots=available_metrics) %>%
    #     filter(randomization == "None") %>%
    #     mutate(yval = cent_between * log(nlinks) * log(nrows + ncols)) %>%
    #     gather("xvar", "xval", ncols, nrows, nlinks) %>%
    #     group_by(randomization, type, xvar) %>%
    #     do(lm(xval~yval, data=.) %>% tidy() %>% filter(term == "yval")) %>%
    #     ungroup() %>%
    #     mutate(direction = sign(estimate),
    #            sig=ifelse(p.value < 0.001, "***",
    #                       ifelse(p.value < 0.01, "**",
    #                              ifelse(p.value < 0.05, "*", "")))) %>%
    #     mutate(dirsig = ifelse(sig == "", 0, direction)) %>%
    #     select(xvar, type, dirsig) %>%
    #     spread(xvar, dirsig) %>%
    #     summarise_at(.vars=vars(-type), .funs=list(sum))
    
    acceptable_metrics <- zzz %>%
        group_by(Variable2) %>%
        count(Variable1, wt=direction) %>%
        spread(Variable2, n) %>%
        mutate(sum=ncols + nlinks + nrows,
               total = abs(ncols) + abs(nlinks) + abs(nrows)) %>%
        filter(total <= MAX_SIG) %>%
        .$Variable1
    Metrics <- available_metrics[acceptable_metrics]
}
Metrics <- available_metrics[acceptable_metrics]

#### subsetting to equal sample sizes for fitting
if ("" == wemp) {
    to_fit_df <- results_df %>%
        filter(randomization == "None",
               type != "biogeography",
               type != "antagonism",
               type != "mutualism") %>%
        group_by(type) %>%
        sample_n(N_FIT) %>% 
        ungroup()
} else {
    to_fit_df <- results_df %>%
        filter(randomization == "None",
               type != "biogeography") %>%
        group_by(type) %>%
        sample_n(N_FIT) %>% 
        ungroup()
}

#### run the PCA
pca_df <- to_fit_df %>% transmute_(.dots=Metrics)
pca_results <- prcomp(pca_df, center=TRUE, scale=TRUE)

#### TABLE OF PCA LOADINGS ####
pca_results$rotation %>% rbind(pca_results$sdev^2 / sum(pca_results$sdev^2) * 100) %>% xtable(digits=3) %>% print()

#### fitting vs testing data plot ####
## plot the fit data
g_fit <- ggbiplot(pca_results, groups=to_fit_df$type, subgroups=to_fit_df$type, ellipse=TRUE, var.axes=FALSE, alpha=0.35) +
    scale_colour_manual(name="Network Type", values=cols) + scale_shape_manual(name="Network Type", values=shps) +
    coord_equal()
## add facetting variable to data
g_fit$data$fitvtest <- "Fitting Data"

## plot the test data 
## get the coordinates
if ("" == wemp) {
    to_test_df <- results_df %>%
        filter(randomization == "None",
               type != "biogeography",
               type != "antagonism",
               type != "mutualism")
} else {
    to_test_df <- results_df %>%
        filter(randomization == "None",
               type != "biogeography")
}

to_test_df <- to_test_df %>%
    anti_join(to_fit_df) %>%
    mutate_(.dots=Metrics) %>%
    ## center and scale according to the fitting data
    mutate_at(names(Metrics),
              funs((. - (pca_df %>% summarise_all(mean))$.) /
                       (pca_df %>% summarise_all(sd))$.))
## transform into PCA space
to_test_df <- to_test_df %>% 
    rowwise() %>%
    ## the ggbiplot code divides the coordinates by the standard deviation
    ## in the pca results. I don't know why, but to make the points line up
    ## with the ellipses, we do so here as well.
    do(as.numeric(.[names(Metrics)]) %*% pca_results$rotation %>%
           sweep(2, pca_results$sdev, '/') %>% 
           tbl_df()) %>% 
    ungroup() %>%
    bind_cols(to_test_df) %>%
    mutate(fitvtest = "Testing Data")

## remove the points from the previous plot before adding the new ones
g_fitvtest <- g_fit  +
    geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.35, data=to_test_df) +
    facet_wrap(~fitvtest, nrow=1) +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          strip.text=element_text(size=14),
          legend.text=element_text(size=12))
if (is.na(MAX_SIG)) {
    ggsave(g_fitvtest, width=1.5*BIPLOT_WIDTH, height=0.75*BIPLOT_HEIGHT + 0.5,
           filename=str_c("../Figures/FitvTest_Main_", wemp,
                          str_replace(basename(results_file), "\\.csv", ".pdf")))
} else {
    ggsave(g_fitvtest, width=2*BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
           filename=str_c("../Figures/FitvTest_", wemp, str_c("lt", MAX_SIG), "_",
                          str_replace(basename(results_file), "\\.csv", ".pdf")))
}

#### plot the randomized webs ####
## get the coordinates
if ("" == wemp) {
    full_df <- results_df %>%
        filter(type != "biogeography",
               type != "antagonism",
               type != "mutualism")
} else {
    full_df <- results_df %>% filter(type != "biogeography")
    
}
full_df <- full_df %>%
    mutate_(.dots=Metrics) %>%
    ## center and scale according to the fitting data
    mutate_at(names(Metrics),
              funs((. - (pca_df %>% summarise_all(mean))$.) /
                       (pca_df %>% summarise_all(sd))$.))
## transform into PCA space
full_df <- full_df %>% 
    rowwise() %>%
    ## the ggbiplot code divides the coordinates by the standard deviation
    ## in the pca results. I don't know why, but to make the points line up
    ## with the ellipses, we do so here as well.
    do(as.numeric(.[names(Metrics)]) %*% pca_results$rotation %>%
           sweep(2, pca_results$sdev, '/') %>% 
           tbl_df()) %>%
    ungroup() %>%
    bind_cols(full_df)
## order for nicer plotting
full_df$randomization <- factor(full_df$randomization,
                                levels=c("None", "Erdos-Renyi", "Configuration model"),
                                labels=c("Empirical", "Erdős-Rényi", "Configuration"))
## remove the fitting data points before adding the new ones
g_rand <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.35, data=full_df) +
    facet_wrap(~randomization, ncol=1) +
    theme(legend.position="bottom",
          legend.title=element_blank())
if (is.na(MAX_SIG)) {
    cairo_pdf(filename=str_c("../Figures/RandComp_Main_", wemp, 
                             str_replace(basename(results_file), "\\.csv", ".pdf")),
              width=BIPLOT_WIDTH, height=2.55*BIPLOT_HEIGHT + 0.5)
    plot(g_rand)
    dev.off()
} else {
    cairo_pdf(filename=str_c("../Figures/RandComp_", wemp, str_c("lt", MAX_SIG), "_",
                             str_replace(basename(results_file), "\\.csv", ".pdf")),
              width=BIPLOT_WIDTH, height=2.55*BIPLOT_HEIGHT + 0.5)
    plot(g_rand)
    dev.off()
}

#### empirical overlay ####
emp_df <- results_df %>%
    filter(type == "antagonism" | type == "mutualism") %>% #| type == "biogeography") %>%
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>%
    filter(randomization == "None") %>% 
    ## standardize feature1 between ecological networks
    mutate(feature1 = tolower(feature1),
           feature1 = ifelse(feature1 %in% c("antagonism", "mutualism"),
                             tolower(feature2), feature1),
           feature1 = str_replace_all(feature1, c("host-parasite" = "parasitism",
                                                  "plant-herbivore" = "herbivory",
                                                  "anemonefish" = "anemone-fish",
                                                  "seeddispersal" = "seed dispersal"))) %>%
    mutate(randomization = "Empirical") %>%
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

## remove the fitting data points before adding the new ones
emp_legend <- (ggplot(emp_df) +
                   geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), size=2, data=emp_df) +
                   geom_path(aes(x=PC1, y=PC2, color=type, group=type), data=emp_df %>%
                                 group_by(type) %>%
                                 do(ellipse_path(.)) %>% 
                                 ungroup()) +
                   scale_shape_manual(values=c("antagonism"    = 4,
                                               "mutualism"     = 3)) +
                   scale_colour_manual(values=c("antagonism"   = "#e07699",
                                                "mutualism"    = "#95a7d5")) +
                   theme(legend.position="bottom",
                         legend.title=element_blank())) %>%
    extract_legend()
g_emp <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.75, size=2, data=emp_df %>% arrange(desc(type))) +
    geom_path(aes(x=PC1, y=PC2, color=type, group=type), data=emp_df %>%
                  group_by(type) %>%
                  do(ellipse_path(.)) %>%
                  ungroup()) +
    scale_colour_manual(values=c("crimes"                  = "grey50",
                                 "antagonism"              = "#e07699",
                                 "mutualism"               = "#95a7d5",
                                 "microbiome"              = "grey50",
                                 "actor collaboration"     = "grey50",
                                 "biogeography"            = "#79c8a5",
                                 "legislature"             = "grey50",
                                 "authorship"              = "grey50"), guide=FALSE) +
    # coord_cartesian(xlim=c(-1.25,0.75), ylim=c(-1.25,0.75), expand=c(0,0)) +
    theme(legend.position="bottom",
          legend.title=element_blank())

if ("" == wemp) {
    if (is.na(MAX_SIG)) {
        ggsave(g_emp %>% replace_legend(emp_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/EmpOverlay_Main_", wemp,
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    } else {
        ggsave(g_emp %>% replace_legend(emp_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/EmpOverlay_", wemp, str_c("lt", MAX_SIG), "_",
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    }
}

#### PCA table of biplots ####
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
# emp_df <- results_df %>%
#     # bind_rows(read_csv("../Results/SimpleDemo_old_results.csv")) %>%
#     filter(type == "mutualism" | type == "antagonism") %>%
#     filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>%
#     filter(randomization == "None") %>% 
#     mutate(randomization = "Empirical") %>%
#     mutate_(.dots=Metrics) %>%
#     mutate_at(names(Metrics),
#               funs((. - (pca_df %>% summarise_all(mean))$.) /
#                        (pca_df %>% summarise_all(sd))$.))
# ## transform into PCA space
# emp_df <- emp_df %>%
#     rowwise() %>%
#     do(as.numeric(.[names(Metrics)]) %*% pca_results$rotation %>%
#            sweep(2, pca_results$sdev, '/') %>% 
#            tbl_df()) %>%
#     ungroup() %>% 
#     bind_cols(emp_df)
emp_pca_table <- pca_table_plot(emp_df %>% select(type, matches("PC\\d+")))
if (is.na(MAX_SIG)) {
    ggsave(emp_pca_table %>% replace_legend(emp_legend), width=2.5*BIPLOT_WIDTH, height=2.5*BIPLOT_HEIGHT + 0.5,
           filename=str_c("../Figures/Emp_PCA_Table_Main_", wemp,
                          str_replace(basename(results_file), "\\.csv", ".pdf")))
} else {
    ggsave(emp_pca_table %>% replace_legend(emp_legend), width=2.5*BIPLOT_WIDTH, height=2.5*BIPLOT_HEIGHT + 0.5,
           filename=str_c("../Figures/Emp_PCA_Table_", str_c("lt", MAX_SIG), "_",
                          str_replace(basename(results_file), "\\.csv", ".pdf")))
}


#### crime focal plot ####
if ("" == wemp) {
    city_legend <- (ggplot(full_df) +
                        geom_point(aes(x=PC1, y=PC2, colour=feature1, shape=feature1),
                                   data=full_df %>% filter(randomization == "Empirical", type == "crimes")) +
                        geom_path(aes(x=PC1, y=PC2, color=feature1, group=feature1),
                                  data=full_df %>%
                                      filter(randomization == "Empirical", type == "crimes") %>%
                                      group_by(feature1) %>%
                                      do(., ellipse_path(.))) +
                        facet_wrap(~randomization, ncol=1) +
                        scale_shape_manual(name="City",
                                           values=c("Chicago"                 = 0,
                                                    "Denver"                  = 1,
                                                    "Minneapolis"             = 2,
                                                    "San Francisco"           = 5,
                                                    "Washington DC"           = 6)) +
                        scale_colour_manual(name="City",
                                            values=c("Chicago"                = "#1b9e77",
                                                     "Denver"                 = "#7570b3",
                                                     "Minneapolis"            = "#e7298a",
                                                     "San Francisco"          = "#d95f02",
                                                     "Washington DC"          = "#66a61e")) +
                        theme(legend.position="bottom",
                              legend.title=element_blank())) %>%
        extract_legend()
    
    g_crimes <- g_fit %>% remove_geom("GeomPoint") +
        geom_point(aes(x=PC1, y=PC2, colour=feature1, shape=feature1), alpha=0.35,
                   data=full_df %>% filter(randomization == "Empirical", type == "crimes")) +
        geom_path(aes(x=PC1, y=PC2, color=feature1, group=feature1), data=full_df %>%
                      filter(randomization == "Empirical", type == "crimes") %>%
                      group_by(feature1) %>%
                      do(ellipse_path(.)) %>% 
                      ungroup()) +
        facet_wrap(~randomization, ncol=1) +
        scale_colour_manual(values=c("crimes"                  = "grey50",
                                     "microbiome"              = "grey50",
                                     "actor collaboration"     = "grey50",
                                     "antagonism"              = "grey50",
                                     "mutualism"               = "grey50",
                                     "legislature"             = "grey50",
                                     "authorship"              = "grey50",
                                     "Chicago"                 = "#1b9e77",
                                     "Denver"                  = "#7570b3",
                                     "Minneapolis"             = "#e7298a",
                                     "San Francisco"           = "#d95f02",
                                     "Washington DC"           = "#66a61e")) +
        scale_shape_manual(values=c("Chicago"                 = 0,
                                    "Denver"                  = 1,
                                    "Minneapolis"             = 2,
                                    "San Francisco"           = 5,
                                    "Washington DC"           = 6)) +
        # coord_cartesian(ylim=c(-0.75,0.5), xlim=c(0.85,1.75), expand=c(0,0)) +
        theme(legend.position="bottom",
              legend.title=element_blank())
    
    if (is.na(MAX_SIG)) {
        ggsave(g_crimes %>% replace_legend(city_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/CrimeFocus_sub_Main_",
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    } else {
        ggsave(g_crimes %>% replace_legend(city_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/CrimeFocus_", str_c("lt", MAX_SIG), "_",
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    }
}

# g_fit %>% remove_geom("GeomPoint") +
#     geom_point(aes(x=PC1, y=PC2,
#                    colour=str_c(feature2, ifelse(is.na(feature3), "na", feature3)),
#                    shape=str_c(feature2, ifelse(is.na(feature3), "na", feature3))),
#                alpha=0.75, size=2, 
#                data=full_df %>% filter(randomization == "Empirical", type == "legislature")) +
#     geom_path(aes(x=PC1, y=PC2,
#                   color=str_c(feature2, ifelse(is.na(feature3), "na", feature3)),
#                   group=str_c(feature2, ifelse(is.na(feature3), "na", feature3))),
#               data=full_df %>%
#                   filter(randomization == "Empirical", type == "legislature") %>%
#                   group_by(feature2, feature3) %>%
#                   do(ellipse_path(.)) %>% 
#                   ungroup()) +
#     scale_colour_discrete() + scale_shape_discrete() +
#     theme(legend.position="bottom",
#           legend.title=element_blank())


# #### Movie Genre subgroups
# g_fit %>% remove_geom(c("GeomPoint", "GeomPath")) +
#     geom_point(aes(x=PC1, y=PC2, colour=feature1), alpha=0.35,
#                data=full_df %>% filter(randomization == "Empirical", type == "actor collaboration")) +
#     geom_path(aes(x=PC1, y=PC2, color=feature1, group=feature1), data=full_df %>%
#                   filter(randomization == "Empirical", type == "actor collaboration") %>%
#                   group_by(feature1) %>%
#                   do(ellipse_path(.)) %>% 
#                   ungroup()) +
#     facet_wrap(~randomization, ncol=1) + scale_color_discrete()


#### Ecological subgroups
if ("" == wemp) {
    sub_emp_legend <- (ggplot(emp_df) +
                           geom_point(aes(x=PC1, y=PC2, colour=feature1, shape=feature1), size=2, data=emp_df) +
                           geom_path(aes(x=PC1, y=PC2, color=feature1, group=feature1), data=emp_df %>%
                                         filter(feature1 != "anemone-fish") %>% 
                                         group_by(feature1) %>%
                                         do(ellipse_path(.)) %>% 
                                         ungroup()) +
                           scale_shape_manual(values=c("anemone-fish"    = 0,
                                                       "ant-plant"       = 1,
                                                       "pollination"     = 2,
                                                       "seed dispersal"  = 3,
                                                       "bacteria-phage"  = 4,
                                                       "herbivory"       = 5,
                                                       "host-parasitoid" = 6,
                                                       "parasitism"      = 7)) +
                           scale_colour_manual(values=c("anemone-fish"    = "#666666",
                                                        "ant-plant"       = "#7570b3",
                                                        "pollination"     = "#79c8a5",
                                                        "seed dispersal"  = "#1b9e77",
                                                        "bacteria-phage"  = "#e07699",
                                                        "herbivory"       = "#95a7d5",
                                                        "host-parasitoid" = "#a6761d",
                                                        "parasitism"      = "#d95f02")) +
                           theme(legend.position="bottom",
                                 legend.title=element_blank())) %>%
        extract_legend()
    g_sub_emp <- g_fit %>% remove_geom("GeomPoint") +
        geom_point(aes(x=PC1, y=PC2, colour=feature1, shape=feature1), alpha=0.5, size=2, data=emp_df) +
        geom_path(aes(x=PC1, y=PC2, color=feature1, group=feature1), data=emp_df %>%
                      filter(feature1 != "anemone-fish") %>%
                      group_by(feature1) %>%
                      do(ellipse_path(.)) %>%
                      ungroup()) +
        scale_colour_manual(values=c("crimes"              = "grey50",
                                     "microbiome"          = "grey50",
                                     "actor collaboration" = "grey50",
                                     "legislature"         = "grey50",
                                     "authorship"          = "grey50",
                                     "antagonism"          = "grey50",
                                     "mutualism"           = "grey50",
                                     "anemone-fish"        = "#666666",
                                     "ant-plant"           = "#7570b3",
                                     "pollination"         = "#79c8a5",
                                     "seed dispersal"      = "#1b9e77",
                                     "bacteria-phage"      = "#e07699",
                                     "herbivory"           = "#95a7d5",
                                     "host-parasitoid"     = "#a6761d",
                                     "parasitism"          = "#d95f02")) +
        scale_shape_manual(values=c("anemone-fish"        = 0,
                                    "ant-plant"           = 1,
                                    "pollination"         = 2,
                                    "seed dispersal"      = 3,
                                    "bacteria-phage"      = 4,
                                    "herbivory"           = 5,
                                    "host-parasitoid"     = 6,
                                    "parasitism"          = 7)) +
        theme(legend.position="bottom",
              legend.title=element_blank())
    
    if (is.na(MAX_SIG)) {
        ggsave(g_sub_emp %>% replace_legend(sub_emp_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/EmpOverlay_Subtypes_Main_",
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    } else {
        ggsave(g_sub_emp %>% replace_legend(sub_emp_legend), width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
               filename=str_c("../Figures/EmpOverlay_Subtypes_", str_c("lt", MAX_SIG), "_",
                              str_replace(basename(results_file), "\\.csv", ".pdf")))
    }
}

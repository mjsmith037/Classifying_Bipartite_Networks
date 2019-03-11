library(scales)
library(stringr)
library(magrittr)
library(tidyverse)

source("ggbiplot.R")
theme_set(theme_bw())

source("ancillary_functions.R")

BIPLOT_HEIGHT <- 6
BIPLOT_WIDTH <- 6

set.seed(0) # to replicate sampling

#### KEY PARAMETERS
N_FIT <- 100    # size of subsets used for fitting
MIN_SIZE <- 5   # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)
wemp <- FALSE   # Should ecological data be included in the fitting?

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
results_df <- read_csv("../../Results/SimpleDemo_results.csv", col_types="cccddddddd") %>%
  # add in the metadata
  left_join(read_csv("../../Data/Metadata.csv", col_types="cccccccccccccidii"),
            by=c("type", "name")) %>%
  ## clean up some names
  mutate(type = str_replace_all(type, c("actorcollaboration" = "actor collaboration"))) %>%
  mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature_1), type)) %>%
  ## remove very small networks
  filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE)

#### specify/calculate the variables to consider in PCA
Metrics <- quos(
  "CM_Ratio"           = l1_cm / l1,
  "ER_Ratio"           = l1_er / l1,
  "MP_Ratio"           = ((1 + sqrt(ncols / nrows)) *
                             sqrt(nlinks / ncols *
                                    (1 - nlinks / nrows / ncols))) / l2
)

#### subsetting to equal sample sizes for fitting
if (wemp) {
  to_fit_df <- results_df %>%
    filter(randomization == "empirical",
           type != "biogeography") %>% # always omit biogeography because too few networks
    group_by(type) %>%
    sample_n(N_FIT) %>%
    ungroup()
} else {
  to_fit_df <- results_df %>%
    filter(randomization == "empirical",
           type != "biogeography", # always omit biogeography because too few networks
           type != "antagonism",
           type != "mutualism") %>%
    group_by(type) %>%
    sample_n(N_FIT) %>%
    ungroup()
}

#### run the PCA
pca_df <- to_fit_df %>% transmute(!!! Metrics)
pca_results <- prcomp(pca_df, center=TRUE, scale=TRUE)

#### TABLE OF PCA LOADINGS ####
pca_results$rotation %>% rbind(pca_results$sdev^2 / sum(pca_results$sdev^2) * 100) %>% print()

#### fitting vs testing data plot ####
## plot the fit data
g_fit <- ggbiplot(pca_results, groups=to_fit_df$type, subgroups=to_fit_df$type, ellipse=TRUE, var.axes=FALSE, alpha=0.35) +
  scale_colour_manual(name="Network Type", values=cols) + scale_shape_manual(name="Network Type", values=shps) +
  coord_equal()

#### plot the randomized webs ####
## get the coordinates
if (wemp) {
  full_df <- results_df %>% filter(type != "biogeography")
} else {
  full_df <- results_df %>%
    filter(type != "biogeography",
           type != "antagonism",
           type != "mutualism")
}
full_df <- full_df %>%
  ## center and scale according to the fitting data
  bind_cols(transmute(., !!! Metrics) %>%
              subtract(pca_df %>% summarise_all(mean) %>% uncount(nrow(full_df))) %>%
              divide_by(pca_df %>% summarise_all(sd) %>% uncount(nrow(full_df))) %>%
              ## transform into PCA space
              rowwise() %>%
              do(as.numeric(.) %>%
                   multiply_by_matrix(pca_results$rotation) %>%
                   ## the ggbiplot code divides the coordinates by the standard deviation
                   ## in the pca results. I don't know why, but to make the points line up
                   ## with the ellipses, we do so here as well.
                   divide_by(pca_results$sdev) %>%
                   tbl_df()) %>%
              ungroup())
## order for nicer plotting
full_df$randomization <- factor(full_df$randomization,
                                levels=c("empirical", "erdos-renyi", "configuration"))
## remove the fitting data points before adding the new ones
g_rand <- g_fit %>% remove_geom("GeomPoint") +
  geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.35, data=full_df) +
  facet_wrap(~randomization, ncol=1) +
  theme(legend.position="bottom",
        legend.title=element_blank())
ggsave(filename="../../Figures/pca_simple_demonstration.pdf",
       width=BIPLOT_WIDTH, height=2.55*BIPLOT_HEIGHT + 0.5)

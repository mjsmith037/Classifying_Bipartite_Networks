library(scales)
library(stringr)
library(magrittr)
library(tidyverse)

source("ggbiplot.R")
theme_set(theme_bw())

source("ancillary_functions.R")

BIPLOT_HEIGHT <- 6
BIPLOT_WIDTH <- 6

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

set.seed(0) # to replicate sampling

#### KEY PARAMETERS ####
N_FIT <- 100 # size of subsets used for fitting
MIN_SIZE <- 5 # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)
wEco <- TRUE

#### Read in the data ####
# specifying column types (while not neccessary) speeds reading in of results files
# specify these columns to be doubles to avoid errors in reading very large integers
column_specs <- cols(
  randomization = col_character(),
  name = col_character(),
  type = col_character(),
  l1 = col_double(), l2 = col_double(), l3 = col_double(),
  ipr1 = col_double(), ipr2 = col_double(), ipr3 = col_double(),
  l1_cm = col_double(), l1_er = col_double(), l1_reg = col_double(),
  l2_mp = col_double(), l3_mp = col_double(),
  l1_lap = col_double(), l2_lap = col_double(), l3_lap = col_double(),
  alg_conn = col_double(),
  clustering_c = col_double(), clustering_r = col_double(),
  deg_assort = col_double(),
  Q = col_double(), N.olap = col_double(), N.temp = col_double(), N.nodf = col_double(),
  H2 = col_double(), H3 = col_double(), H4 = col_double(), H17 = col_double(),
  cent_between = col_double(), cent_close = col_double(), cent_eigen = col_double(),
  diam = col_double(), mean_path_length = col_double(),
  deg_het_row = col_double(), deg_het_col = col_double()
)

# compile all individual results into a single table of metrics (takes a little while)
full_results <- Sys.glob("../../Results/*/*/*.csv") %>%
  lapply(read_csv, progress=FALSE, col_types=column_specs) %>%
  bind_rows() %>%
  # add in the metadata
  left_join(read_csv("../../Data/Metadata.csv", col_types="cccccccccccccidii"), by=c("type", "name")) %>%
  ## clean up some names
  mutate(type = str_replace_all(type, c("actorcollaboration" = "actor collaboration"))) %>%
  mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature_1), type)) %>%
  ## remove very small networks
  filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE)

#### specify/calculate the variables to consider in PCA ####
Metrics <- quos(
  CM_Ratio           = l1_cm / l1,
  ER_Ratio           = l1_er / l1,
  MP_Ratio           = ((1 + sqrt(ncols / nrows)) *
                          sqrt(nlinks / ncols * (1 - connectance))) / l2,
  Reg_Ratio          = l1 / l1_reg,
  Lap_Ratio          = l1 / l1_lap,
  Gap_Ratio_21       = (l1 - l2) / l1,
  Gap_Ratio_32       = (l2 - l3) / l2,
  Lap_Gap_Ratio_21   = (l1_lap - l2_lap) / l1_lap,
  Lap_Gap_Ratio_32   = (l2_lap - l3_lap) / l2_lap,
  Gap_Lap_Gap_Ratio  = (l1_lap - l2_lap) / (l1 - l2),
  l2l3_Ratio         = l3 / l2,
  ipr1               = ipr1,
  ipr2               = ipr2,
  ipr3               = ipr3,
  Alg_Connectance    = alg_conn,
  cent_between       = cent_between,
  cent_close         = cent_close,
  cent_eigen         = cent_eigen,
  av_short_path      = mean_path_length,
  diam               = diam,
  log_H2             = log(H2),
  log_H3             = log(H3),
  log_H4             = log(H4),
  log_H17            = log(H17 + 1), # +1 because not all webs have an H17 subgraph
  H4H2               = H4 / H2,
  H17H4              = H17 / H4,
  deg_het_row        = deg_het_row,
  deg_het_col        = deg_het_col,
  deg_assort         = deg_assort
  ## These metrics are not present for a substatntial (> 10%) fraction of networks
  # Q                  = Q,
  # N.olap             = N.olap,
  # N.temp             = N.temp,
  # N.nodf             = N.nodf,
  # clustering_r       = clustering_r / (1 - (1 - (connectance)^2)^(ncols - 2)),
  # clustering_c       = clustering_c / (1 - (1 - (connectance)^2)^(nrows - 2)),
)

#### subsetting to equal sample sizes for fitting ####
if (wEco) {
  to_fit_df <- full_results %>%
    filter(randomization == "empirical",
           type != "biogeography", # always omit biogeography because few networks
           type != "antagonism",
           type != "mutualism") %>%
    group_by(type) %>%
    sample_n(N_FIT) %>%
    ungroup()
} else {
  to_fit_df <- full_results %>%
    filter(randomization == "empirical",
           type != "biogeography") %>% # always omit biogeography because few networks
    group_by(type) %>%
    sample_n(N_FIT) %>%
    ungroup()
}

#### run the PCA ####
pca_df <- to_fit_df %>% transmute(!!! Metrics)
pca_results <- prcomp(pca_df, center=TRUE, scale=TRUE)

#### Plotting ####
# plot the fit data
g_fit <- ggbiplot(pca_results, groups=to_fit_df$type, subgroups=to_fit_df$type,
                  ellipse=TRUE, var.axes=FALSE, alpha=0.35) +
  scale_colour_manual(name="Network Type", values=cols) +
  scale_shape_manual(name="Network Type", values=shps) +
  coord_equal()
# add facetting variable to data
g_fit$data$fitvtest <- "Fitting Data"

# plot the test data
if (wEco) {
  to_test_df <- full_results %>%
    filter(randomization == "empirical",
           type != "biogeography", # always omit biogeography because few networks
           type != "antagonism",
           type != "mutualism")
} else {
  to_test_df <- full_results %>%
    filter(randomization == "empirical",
           type != "biogeography") # always omit biogeography because few networks
}

to_test_df <- to_test_df %>%
  anti_join(to_fit_df) %>%
  bind_cols(transmute(., !!! Metrics) %>%
              # center and scale according to the fitting data
              subtract(pca_df %>% summarise_all(mean) %>% uncount(nrow(to_test_df) - nrow(to_fit_df))) %>%
              divide_by(pca_df %>% summarise_all(sd) %>% uncount(nrow(to_test_df) - nrow(to_fit_df))) %>%
              # transform into PCA space
              rowwise() %>%
              do(as.numeric(.) %>%
                   multiply_by_matrix(pca_results$rotation) %>%
                   # the ggbiplot code divides the coordinates by the standard deviation
                   # in the pca results. I don't know why, but to make the points line up
                   # with the ellipses, we do so here as well.
                   divide_by(pca_results$sdev) %>%
                   tbl_df()) %>%
              ungroup()) %>%
  mutate(fitvtest = "Testing Data")

# remove the points from the previous plot before adding the new ones
g_fitvtest <- g_fit  +
  geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.35, data=to_test_df) +
  facet_wrap(~fitvtest, nrow=1) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        strip.text=element_text(size=14),
        legend.text=element_text(size=12))
# ggsave(g_fitvtest, width=1.5*BIPLOT_WIDTH, height=0.75*BIPLOT_HEIGHT + 0.5,
#        filename=str_c("../../Figures/FitvTest_",
#                       ifelse(wEco, "FitwithEcoData", "FitwithoutEcoData"), ".pdf"))

# get the principal coordinates for all of the data
full_df <- full_results %>%
  # order for nicer plotting
  mutate(randomization = factor(randomization,
                                levels=c("empirical", "erdos-renyi", "configuration"))) %>%
  # center and scale according to the fitting data
  bind_cols(transmute(., !!! Metrics) %>%
              subtract(pca_df %>% summarise_all(mean) %>% uncount(nrow(full_results))) %>%
              divide_by(pca_df %>% summarise_all(sd) %>% uncount(nrow(full_results))) %>%
              # transform into PCA space
              rowwise() %>%
              do(as.numeric(.) %>%
                   multiply_by_matrix(pca_results$rotation) %>%
                   # the ggbiplot code divides the coordinates by the standard deviation
                   # in the pca results. I don't know why, but to make the points line up
                   # with the ellipses, we do so here as well.
                   divide_by(pca_results$sdev) %>%
                   tbl_df()) %>%
              ungroup())

## Comparison to randomizations (Figure 1)
source("randomization_comparison_figure.R")
rand_comp_plot <- get_randomization_comparison_plot(g_fit, wEco, full_df)
ggsave(rand_comp_plot, width=BIPLOT_WIDTH, height=2.55*BIPLOT_HEIGHT + 0.5,
       filename=str_c("../../Figures/RandomizationComparison_",
                      ifelse(wEco, "FitwithEcoData", "FitwithoutEcoData"), ".pdf"))

## Crime sub-structure (Figure 2)
source("crime_substructure_figure.R")
crime_substructure_plot <- get_crime_substructure_plot(g_fit, wEco, full_df)
ggsave(crime_substructure_plot, width=BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
       filename=str_c("../../Figures/CrimeSubstructure_",
                      ifelse(wEco, "FitwithEcoData", "FitwithoutEcoData"), ".pdf"))

## Ecological data overlay (Figure 3)
# Note that this plot only makes sense if ecological data was not used in the fitting
source("ecological_overlay_figure.R")
eco_overlay_plot <- get_ecological_overlay_plot(g_fit, full_df)
ggsave(eco_overlay_plot, width=1.75*BIPLOT_WIDTH, height=BIPLOT_HEIGHT + 0.5,
       filename="../../Figures/EcologicalOverlay_FitwithoutEcoData.pdf")


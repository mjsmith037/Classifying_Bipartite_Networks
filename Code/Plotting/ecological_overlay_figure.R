#### This script creates a plot where ecological data is overlaid into
#### an established principal component space
require(tidyverse)
require(cowplot)
theme_set(theme_bw()) # have to reset theme to undo cowplot's retheming

source("ancillary_functions.R")

eco_shapes <- c("antagonism"    = 4,
                "mutualism"     = 3)
eco_colours <- c("antagonism"   = "#e07699",
                 "mutualism"    = "#95a7d5")

sub_eco_shapes <- c("anemone-fish"    = 0,
                    "ant-plant"       = 1,
                    "pollination"     = 2,
                    "seed dispersal"  = 3,
                    "bacteria-phage"  = 4,
                    "herbivory"       = 5,
                    "host-parasitoid" = 6,
                    "parasitism"      = 7)
sub_eco_colours <- c("anemone-fish"    = "#666666",
                     "ant-plant"       = "#7570b3",
                     "pollination"     = "#79c8a5",
                     "seed dispersal"  = "#1b9e77",
                     "bacteria-phage"  = "#e07699",
                     "herbivory"       = "#95a7d5",
                     "host-parasitoid" = "#a6761d",
                     "parasitism"      = "#d95f02")

## remove the fitting data points before adding the new ones
get_ecological_legend <- function(eco_data) {
  extract_legend(ggplot(eco_data) +
                   geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), size=2) +
                   scale_shape_manual(values=eco_shapes) +
                   scale_colour_manual(values=eco_colours) +
                   theme(legend.position="bottom",
                         legend.title=element_blank()))
}

get_subecological_legend <- function(eco_data) {
  extract_legend(ggplot(eco_data) +
                   geom_point(aes(x=PC1, y=PC2, colour=feature_2, shape=feature_2), size=2) +
                   scale_shape_manual(values=sub_eco_shapes) +
                   scale_colour_manual(values=sub_eco_colours) +
                   theme(legend.position="bottom",
                         legend.title=element_blank()))
}

get_ecological_overlay_plot <- function(g_fit, full_df) {
  eco_data <- full_df %>%
    filter(type == "antagonism" | type == "mutualism") %>%
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>%
    filter(randomization == "empirical",
           feature_1 %in% c("antagonism", "mutualism"))

  g_emp <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.75, size=2, data=eco_data) +
    geom_path(aes(x=PC1, y=PC2, color=type, group=type), data=eco_data %>%
                group_by(type) %>%
                do(ellipse_path(.)) %>%
                ungroup()) +
    scale_colour_manual(values=c("crimes"                  = "grey50",
                                 "microbiome"              = "grey50",
                                 "actor collaboration"     = "grey50",
                                 "legislature"             = "grey50",
                                 "authorship"              = "grey50",
                                 eco_colours), guide=FALSE) +
    theme(legend.position="bottom",
          legend.title=element_blank())

  g_sub_emp <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=feature_2, shape=feature_2), alpha=0.5, size=2, data=eco_data) +
    geom_path(aes(x=PC1, y=PC2, color=feature_2, group=feature_2), data=eco_data %>%
                filter(feature_2 != "anemone-fish") %>%
                group_by(feature_2) %>%
                do(ellipse_path(.)) %>%
                ungroup()) +
    scale_colour_manual(values=c("crimes"              = "grey50",
                                 "microbiome"          = "grey50",
                                 "actor collaboration" = "grey50",
                                 "legislature"         = "grey50",
                                 "authorship"          = "grey50",
                                 "antagonism"          = "grey50",
                                 "mutualism"           = "grey50",
                                 sub_eco_colours), guide=FALSE) +
    scale_shape_manual(values=sub_eco_shapes) +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  eco_overlay_plot <- plot_grid(replace_legend(g_emp, get_ecological_legend(eco_data)),
                                replace_legend(g_sub_emp, get_subecological_legend(eco_data)),
                                align="hv", axis="lr")

  # return the plot
  return(eco_overlay_plot)
}

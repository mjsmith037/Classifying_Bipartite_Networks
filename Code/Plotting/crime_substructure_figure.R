require(stringr)
require(tidyverse)

source("ancillary_functions.R")

#### crime focal plot ####
crime_colours <- c("Chicago"                 = "#1b9e77",
                   "Denver"                  = "#7570b3",
                   "Minneapolis"             = "#e7298a",
                   "San Francisco"           = "#d95f02",
                   "Washington DC"           = "#66a61e")
crime_shapes <- c("Chicago"                 = 0,
                  "Denver"                  = 1,
                  "Minneapolis"             = 2,
                  "San Francisco"           = 5,
                  "Washington DC"           = 6)
# a function to get a legend that only contains the cities from the crime data
get_city_legend <- function(crime_data) {
  extract_legend(ggplot(crime_data) +
                   geom_point(aes(x=PC1, y=PC2, colour=feature_1, shape=feature_1)) +
                   geom_path(aes(x=PC1, y=PC2, color=feature_1, group=feature_1),
                             data=group_by(crime_data, feature_1) %>% do(ellipse_path(.))) +
                   facet_wrap(~randomization, ncol=1) +
                   scale_shape_manual(name="City", values=crime_shapes) +
                   scale_colour_manual(name="City", values=crime_colours) +
                   theme(legend.position="bottom",
                         legend.title=element_blank()))
}
# plot the crime data colored by city and with a subsetted legend
get_crime_substructure_plot <-function(g_fit, wEco, full_df) {
  # extract the relevant information
  crime_data <- full_df %>%
    mutate(feature_1 = str_replace_all(str_to_title(feature_1),
                                       "Washington Dc", "Washington DC")) %>%
    filter(randomization == "empirical", type == "crimes")
  # the base plot
  g_crimes <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=feature_1, shape=feature_1), alpha=0.35,
               data=filter(crime_data, randomization == "empirical", type == "crimes")) +
    geom_path(aes(x=PC1, y=PC2, color=feature_1, group=feature_1), data=crime_data %>%
                group_by(feature_1) %>%
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
                                 crime_colours)) +
    scale_shape_manual(values=crime_shapes) +
    theme(legend.position="bottom",
          legend.title=element_blank())
  # replace the legend
  crime_substructure_plot <- replace_legend(g_crimes, get_city_legend(crime_data))
  # return the plot
  return(crime_substructure_plot)
}

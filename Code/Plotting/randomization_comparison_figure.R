#### This function creates a plot comparing a set of empirical networks to their
#### randomizations when overlaid into an established principal component space
get_randomization_comparison_plot <- function(g_fit, wEco, full_df) {
  # remove the fitting data points before adding the new ones
  rand_comp_plot <- g_fit %>% remove_geom("GeomPoint") +
    geom_point(aes(x=PC1, y=PC2, colour=type, shape=type), alpha=0.35, data=full_df) +
    facet_wrap(~randomization, ncol=1) +
    theme(legend.position="bottom",
          legend.title=element_blank())
  # return the plot
  return(rand_comp_plot)
}

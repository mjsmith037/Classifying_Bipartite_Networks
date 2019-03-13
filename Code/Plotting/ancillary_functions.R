## This function takes in a ggplot object and returns the object with all geoms
## of a specified type(s) removed
remove_geom <- function(ggplot2_object, geom_type) {
    ## Nullify layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
        if (class(x$geom)[1] %in% geom_type) NULL else x
    })
    ## Delete the unwanted layers.
    layers <- layers[!unlist(lapply(layers, is.null))]
    ## and replace them in the original ggplot object before returning
    ggplot2_object$layers <- layers
    return(ggplot2_object)
}

## This function takes a ggplot object and replaces the legend with a provided
## alternative
replace_legend <- function(plot_obj, new_legend) {
    tmp <- ggplot_gtable(ggplot_build(plot_obj))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]] <- new_legend
    return(tmp)
}

## this function takes a ggplot object and returns the grob that contains the
## legend (for insertion into another plot)
extract_legend <- function(plot_obj) {
    tmp <- ggplot_gtable(ggplot_build(plot_obj))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

## given a dataframe and the names of two numerical columns, this function
## returns points along a path forming an ellipse containing approximately
## 68% of the points (assuming a multivariate normal distribution)
ellipse_path <- function(x, one="PC1", two="PC2", num_points=100) {
    theta <- c(seq(-pi, pi, length=num_points/2), seq(pi, -pi, length=num_points/2))
    base_circle <- cbind(cos(theta), sin(theta))
    if(nrow(x) <= 2) return(NULL)
    sigma <- x %>% select_(one, two) %>% filter_all(all_vars(is.finite(.))) %>% var()
    mu <- x %>% select(one, two) %>% filter_all(all_vars(is.finite(.))) %>%
        summarise_all("mean") %>% as.numeric()
    ed <- 0.68 %>% qchisq(df=2) %>% sqrt()
    return(data.frame(sweep(base_circle %*% chol(sigma) * ed, 2, mu, FUN = '+')))
}

## this function takes in a dataframe of pca results and cleans up the data for
## a facetted ggplot of the various biplots of PC combinations
pca_table_plot_data <- function(dat, n_pcs) {
    zzz <- dat %>%
        select(matches(str_c("PC", 1:n_pcs, "$", collapse="|"))) %>%
                   gather("Axis", "Value", everything())
    zzz %>%
        group_by(Axis) %>%
               do(lapply(unique(zzz$Axis),
                         function(xx)
                             bind_cols(., filter(zzz, Axis == xx))) %>%
                      bind_rows())
}

library(stringr)
library(tidyverse)


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
results_df <- read_csv("../Results/full_results.csv",
                       ## specify type for these columns because they vary between networks or cause errors in parsing
                       col_types=cols(nrows = col_double(),
                                      feature1 = col_character(),
                                      feature2 = col_character(),
                                      feature3 = col_character(),
                                      feature4 = col_character()))
results_df <- results_df %>%
  ## clean up some names
  mutate(type = tolower(type),
         type = str_replace(type, "movies", "actor collaboration"),
         type = str_replace(type, "ecologicalinteractions", "ecological interactions")) %>%
  mutate(type = ifelse(type == "ecological interactions", tolower(feature1), type)) %>%
  mutate(type = str_replace(type, "fungi|islands|mountains", "biogeography"),
         type = str_replace(type, "anemonefish|ant-plant|pollination|seeddispersal", "mutualism"),
         type = str_replace(type, "bacteria-phage|host-parasitoid|parasitism|plant-herbivore", "antagonism"))

toplot <- results_df %>%
    filter(randomization == "Erdos-Renyi") %>%
    mutate(p = nlinks / (nrows * ncols)) %>%
    mutate(l2_mp = (1 + sqrt(nrows / ncols)) * sqrt(ncols * p * (1-p))) %>%
    select(l1, l2, l1_cm, l1_er, l2_mp, type)


plqualityl1cm <- toplot %>% 
    ggplot(aes(l1, l1_cm, colour = type)) + 
    geom_point(alpha = 0.25) + 
    scale_x_continuous(expression(lambda[1]), trans = "log2") + 
    scale_y_continuous(expression(lambda[1]^displaystyle(cm)), trans = "log2") + 
    theme_bw() +
    scale_colour_manual(name="Network Type", values=cols) +
    scale_shape_manual(name="Network Type", values=shps)
plqualityl1cm

plqualityl1er <- toplot %>% 
    ggplot(aes(l1, l1_er, colour = type)) + 
    geom_point(alpha = 0.25) + 
    scale_x_continuous(expression(lambda[1]), trans = "log2") + 
    scale_y_continuous(expression(lambda[1]^displaystyle(er)), trans = "log2") + 
    theme_bw() +
    scale_colour_manual(name="Network Type", values=cols) +
    scale_shape_manual(name="Network Type", values=shps)
plqualityl1er

plqualityl2 <- toplot %>% 
    ggplot(aes(l2, l2_mp, colour = type)) + 
    geom_point(alpha = 0.25) + 
    scale_x_continuous(expression(lambda[2]), trans = "log2") + 
    scale_y_continuous(expression(lambda[2]^displaystyle(mp)), trans = "log2") + 
    theme_bw() +
    scale_colour_manual(name="Network Type", values=cols) +
    scale_shape_manual(name="Network Type", values=shps)
plqualityl2

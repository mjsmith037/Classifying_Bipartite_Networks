library(xtable)
library(igraph)
library(stringr)
library(tidyverse)
library(broom)

#### KEY PARAMETERS
MIN_SIZE <- 5 # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)
results_file <- "../../results/SimpleDemo_results.csv"

read_csv(results_file) %>%
  filter(randomization == "empirical") %>%
  # add in the metadata
  left_join(read_csv("../../Metadata.csv", col_types="cccccccccccccidii"), by=c("type", "name")) %>%
  ## clean up some names
  mutate(type = str_replace_all(type, c("actorcollaboration" = "actor collaboration"))) %>%
  mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature_1), type)) %>%
  mutate(type = str_to_title(type)) %>%
  ## remove very small networks
  filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>%
  select(type, nrows, ncols, nlinks, connectance, l1) %>%
  distinct(l1, .keep_all=TRUE) %>%
  left_join(count(., type), by="type") %>%
  mutate(type=str_c(type, " (", n, ")")) %>%
  select(-n, -l1) %>%
  gather("Metric", "value", -type) %>%
  group_by(type, Metric) %>%
  summarise_all(list("Minimum"=~min(., na.rm=TRUE),
                     "Maximum"=~max(., na.rm=TRUE),
                     "Mean"=~round(mean(., na.rm=TRUE), 3),
                     "Median"=~median(., na.rm=TRUE))) %>%
  mutate(Metric = factor(Metric,
                         levels=c("connectance", "nlinks", "nrows", "ncols"),
                         labels=c("Connectance", "Links", "Rows", "Cols"))) %>%
  arrange(type, Metric) %>%
  group_by(type) %>%
  mutate(Type = ifelse(duplicated(type), "",
                       paste0('\\multirow{', n(), '}*{', type, '}'))) %>%
    ungroup() %>%
    select(-type) %>%
    select(Type, everything()) %>%
    mutate_at(vars(-Type, -Metric), format, format="g", digits=1, nsmall=4,
              drop0trailing=TRUE, scientific=FALSE, big.mark=",") %>%
    xtable(caption="Summary statistics on the size and fill of collected networks.",
           label="tab:DataSummary") %>%
    print(include.rownames=FALSE,
          hline.after=c(-1,0,1:(length(unique(.$Type)) - 1) * 4),
          caption.placement="top",
          sanitize.text.function=function(xx) return(xx),
          file="../../results/tables/network_summary_stats.tex")

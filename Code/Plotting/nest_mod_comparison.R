library(stringr)
library(tidyverse)
library(broom)

theme_set(theme_bw())

cols <- c("antagonism" = "#e07699", "mutualism" = "#95a7d5")

nest_mod_df <- read_csv("../Results/nest_mod_results.csv")

MIN_SIZE <- 5
results_file <- "../Results/SimpleDemo_results.csv"
filtered_webs <- read_csv(results_file,
                       ## specify type for these columns because they vary between networks or cause errors in parsing
                       col_types=cols(nrows = col_double(),
                                      feature1 = col_character(),
                                      feature2 = col_character(),
                                      feature3 = col_character(),
                                      feature4 = col_character())) %>%
    bind_rows(read_csv("../Results/SimpleDemo_old_results.csv")) %>% 
    mutate(type = tolower(type)) %>%
    filter(type %in% c("ecologicalinteractions", "antagonism", "mutualism")) %>%
    mutate(type = ifelse(type == "ecologicalinteractions", tolower(feature1), type)) %>% 
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE, randomization == "None") %>% 
    select(file=name, type)

## normalize a vector to be between 0 and 1
normalize <- function(x) {
    y <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    return(sign(y) * log(abs(y) + 1))
}

nest_mod_df %>%
    inner_join(filtered_webs) %>%
    na.omit() %>%
    group_by(rand) %>%
    mutate_at(vars(contains("N."), Q), normalize) %>%
    ungroup() %>%
    gather("metric", "value", contains("N."), Q) %>%
    mutate(rand = factor(rand, levels=c("ER", "CM"),
                         labels=c("Erdos-Renyi", "Configuration Model"))) %>%
    ggplot() +
    aes(x=metric, y=value, fill=type) +
    geom_boxplot(aes(colour=type), fill=NA) +
    geom_boxplot(outlier.colour=NA) +
    scale_fill_manual(values=cols) + scale_colour_manual(values=cols) +
    ylab("Z-Score of Metric Value") +
    facet_wrap(~rand, ncol=1, scales="free_y") +
    theme(legend.position="bottom",
          legend.title=element_blank(),
          axis.title.x=element_blank())

ggsave(width=6, height=5, filename="../Figures/ThebaultFontaine_Main.pdf")

test_results <- nest_mod_df %>%
    inner_join(filtered_webs) %>%
    gather("metric", "value", contains("N."), Q) %>%
    group_by(rand, metric) %>%
    do(., t.test(value~type, data=.) %>% tidy()) %>%
    mutate(sig=ifelse(p.value < 0.05, "*", ""))

sample_size <- nest_mod_df %>%
    inner_join(filtered_webs) %>%
    gather("metric", "value", contains("N."), Q) %>%
    group_by(rand, metric, type) %>%
    tally() %>% 
    spread(type, n)

left_join(test_results, sample_size) %>%
    mutate(p.value       = format(p.value, digits=3),
           mu_antagonism = str_c(round(estimate, 4), " (", antagonism, ")"),
           mu_mutualism  = str_c(round(estimate1, 4), " (", mutualism, ")")) %>% 
    select("Null Model"=rand, "Metric"=metric,
            mu_antagonism, mu_mutualism,
            "degrees of freedom"=parameter, "P-Value"=p.value)
    
    

library(tidyverse)
library(stringr)

ind_res <- list.files("../../results/individual_results", recursive=TRUE, full.names=TRUE)
lapply(ind_res, function(res) read_csv(res, col_types=list("feature1" = col_character(),
                                                           "feature2" = col_character(),
                                                           "feature3" = col_character(),
                                                           "feature4" = col_character())) %>%
    mutate_at(vars(matches("l\\d")), round, digits=6)) %>% bind_rows() %>% write_csv("../../results/full_results.csv")

new_res <- ind_res[str_detect(ind_res, "/data")]
lapply(new_res, function(res) read_csv(res, col_types=list("feature1" = col_character(),
                                                           "feature2" = col_character(),
                                                           "feature3" = col_character(),
                                                           "feature4" = col_character())) %>%
    mutate_at(vars(matches("l\\d")), round, digits=6)) %>% bind_rows() %>% write_csv("../../results/only_new-data_results.csv")

old_res <- ind_res[str_detect(ind_res, "/old-data")]
lapply(old_res, function(res) read_csv(res, col_types=list("feature1" = col_character(),
                                                           "feature2" = col_character(),
                                                           "feature3" = col_character(),
                                                           "feature4" = col_character())) %>%
    mutate_at(vars(matches("l\\d")), round, digits=6)) %>% bind_rows() %>% write_csv("../../results/only_old-data_results.csv")


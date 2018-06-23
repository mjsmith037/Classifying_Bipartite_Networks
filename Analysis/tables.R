library(xtable)
library(igraph)
library(stringr)
library(tidyverse)
library(broom)

#### KEY PARAMETERS
MIN_SIZE <- 5 # minimum size of networks included in analysis (both nrows and ncols >= MIN_SIZE)
results_file <- "../Results/SimpleDemo_results.csv"

results_df <- read_csv(results_file,
                       ## specify type for these columns because they vary between networks or cause errors in parsing
                       col_types=cols(nrows = col_double(),
                                      # H2 = col_double(),
                                      # H3 = col_double(),
                                      # H4 = col_double(),
                                      feature1 = col_character(),
                                      feature2 = col_character(),
                                      feature3 = col_character(),
                                      feature4 = col_character())) %>%
    ## clean up some names
    mutate(type = tolower(type),
           type = str_replace(type, "movies", "actor collaboration"),
           type = ifelse(type == "ecologicalinteractions", tolower(feature1), type))

results_df %>%
    bind_rows(read_csv("../Results/SimpleDemo_old_results.csv")) %>%
    filter(nrows >= MIN_SIZE, ncols >= MIN_SIZE) %>% 
    filter(randomization == "None") %>%
    select(type, nrows, ncols, nlinks, connectance, l1) %>%
    mutate(type = str_to_title(type)) %>%
    distinct(l1, .keep_all=TRUE) %>% 
    group_by(type) %>% 
    mutate(n=n()) %>% 
    ungroup() %>%
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
          file="../Results/tables/network_summary_stats.tex")
    

################## ECOLOGICAL INTERACTION NETWORK REFERENCES ###################
all_files <- c(list.files("../Data/data/ecology", full.names=TRUE),
               list.files("../Data/old-data", recursive=TRUE, full.names=TRUE))
tmp <- lapply(all_files, function(ff) {
    load(ff)
    if (str_detect(ff, "SpectralDifferences")) {
        return(data_frame(BibTeX=c(as.character(metadata$Subtype), str_c(metadata$Subtype, "_data"))) %>%
                   mutate(Type=ifelse(metadata$Type == "Pollination", "mutualism", "antagonism"),
                          Subtype=metadata$Type,
                          Rows=metadata$N_Rows, Cols=metadata$N_Cols, Links=metadata$N_Links,
                          Connectance=Links / Rows / Cols))
    } else if (str_detect(ff, "BipartiteMotifs")) {
        if ("BibTeX" %in% names(Data)) BibTeX <- Data$BibTeX else BibTeX <- ""
        return(data_frame(BibTeX=BibTeX) %>%
                   mutate(Type=Data$Type, Subtype=Data$SubType,
                          Rows=Data$NRows, Cols=Data$NCols, Links=Data$NLinks,
                          Connectance=Links / Rows / Cols))
    } else {
        BibTeX <- str_c(str_extract(metadata$citation, "\\w+"), "_", str_extract(metadata$citation, "\\d\\d\\d\\d"))
        if (str_detect(ff, "Data Dryad")) BibTeX <- c(BibTeX, str_c(BibTeX, "_data"))
        return(data_frame(BibTeX=BibTeX) %>%
                   mutate(Type=metadata$feature_1, Subtype=ifelse(metadata$feature_1 == "Mutualism", "pollination", "parasitism"),
                          Rows=metadata$nrows, Cols=metadata$ncols, Links=metadata$n_links,
                          Connectance=Links / Rows / Cols))
    }
}) %>% bind_rows() %>% mutate(Type=tolower(Type), Subtype=tolower(Subtype)) %>% 
    filter(Type %in% c("antagonism", "mutualism"), Rows>5, Cols>5)

table(tmp$Subtype)

for (tt in unique(tmp$Subtype)) {
    print(tt)
    print(str_c("[", str_c("@", tmp %>% filter(Subtype == tt) %>% .$BibTeX %>% unique(), collapse=";"), "]"))
}

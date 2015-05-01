# read in the source csv file and turn it into a nice R object.
library(dplyr)
library(tidyr)


cs_uncorrelated <- tbl_df(read.csv("data-raw/CS_uncorrelated.csv",  
                                   stringsAsFactors = FALSE)) %>%
  select(-sero, -uno)

# variables are grouped into five alternatives.
attr_names <- paste(rep(letters[1:5], each = 5), 
                    rep(1:5, 5), sep = "_")
alt_names <- paste("choice", 1:5, sep = "_")

names(cs_uncorrelated) <- c(attr_names, alt_names, "id")


cs_uncorrelated <- cs_uncorrelated %>%
  gather(var, value, a_1:choice_5) %>%
  separate(var, c("var", "alt"))  %>%
  spread(var, value) %>%
  select(id, alt, choice, a, b, c, d, e)

devtools::use_data(cs_uncorrelated, overwrite = TRUE)

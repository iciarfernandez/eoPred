## code to prepare `DATASET` dataset goes here

mod <- readRDS(file="/Users/iciar/Documents/mod.rds")

usethis::use_data(mod,
                  internal = TRUE,
                  overwrite = TRUE,
                  compress = "bzip2")

##\This script takes data from STAGdb  accession numbers,
# and queries the associated metadata to add to the genotypes table


#load libraries
library(XML)
library(RCurl)
library(rlist)
library(rvest)
library(stringr)
library(tidyverse)
library(RMariaDB)
library(data.table)

#download table of all mlg clonal ids
data <- readHTMLTable(getURL("https://coralsnp.science.psu.edu/reports/genotypes/all",.opts = list(ssl.verifypeer = FALSE)))[[1]]

#attach the genotype db id to the table
webURL <- 'https://coralsnp.science.psu.edu/reports/genotypes/all'
page<-read_html(webURL)
links<-page %>% html_nodes("a") %>% html_attr("href")
data[,5] <-str_match(links,"genotype_id=(\\w+?)&")[,2]

#select only coloumns of interest and filter to only cervicornis
data <- data %>% 
  select(b_mlg_id = `Coral MLG Clonal ID`,
         species = `Genetic Coral Species Call`,
         b_geno_id = V5) %>%
  filter(species == 'A.cervicornis') %>%
  select(-species)
  
rm(page, links, webURL)

#grab microsatelite id
data <- split(data, seq(nrow(data)))
data <- lapply(data_list, function(x) {
  x %>% 
    mutate(microsat_id = 
                   paste(grep("^C[[:digit:]]",
                              unique(readHTMLTable(getURL(paste0("https://coralsnp.science.psu.edu/reports/samples/with_genotype?genotype_id=",b_geno_id,
                                                                 "&coral_mlg_clonal_id=",b_mlg_id),
                                                          .opts = list(ssl.verifypeer = FALSE)))[[1]][,22]),
                              value = T),
                         collapse = ", "),
           local_name = paste(unique(readHTMLTable(getURL(paste0("https://coralsnp.science.psu.edu/reports/samples/with_genotype?genotype_id=",b_geno_id,
                                                                      "&coral_mlg_clonal_id=",b_mlg_id),
                                                               .opts = list(ssl.verifypeer = FALSE)))[[1]][,6]),
                              collapse = ", "))
}) %>% rbindlist() %>% separate_rows(microsat_id)


#work with the genotypes in my db
db <- dbConnect(RMariaDB::MariaDB(), host = "localhost", user = "acdc", password = "Acervicornis",
                dbname='acdc', timeout = -1)


genos <- dbGetQuery(db, "SELECT geno_name, mlg_id, mlg_geno_id, microsat_id, alt_names FROM genotypes
                        WHERE microsat_id IS NOT NULL")

#check for microsat overlap
microsat_overlap <- genos %>%
  left_join(data, by='microsat_id') %>%
  filter((mlg_id.x != mlg_id.y) | is.na(mlg_id.x) & !is.na(mlg_id.y))

#check for mlg_id overlap
mlgID_overlap <- genos %>%
  left_join(data, by='mlg_id') %>%
  filter((mlg_geno_id.x != mlg_geno_id.y) & !is.na(mlg_id))

#check for matches in STAGdb that are not in our db
genos <- dbGetQuery(db, "SELECT geno_name, mlg_id, mlg_geno_id, microsat_id, alt_names, source_lat, source_long FROM genotypes
                        WHERE mlg_id IS NULL AND
                        source_lat IS NOT NULL")

#grab the lat long from STAGdb
stagdb_coords <- lapply(split(data, seq(nrow(data))), function(x) {
  webURL <- paste0("https://coralsnp.science.psu.edu/reports/samples/with_genotype?genotype_id=",x$mlg_geno_id,
                   "&coral_mlg_clonal_id=",x$mlg_id)
  page<-read_html(webURL)
  links<-page %>% html_nodes("a") %>% html_attr("href")
  colony_url <- paste0("https://coralsnp.science.psu.edu/reports/colonies/of_sample?colony_id=",
                       str_match(links[4],"colony_id=(\\w+?)&")[,2])
  coord <- readHTMLTable(getURL(colony_url,
                             .opts = list(ssl.verifypeer = FALSE)))[[1]][,1:2]
  
  x %>% 
    mutate(source_lat = paste(unique(coord$Latitude), collapse = ", "),
           source_long = paste(unique(coord$Longitude), collapse = ", "))
}) %>%
  bind_rows() %>%
  separate_rows(source_lat, source_long) %>%
  filter(!source_lat %in% c('9.000000','') & !source_long %in% c('9.000000','')) %>%
  mutate(lat_fixed = ifelse(as.numeric(source_lat)<26.5,
                            source_lat, source_long),
         long_fixed = ifelse(as.numeric(source_long)>26.5,
                             source_long, source_lat)) %>%
  mutate(source_lat = as.numeric(lat_fixed),
         source_long = -1*as.numeric(long_fixed)) %>%
  select(-c(lat_fixed,long_fixed))

#does the filtered genos df match any of the coords from galaxy?
#rounded off to deal with less precise coords
stagdb_coords_sub <- stagdb_coords %>%
  mutate(source_lat_rounded = round(source_lat,3),
         source_lon_rounded = round(source_long,3)) %>%
  filter(source_lat_rounded %in% round(genos$source_lat,3) |
           source_lon_rounded %in% round(genos$source_long,3)) %>%
  select(-c(source_lat_rounded, source_lon_rounded))
  




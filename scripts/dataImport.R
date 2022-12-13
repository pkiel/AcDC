#load libraries --------------------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)
library(RMariaDB)

#establish database conection ------------------------------------------
#if running on Patrick's Mac
db <- dbConnect(RMariaDB::MariaDB(),
                host = "localhost",
                user = "acdc",
                password = "Acervicornis",
                dbname='acdc',
                timeout = -1)

#functions -------------------------------------------------------------
#add genotype id
geno2ID <- function(x) {
  geno_IDs <- dbGetQuery(db, "SELECT
           geno_id, geno_name FROM genotypes")
  x %>%
    left_join(geno_IDs,by = "geno_name") %>%
    select(-geno_name) %>%
    drop_na()
}

#replace trait with it's corresponding id
trait2ID <- function(x) {
  trait_IDs <- dbGetQuery(db, "SELECT
           trait_id, trait_name AS trait FROM traits")
  x %>%
    left_join(trait_IDs,by = "trait") %>%
    select(-trait) %>%
    drop_na()
}

#replace trait with it's corresponding id
unit2ID <- function(x) {
  std_IDs <- dbGetQuery(db, "SELECT
           standard_id, unit FROM standards")
  x %>%
    left_join(std_IDs,by = "unit") %>%
    select(-unit) %>%
    drop_na()
}

#replace location with it's corresponding id
location2ID <- function(x) {
  loc_IDs <- dbGetQuery(db, "SELECT
           location_id, location_name FROM locations")
  x %>%
    left_join(loc_IDs,by = "location_name") %>%
    select(-location_name) %>%
    drop_na()
}

method2ID <- function(x) {
  #gather list of methods from db
  method_IDs <- dbGetQuery(db, "SELECT
           method_id, method FROM methods")
  #send new method names to db
  y <- x %>%
    left_join(method_IDs,by = "method") %>%
    filter(is.na(method_id)) %>%
    distinct(method)
  if(nrow(y)>0) {
    dbWriteTable(db,
                 name = "methods",
                 value = y,
                 overwrite = FALSE,
                 append = TRUE,
                 row.names = FALSE)}
  #retrieve updated method table
  method_IDs <- dbGetQuery(db, "SELECT
           method_id, method FROM methods")
  #append method_id to df
  x %>%
    left_join(method_IDs,by = "method") %>%
    select(-method) %>%
    drop_na()
}

masterID <- function(x) {
  x %>% 
    geno2ID() %>%
    trait2ID() %>%
    unit2ID() %>%
    method2ID()
}

dataset2ID <- function(x) {
  dataset_IDs <- dbGetQuery(db, "SELECT
           dataset_id, datafile_name FROM datasets")
  x %>%
    left_join(dataset_IDs,by = "datafile_name") %>%
    select(-datafile_name) %>%
    drop_na()
}
#Data -----------------------------------------------------------------
#load data
file <- file.choose()
data <- read_csv(file)

#create a list of observations with each in its own df
observations <- lapply(1:nrow(data), function(x) data[x,])
rm(data)

#gather each measurement into rows
#append specific unit and method to the corresponding traits
measurements <- lapply(observations, function(x) {
  x %>%
    mutate_all(as.character) %>%
    select(!c(location_name)) %>%
    pivot_longer(cols = !geno_name,
                 names_to = "trait",
                 values_to = "value",
                 values_drop_na = T) %>%
  #add in the units and methods
  #add additional conditional statements to match the traits present in dataset
  #first , e.g.:trait=='mass', 'g'
  mutate(unit = case_when(trait == 'tag' ~ 'text',
                          trait == 'date' ~ 'date',
                          trait == 'mass' ~ 'g',
                          trait == 'TLE' ~ 'cm',
                          TRUE ~ NA),
        method = case_when(trait == 'tag' ~ 'tag',
                           trait == 'date' ~ 'date',
                           trait == 'mass' ~ 'buoyant weight',
                           trait == 'TLE' ~ 'ruler'))
})

#change dataset,unit,method,trait,genotype to appropriate IDs
#change dataset,unit,method,trait,genotype to appropriate IDs
IDtrial <- bind_rows(measurements) %>%
  mutate(tst = row_number())
catch_errors <- IDtrial %>%
  masterID()
lost <- IDtrial %>%
  filter(!tst %in% catch_errors$tst)
if(nrow(lost)>0) {
  #catches data that is not converting to an appropriate id
  #review dataset, likely spelling error or using not accepted trait, standard, method terms
  View(lost)
} else {
  rm(IDtrial, catch_errors, lost)
  measurements <- lapply(measurements, masterID) 
}

#remove measurement level data from observations
observations <- lapply(observations, function(x){
  x %>%
    select(geno_name, datafile_name, location_name) %>%
    mutate_all(as.character) %>%
    #add the short hand datafile_name
    mutate(datafile_name = )
})
#apply genoID to observations
observations_test <-lapply(observations, geno2ID)
#apply locationID to locations
observations_test <-lapply(observations_test, location2ID)
#apply datasetID to observations
observations_test <- lapply(observations_test, dataset2ID)
observations_lost <- mapply(function(x,y) {
  x <- x %>% mutate(tst = row_number())
  y <- y %>% mutate(tst = row_number()) %>%
    filter(!tst %in% x$tst) %>%
    select(-tst)
  return(y)
}, observations_test, observations, SIMPLIFY = F) %>%
  rbindlist()

if(nrow(observations_lost)>0) {
  #catches data that is not converting to an appropriate id
  #review dataset, likely spelling error or using not accepted genet name
  View(observations_lost)
} else {
  observations <- observations_test
  rm(observations_lost, observations_test)
}

rm(geno2ID, masterID, method2ID, trait2ID, unit2ID, location2ID)

#remove empty dfs from lists
measurements <- measurements[sapply(measurements, function(x) dim(x)[1]) > 0]
observations <- observations[sapply(observations, function(x) dim(x)[1]) > 0]

#assign observation_id to each row
id_start <- tail(arrange(dbReadTable(db, "observations"), observation_id),
                 1)[,1] # retrieves last observation id
append_obsID <- function(x,y) {
  x %>% mutate(observation_id = id_start + y) #adds index# to df to last obs_id for new obs_id
}
measurements <- mapply(append_obsID, measurements, seq_along(measurements), SIMPLIFY = F)
rm(id_start, append_obsID)

#push observation dfs to observation table
  dbWriteTable(db,
               name = "observations",
               value = rbindlist(observations),
               overwrite = FALSE,
               append = TRUE,
               row.names = FALSE)

#push measurements dfs to observation table
  dbWriteTable(db,
               name = "measurements",
               value = rbindlist(measurements),
               overwrite = FALSE,
               append = TRUE,
               row.names = FALSE)
  
  
dbDisconnect(db)

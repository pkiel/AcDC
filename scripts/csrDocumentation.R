#This script takes contributor supplied Coral Sample Registry (CSR) accession numbers,
# and queries the associated metadata to add to the genotypes table

#load libaries
library(tidyverse)
library(RMariaDB)
library(httr)
library(jsonlite)
library(data.table)
library(lubridate)

#must have a sql user with -u=acdc, -p=Acervicornis, alter as needed
db <- dbConnect(RMariaDB::MariaDB(), host = "localhost", user = "acdc", password = "Acervicornis",
                dbname='acdc', timeout = -1)

#create token to connect to Coral Sample Registry database, lasts 12 hours
#this is a private token with url, user and password redacted here for protection
#please contact patrick.kiel@noaa.gov for more information
token <- content(httr::POST(url = "",#redacted
           add_headers(Authorization = ""), #redacted
           body = list(grant_type = "password", username = "",#redacted
                       password = ""),#redacted
           encode = "form"))[[1]]

#example where I pull all contributed CSR in the database and query the metadata
genos <- dbGetQuery(db, "SELECT geno_id, CSR_accession FROM genotypes
                    WHERE CSR_accession IS NOT NULL")

#vector of accession #'s to query from the CSR database
accession <- genos$CSR_accession

#function to create dataframe from the API call
grabGenos <- function(x) {
  call <- paste("",#redacted, same as the URL above
                x,
                sep = "")
  data <- as.data.frame(fromJSON(content(httr::GET(call, add_headers(Authorization = paste("Bearer", token))), "text"), flatten = TRUE))
  data %>% mutate(CSR_accession = x, .before = "type") %>%
    rename(geno_name = localGenotypeName) %>%
    mutate(source_long = as.numeric(str_sub(longitude,1,-2)),
           source_lat = as.numeric(str_sub(latitude,1,-2)),
           source_date = date(collectionDate),
           sexual_recruit = ifelse(type == "SEXUAL_RECRUIT", 1, 0)) %>%
    select(-type)
}

#query accession #'s supplied above
metadata <- mapply(grabGenos, accession, SIMPLIFY = F) %>% bind_rows()

#example of generated dataframe from above on 10 random CSR accession #s
metadata <- structure(list(geno_id = c(1:5,10:14),
                          CSR_accession = c("5c98b2aa-b95c-d02c-1a8e-9732b437adf4", 
                                   "c659553a-e522-d6a8-957d-ad46ea30d52a", 
                                   "10051424-b078-44c9-3738-856e93701356", 
                                   "60c3c98c-e6f7-5779-6d64-a818f37aa854", 
                                   "943c45b0-5db6-e5dc-fe6a-ee9f6edc1d67", 
                                   "0ccf8581-fe36-a3a9-1d32-eeb28bf6e7c2", 
                                   "1aa6c453-ff50-9f93-9c28-73e1715fb7cd", 
                                   "da4b435b-f7e3-a001-b9d9-4f8a8a8d4374", 
                                   "15248325-efd0-9a0a-6c45-669a867faa72", 
                                   "4d152000-3436-340b-815d-cb8673f625e4"),
                           sexual_recruit = rep(0,10), 
                          geno_name= c("Acerv 1 (A-AC)", "Acerv 2 (B-AC)", 
                                                 "Acerv 3 (C-AC)", "Acerv 4 (D-AC)", 
                                                 "Acerv 5 (E-AC)", "BC-1", 
                                                 "BC-11", "BC-8A", 
                                                 "BC-8B", "BC-E"),
                           source_long = c(-80.1926, -80.19882, -80.20057,
                                         -80.1854, -80.1854, -80.08962,
                                         -80.08962, -80.08933, -80.08933,
                                         -80.08747),
                           source_lat = c(25.3391, 25.3129, 25.30838,
                                        25.31895, 25.31895, 26.1763,
                                        26.1763, 26.18298, 26.18298, 26.19802),
                           origin_nursery = rep("University of Miami",10),
                           source_date = structure(c(14245, 14245, 14245, 14245, 
                                                        14245, 17474, 17474, 17474, 17474,
                                                        17721), class = "Date")),
                      row.names = c(NA,-10L), class = "data.frame")

updateStatement <- function(x) {
  
  paste0("UPDATE `acdc`.`genotypes` SET `sexual_recruit` = '",x[3],"', ",
                          "`geno_name` = '", x[4], "', ",
                          "`source_long` = '", x[5], "', ",
                          "`source_lat` = '", x[6], "', ",
                          "`origin_nursery` = '", x[7], "', ",
                          "`source_date` = '", x[8], "' ",
                          "WHERE (`CSR_accession` ='", x[2], "');")
}

#send multiple updates
dbSendQueries <- function(con,sql){
  dummyfunction <- function(sql,con){
    dbSendQuery(con,sql)
  }
  lapply(sql, dummyfunction, con)
}

#update the genotypes table with the metadata
dbSendQueries(db, sql = apply(metadata, 1, updateStatement))

dbDisconnect(db)


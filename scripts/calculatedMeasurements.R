#libraries
library(tidyverse)
library(data.table)
library(RMariaDB)
library(lubridate)
library(rstudioapi)

#set wd
setwd(dirname(getActiveDocumentContext()$path))

#db connection
#if running on Patrick's Mac
db <- dbConnect(RMariaDB::MariaDB(),
                host = "localhost",
                user = "acdc",
                password = "Acervicornis",
                dbname='acdc',
                timeout = -1)

#functions -------------------------------------------------------------
#for bleaching trait rates standardize to # of heating days
heating_days <- function(x) {
  #Degree heating days, time integrated measure of bleaching stress
  #regional limit taken as 30.5C, D.P. Manzello et al.(2007) 10.1016/j.marpolbul.2007.08.009
  x = x %>% 
    mutate(temperature_value = as.numeric(temperature_value),
           dhd = temperature_value-30.5) %>%
    filter(temperature_value > 30.5) 
  sum(x$dhd)
} 

#calculate traits
calculatedTraits <- function(corals) { 
  #create day col and remove empty cols
  corals <- lapply(corals, function(x){
    x <- x %>%
      mutate(month_value = ifelse(has_name(x, 'month_value'), as.integer(month_value),
                                     NA),
             date_value = parse_date(date_value, format="%Y-%m-%d"),
             days_value = ifelse(has_name(x, 'days_value'), as.integer(days_value),
                                    NA)) %>%
      arrange(month_value, date_value, days_value) %>%
      mutate(day = ifelse(!is.na(month_value),
                          (month_value-first(month_value))*30, #months * 30 days per month
                          ifelse(!is.na(date_value),
                                 as.numeric(date_value - first(date_value)),
                                 ifelse(!is.na(days_value),
                                        days_value-first(days_value), NA))))
    x <- Filter(function(y)!all(is.na(y)), x)
  })
  
  #return this df
  corals <- lapply(corals, function(x) {
    #create new data filters
    new_filters <- c(NA)
    new_filters <- if(has_name(x,'stress hardened_value')) {
      append(new_filters, "Stress Hardening")
    } else {new_filters}
    new_filters <- if(has_name(x,'temperature_value')) {
      z <- x$temperature_value
      z <- z[!is.na(z)]
      if(any(z>30.5)) { 
        append(new_filters, "Bleaching Stress")
      } else {append(new_filters, "ambient temperature")}
    } else {append(new_filters, "ambient temperature")}
    new_filters <- if(has_name(x, 'pH_value')) {
      z <- x$pH_value
      z <- z[!is.na(z)]
      if(any(z<7.9)) {
        append(new_filters, "OA Stress")
      } else {append(new_filters, "ambient pH")}
    } else {append(new_filters, "ambient pH")}
    new_filters <- if(!is.null(x$date_value) & all(!is.na(x$date_value)) & max(x$day, na.rm = T)<=185) {
        z <-x$date_value
      
        #wet season May 1 - Oct 31 (6 months)
        if(all(between(month(z),5,10))) {
          append(new_filters, 'wet season')
        } else if(all(month(z)<5 | month(z)>10)) {
          append(new_filters, 'dry season')
          } else{new_filters}
    }else{new_filters}
    
    #apply season filter to 6 month growth stats when total growth > 6 months
    #remove this filter from all traits except 6 month traits
    if(any(between(x$day,120,240)) & !is.null(x$date_value) & all(!is.na(x$date_value))) {
      z <- c(first(x$date_value),
             x$date_value[which(between(x$day,120,240) & abs(x$day-182)==min(abs(x$day-182)))])
    
      #wet season May 1 - Oct 31 (6 months)
      if(all(between(month(z),5,10))) {
        append(new_filters, 'wet season 6 MONTHS')
      } else if(all(month(z)<5 | month(z)>10)) {
        append(new_filters, 'dry season 6 MONTHS')
      } else{new_filters}
      } else {new_filters}
    
    new_filters <- new_filters[!is.na(new_filters)]
    
   data.frame(
      #coral metadata
      geno_id <- ifelse(has_name(x,'geno_id'), unique(x$geno_id), NA),
      tag <- ifelse(has_name(x, 'tag_value'), unique(x$tag_value), NA),
      dataset_id <- ifelse(has_name(x, 'dataset_id'), unique(x$dataset_id), NA),
      filters <- paste(unique(x$filters), paste(new_filters, collapse="; "), sep = "; "),
      location_id <- ifelse(has_name(x, 'location_id'), unique(x$location_id), NA),
      
      #delta metrics
      growth.cm <- ifelse(has_name(x, 'TLE_value') & all(as.numeric(na.omit(x$TLE_value)) == cummax(as.numeric(na.omit(x$TLE_value)))) & diff(range(cummax(as.numeric(na.omit(x$TLE_value)))))>0,
                          max(as.numeric(x$TLE_value), na.rm = T) - as.numeric(na.omit(x$TLE_value)[1]),
                          NA),
      new.mass <- ifelse(has_name(x, 'mass_value'),
                         max(as.numeric(x$mass_value), na.rm = T) - as.numeric(na.omit(x$mass_value)[1]),
                         NA),

      #calcification days
      calc.days <- ifelse(!is.na(new.mass),
                          as.numeric(x$day[which(as.numeric(x$mass_value) == as.numeric(last(na.omit(x$mass_value))))]-x$day[which(as.numeric(x$mass_value) == as.numeric(first(na.omit(x$mass_value))))]),
                          NA),
      #TLE days
      tle.days <- ifelse(!is.na(growth.cm),
                         as.numeric(x$day[which(as.numeric(x$TLE_value) == as.numeric(last(na.omit(x$TLE_value))))]-x$day[which(as.numeric(x$TLE_value) == as.numeric(first(na.omit(x$TLE_value))))]),
                  NA),
      
      #volume days
      volume.days <- ifelse(has_name(x, 'colony volume_value'),
                        as.numeric(x$day[which(as.numeric(x$`colony volume_value`) == as.numeric(last(na.omit(x$`colony volume_value`))))]-x$day[which(as.numeric(x$`colony volume_value`) == as.numeric(first(na.omit(x$`colony volume_value`))))]),
                        NA),
      
      #intersitital space days
      is.days <-  ifelse(has_name(x, 'interstitial space volume_value'),
                    as.numeric(x$day[which(as.numeric(x$`interstitial space volume_value`) == as.numeric(last(na.omit(x$`interstitial space volume_value`))))]-x$day[which(as.numeric(x$`interstitial space volume_value`) == as.numeric(first(na.omit(x$`interstitial space volume_value`))))]),
                    NA),
      
      #methodology used for colony volumetric
      volume_method <- ifelse(has_name(x,'colony volume_method'),
                                       unique(x$`colony volume_method`),
                                       NA),
                              
      #bleaching scores
      #exponential decay of bleaching metric expressed in intuitive per cent change in bleaching metric per degree heating day
      #where min value represents the final and max denotes the original, following cummmin/cummax filter
      bleachingIPAM <- ifelse(has_name(x, 'maximum photosynthetic yield_value') & diff(range(cummin(as.numeric(na.omit(x$`maximum photosynthetic yield_value`)))))>0,
                              ((min(as.numeric(x$`maximum photosynthetic yield_value`), na.rm = T)/max(as.numeric(x$`maximum photosynthetic yield_value`),na.rm = T))^(1/heating_days(x))-1)*100,
                              NA),
      
      bleachingColorScore <- ifelse(has_name(x, 'color score_value') & diff(range(cummin(as.numeric(na.omit(x$`color score_value`)))))>0,
                                    ((min(as.numeric(x$`color score_value`), na.rm = T)/max(as.numeric(x$`color score_value`),na.rm = T))^(1/heating_days(x))-1)*100,
                                    NA), 
      
      bleachingRscore <- ifelse(has_name(x, 'bleaching R-score_value') & diff(range(cummax(as.numeric(na.omit(x$`bleaching R-score_value`)))))>0,
                                (((max(as.numeric(x$`bleaching R-score_value`), na.rm = T)/min(as.numeric(x$`bleaching R-score_value`),na.rm = T))^(1/heating_days(x))-1)*-100),
                                NA),
      
      #wound healing
      healing_rate <- ifelse(has_name(x,'lesion healing_value'),
                             yes = (max(as.numeric(x$`lesion healing_value`), na.rm = T)/
                                      as.numeric(x %>% filter(as.numeric(x$`lesion healing_value`) == max(as.numeric(x$`lesion healing_value`), na.rm = T)) %>%
                                                   pull(day) %>% first())),
                             no = NA),

      #growth rates
      #capture the first 6 months of growth for corals grown > 6 months
      `6-month linear growth` <- ifelse(!is.na(growth.cm) & any(between(x$day,150,210)),
                                        (as.numeric(x$TLE_value[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]) -
                                           as.numeric(x$TLE_value[1]))/
                                          as.numeric(x$TLE_value[1])/
                                          x$day[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]*182,NA),
      
      #capture the first 12 months of growth for corals grown > 12 months
      annualLinearGrowth <- ifelse(!is.na(growth.cm) & any(between(x$day,330,400)),
                                   (as.numeric(x$TLE_value[which(between(x$day,330,400) & abs(x$day-365)==min(abs(x$day-365)))])-as.numeric(x$TLE_value[1]))/as.numeric(x$TLE_value[1])/
                                     x$day[which(between(x$day,330,400) & abs(x$day-365)==min(abs(x$day-365)))]*365,NA),
      
      #SGR = per cent increase in TLE per day, exponential growth model
      SGR <- ifelse(!is.na(growth.cm) & tle.days>=90,
                    ((max(as.numeric(x$TLE_value), na.rm = T)/as.numeric(x$TLE_value[1]))^(1/tle.days)-1)*100,
                    NA),
      
      #unstandardized calcification
      DailyCalcification <- ifelse(!is.na(new.mass),
                                   (new.mass*1000/calc.days),
                                   NA),
      
      #mg per day standardized per initial gram
      massnormalizedDailyCalcification <- ifelse(!is.na(new.mass),
                                             (new.mass*1000/as.numeric(x$mass_value[1])/calc.days), NA),
      #mg per day standardized per  surface area
      SAnormalizedDailyCalcification <- ifelse((!is.na(new.mass)&has_name(x, 'SA_value')),
                                               (new.mass*1000/min(as.numeric(x$SA_value), na.rm = T)/calc.days), NA),

      
      #colony volumetric growth

      `6-month interstitial space growth` <- ifelse(has_name(x, 'interstitial space volume_value') & all(as.numeric(x$`interstitial space volume_value`) == cummax(as.numeric(x$`interstitial space volume_value`))) & any(between(x$day,150,210)),
                                                    (as.numeric(x$`interstitial space volume_value`[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]) -
                                                       as.numeric(x$`interstitial space volume_value`[1]))/
                                                      as.numeric(x$`interstitial space volume_value`[1])/
                                                      x$day[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]*182,NA),
      
      `6-month colony volumetric growth` <- ifelse(has_name(x, 'colony volume_value') & all(as.numeric(x$`colony volume_value`) == cummax(as.numeric(x$`colony volume_value`))) & any(between(x$day,150,210)),
                                                   (as.numeric(x$`colony volume_value`[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]) -
                                                      as.numeric(x$`colony volume_value`[1]))/
                                                     as.numeric(x$`colony volume_value`[1])/
                                                     x$day[which(between(x$day,150,210) & abs(x$day-182)==min(abs(x$day-182)))]*182,NA),
      
      
      `annual interstitial space growth` <-  ifelse(has_name(x, 'interstitial space volume_value') & all(as.numeric(x$`interstitial space volume_value`) == cummax(as.numeric(x$`interstitial space volume_value`))) &  any(between(x$day,330,400)),
                                                    (as.numeric(x$`interstitial space volume_value`[which(between(x$day,330,400) & abs(x$day-365)==min(abs(x$day-365)))]) -
                                                       as.numeric(x$`interstitial space volume_value`[1]))/
                                                      as.numeric(x$`interstitial space volume_value`[1])/
                                                      x$day[which(between(x$day,330,400) & abs(x$day-365)==min(abs(x$day-365)))]*365,NA),
      
      `annual colony volumetric growth` <- ifelse(has_name(x, 'colony volume_value') & all(as.numeric(x$`colony volume_value`) == cummax(as.numeric(x$`colony volume_value`))) &  any(between(x$day,330,400)),
                                                  (max(as.numeric(x$`colony volume_value`), na.rm = T) -
                                                     as.numeric(x$`colony volume_value`[1]))/
                                                    as.numeric(x$`colony volume_value`[1])/
                                                    volume.days*365, NA),
      
      `interstitial space SGR` <- ifelse(has_name(x, 'interstitial space volume_value') & all(as.numeric(x$`interstitial space volume_value`) == cummax(as.numeric(x$`interstitial space volume_value`))) & is.days>=90,
                                         ((max(as.numeric(x$`interstitial space volume_value`), na.rm = T)/as.numeric(x$`interstitial space volume_value`[1]))^(1/is.days)-1)*100,
                                         NA),
      
      `colony volumetric SGR` <- ifelse(has_name(x, 'colony volume_value') & all(as.numeric(x$`colony volume_value`) == cummax(as.numeric(x$`colony volume_value`))) & volume.days>=90,
                                        ((max(as.numeric(x$`colony volume_value`), na.rm = T)/as.numeric(x$`colony volume_value`[1]))^(1/volume.days)-1)*100,
                                        NA)
    ) %>%
     #remove the days cols
     select(- c(8:11)) %>%
      setNames(., c("geno_id","tag","dataset_id", "filters", "location_id","growth.cm",
                    "new.mass", "volume_method", "bleaching photochemical efficiency", "bleaching color score",
                    "bleaching R-score", "healing rate", "6-month linear growth", "annual linear growth",
                    "specific growth rate", "daily calcification", "mass normalized daily calcification",
                    "SA normalized daily calcification", "6-month interstitial space growth",
                    "6-month colony volumetric growth", "annual interstitial space growth",
                    "annual colony volumetric growth", "interstitial space SGR", "colony volumetric SGR"))
  }) %>% rbindlist() %>%
    filter(across(everything(), ~ !is.infinite(.))) %>%
    #cast long
    pivot_longer(
      `bleaching photochemical efficiency`:`colony volumetric SGR`,
      names_to = "stat",
      values_to = "value",
      values_drop_na = T) %>%
    select(geno_id,trait = stat, value, dataset_id, filters, location_id, volume_method) %>%
    #filter out some data
    filter((value > -50 & trait=='bleaching photochemical efficiency') | trait != 'bleaching photochemical efficiency') %>%
    filter((value > 0 & grepl('growth|SGR|calcification', trait)) | !grepl('growth|SGR|calcification', trait)) %>% #filter out no growth/negative growth values
    #apply methods, units
    mutate(method = case_when(trait %in% c('daily calcification','mass normalized daily calcification','SA normalized daily calcification')   ~ 'buoyant weight',
                              trait %in% c('healing rate','bleaching R-score') ~ 'image analysis',
                              trait == 'bleaching color score' ~ 'Coral Watch Bleaching Card',
                              trait == 'bleaching photochemical efficiency' ~ 'I-PAM',
                              trait %in% c('6-month linear growth','annual linear growth','specific growth rate') ~ 'ruler',
                              volume_method == '3D photogrammetry' ~ '3D photogrammetry',
                              volume_method == 'ellipsoid volume' ~ 'ellipsoid volume'),
           unit = case_when(trait == 'daily calcification' ~ 'mg day^-1',
                            trait == 'mass normalized daily calcification' ~ 'mg day^-1 g^-1',
                            trait == 'SA normalized daily calcification' ~'mg day^-1 cm^-2',
                            trait == 'healing rate' ~ '% healed day^-1 cm^-2',
                            trait %in% c('bleaching color score','bleaching photochemical efficiency','bleaching R-score') ~ '% change heating day^-1',
                            trait %in% c('6-month linear growth','6-month interstitial space growth','6-month colony volumetric growth') ~ 'proportionate growth',
                            trait %in% c('annual linear growth','annual interstitial space growth','annual colony volumetric growth') ~ 'annual proportionate growth',
                            trait %in% c('specific growth rate','interstitial space SGR','colony volumetric SGR') ~ '% change day^-1'),
           .after="trait") %>%
    select(-volume_method)
} 

#replace trait with it's corresponding id
trait2ID <- function(x) {
  trait_IDs <- dbGetQuery(db, "SELECT
           trait_id, trait_name AS trait FROM traits")
  x %>%
    left_join(trait_IDs,by = "trait") %>%
    select(-trait) %>%
    drop_na(trait_id)
}

#replace trait with it's corresponding id
unit2ID <- function(x) {
  std_IDs <- dbGetQuery(db, "SELECT
           standard_id, unit FROM standards")
  x %>%
    left_join(std_IDs,by = "unit") %>%
    select(-unit) %>%
    drop_na(standard_id)
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
    drop_na(method_id)
}

filter2ID <- function(x) {
  #gather list of methods from db
  filter_IDs <- dbGetQuery(db, "SELECT
           filter_id, filters FROM filters")
  #send new method names to db
  y <- x %>%
    left_join(filter_IDs,by = "filters") %>%
    filter(is.na(filter_id)) %>%
    distinct(filters)
  if(nrow(y)>0) {
    dbWriteTable(db,
                 name = "filters",
                 value = y,
                 overwrite = FALSE,
                 append = TRUE,
                 row.names = FALSE)}
  #retrieve updated method table
  filter_IDs <- dbGetQuery(db, "SELECT
           filter_id, filters FROM filters")
  #append method_id to df
  x %>%
    left_join(filter_IDs,by = "filters") %>%
    select(-filters) %>%
    drop_na(filter_id)
}

masterID <- function(x) {
  x %>% 
    trait2ID() %>%
    unit2ID() %>%
    method2ID() %>%
    filter2ID()
}


#grab raw data----
corals <- dbGetQuery(db, "Select 
                              o.geno_id AS geno_id,
                              m.observation_id AS observation_id,
                              o.location_id AS location_id,
                              l.type AS filters,
                              trait_name AS trait,
                              trait_class,
                            	value,
                            	unit,
                            	method,
                            	dataset_id
                              FROM
	                            measurements AS m
                            	INNER JOIN
                            	observations AS o
                            	ON m.observation_id = o.observation_id
                              INNER JOIN
                              standards AS s
                              ON m.standard_id = s.standard_id
                              INNER JOIN
                              traits AS t
                              ON m.trait_id = t.trait_id
                              INNER JOIN
                              locations AS l
                              ON o.location_id = l.location_id
                              INNER JOIN
                              methods as meth
                              ON m.method_id = meth.method_id
                                WHERE (trait_name IN ('TLE', 'mass', 'SA', 'bleaching R-score', 'color score',
                                               'maximum photosynthetic yield', 'lesion healing', 'interstitial space volume',
                                               'lesion area', 'colony volume') OR trait_class = 'contextual')
                                AND m.observation_id IN
              										(SELECT observation_id 
              										FROM measurements m
              										INNER JOIN traits t
              										ON m.trait_id=t.trait_id
              										AND t.trait_name IN ('TLE', 'branches', 'mass', 'SA', 'bleaching R-score', 'color score',
                                               'maximum photosynthetic yield', 'lesion healing', 'interstitial space volume',
                                               'lesion area', 'colony volume'))
                                ORDER BY m.observation_id")

#calculate rates using functions-----
y <- corals %>%
  distinct() %>% 
  select(-c(trait_class)) %>%
  pivot_wider(names_from = trait,
              values_from = c(value, unit, method),
              names_glue = "{trait}_{.value}") %>%
  split(with(., interaction(tag_value,dataset_id,location_id)), drop = TRUE)

corals <- calculatedTraits(y)

rm(y, calculatedTraits,heating_days)

#grab the point measurement data----
x <- dbGetQuery(db, "Select 
                              m.observation_id,
                              o.geno_id AS geno_id,
                              o.location_id,
                              l.type AS filters,
                              trait_name AS trait,
                              trait_class AS class,
                            	value,
                            	unit,
                            	method,
                            	dataset_id
                              FROM
	                            measurements AS m
                            	INNER JOIN
                            	observations AS o
                            	  ON m.observation_id = o.observation_id
                              INNER JOIN
                              standards AS s
                                ON m.standard_id = s.standard_id
                              INNER JOIN
                              locations AS l
                              ON o.location_id = l.location_id
                              INNER JOIN
                              traits AS t
                                ON m.trait_id = t.trait_id
                              INNER JOIN
                              methods as meth
                                ON m.method_id = meth.method_id
                              WHERE (trait_name IN ('lipid density', 'bulk density',
                                                    'dark calcification', 'light calcification', 'respiration', 
                                                    'photosynthesis', 'symbiont density', 'chlorophyll-a density',
                                                    'symbiont chlorophyll-a density', 'tissue dry weight', 'ED50',
                                                    'low credible interval', 'median relative risk', 'high credible interval',
                                                    'polyp density', 'oocyte density', 'oocyte length',
                                                    'oocyte width', 'oocyte volume', 'bundle sperm density',
                                                    'bundle egg density') OR trait_class = 'contextual')
                              AND m.observation_id IN
                            		(SELECT observation_id 
                            		 FROM measurements m
                            		 INNER JOIN traits t
                            		 ON m.trait_id=t.trait_id
                            		 AND t.trait_name IN ('lipid density', 'bulk density',
                                                    'dark calcification', 'light calcification', 'respiration', 
                                                    'photosynthesis', 'symbiont density', 'chlorophyll-a density',
                                                    'symbiont chlorophyll-a density', 'tissue dry weight', 'ED50',
                                                    'low credible interval', 'median relative risk', 'high credible interval',
                                                    'polyp density', 'oocyte density', 'oocyte length',
                                                    'oocyte width', 'oocyte volume', 'bundle sperm density',
                                                    'bundle egg density'))
                             ORDER BY m.observation_id")

#append new filters
x <- x %>%
  distinct() %>% #remove duplicate rows, only caused by when calling SA for some reason
  pivot_wider(names_from = trait,
              values_from = c(value, class),
              names_glue = "{trait}_{.value}") %>%
  group_by(observation_id) %>%
  mutate(across(c(temperature_value,pH_value), ~as.numeric(.))) %>%
  fill(c(pH_value, temperature_value), .direction = "downup") %>%
  mutate(filters = case_when(
      any(pH_value<7.9 & temperature_value>30.5 , na.rm = T) ~ paste(filters, "Bleaching Stress", "OA Stress",  sep="; "),
      any(pH_value<7.9, na.rm = T) ~ paste(filters, "ambient temperature", "OA Stress", sep="; "),
      any(temperature_value>30.5, na.rm = T) ~ paste(filters, "Bleaching Stress", "ambient pH", sep="; "),
      TRUE ~ paste(filters, "ambient temperature", "ambient pH", sep="; "))
  ) %>%
  ungroup() %>%
  #fill in blanks
  mutate(pH_value = ifelse(method!="pH probe",NA,pH_value),
         temperature_value = ifelse(method!="temperature gauge",NA,temperature_value)) %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = -c(observation_id, geno_id,filters,unit,method,dataset_id,location_id),
               names_to = c("trait",".value"),
               values_drop_na = T,
               names_pattern = "(.*)_(.*)") %>%
  filter(class == 'organismal') %>%
  select(-class)

#photosynthesis + respiration
p_r <- x %>% filter(trait %in% c('photosynthesis', 'respiration')) %>%
  pivot_wider(names_from = trait,
              values_from = value) %>%
  mutate(`P:R` = ifelse(!is.na(photosynthesis) & !is.na(respiration),
                        abs(as.numeric(photosynthesis)/as.numeric(respiration)),
                        NA)) %>%
  mutate(across(respiration:`P:R`, as.character)) %>%
  pivot_longer(cols = respiration:`P:R`,
               names_to = "trait",
               values_to = "value") %>%
  mutate(unit = ifelse(trait == 'P:R', 
                       'dimensionless',
                       unit))

x <- x %>% 
  filter(!trait %in% c('photosynthesis', 'respiration')) %>%
  bind_rows(p_r) %>%
  arrange(observation_id) %>%
  select(-observation_id)

rm(p_r)

corals <- corals %>% mutate(across(everything(), as.character))
x <- x %>% mutate(across(everything(), as.character))
corals <- bind_rows(corals, x)
rm(x)

#change everything to ids----
test <- corals %>% mutate(tst = row_number())
test <- masterID(test)
lost <- test %>% filter(!tst %in% test$tst)
if(nrow(lost)>0){View(lost)}
corals <- test %>% select(-tst) %>% drop_na()
rm(test, lost, filter2ID, method2ID, trait2ID, unit2ID, masterID)

#push df to the database-----
#delete all data from the combined_measurements table
dbExecute(db, "DELETE FROM calculatedmeasurements")

#set combined_measurements table auto_increment to 1
dbExecute(db, "ALTER TABLE `acdc`.`calculatedmeasurements` 
AUTO_INCREMENT = 1;")

#write table
dbWriteTable(db,
             name = "calculatedmeasurements",
             value = corals,
             overwrite = FALSE,
             append = TRUE,
             row.names = FALSE)

dbDisconnect(db)


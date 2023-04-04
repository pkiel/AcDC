#head ----
#NOAA Acropora cervicornis Data Coordination hub (AcDC)
#UI for RShiny Database App
#Updated 10/30/2022
#load libraries------
library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyjs)
library(shinyBS)
library(DT)
library(tidyverse)
library(data.table)
library(RMariaDB)
library(glue)
library(leaflet)
library(leaflet.extras)
library(sp)
library(spsComps)
library(openxlsx)

#establish db connection + collect geno names before app begins----
#connect to database w/ query only account
db <- dbConnect(RMariaDB::MariaDB(),
                host = "localhost",
                user = "acdc",
                password = "Acervicornis",
                dbname='acdc',
                timeout = -1)

biomechanicalTraits <- c('bulk density', 'polyp density')
bleachingTraits <- c('bleaching photochemical efficiency', 'bleaching color score',
                     'bleaching R-score', 'ED50', 'symbiont chlorophyll-a density',
                     'chlorophyll-a density', 'symbiont density', 'composite bleaching resistance index')
diseaseTraits <- c('high credible interval', 'low credible interval', 
                   'median relative risk')
growthTraits <- c('6-month linear growth','annual linear growth', 'annual linear growth per branch',
                  'specific growth rate','6-month colony volumetric growth', 
                  '6-month interstitial space growth','annual colony volumetric growth', 
                  'annual interstitial space growth','colony volumetric SGR', 
                  'interstitial space SGR', 'daily calcification', 
                  'mass normalized daily calcification', 
                  'SA normalized daily calcification',
                  'dark calcification', 'light calcification', 'composite growth index')
hostPhysiologyTraits <- c('P:R', 'respiration', 'photosynthesis', 'lipid density', 'tissue dry weight')
reproductionTraits <- c('oocyte density', 'oocyte volume', 'bundle sperm density',
                        'bundle egg density')
woundHealingTraits <- c('healing rate')


datasets <- dbGetQuery(db, "SELECT
                            o.dataset_id AS dataset_id,
                            datafile_name,
                            title,
                            author,
                            email,
                            ORCID,
                            DOI,
                            description,
                            meth_proc,
                            submission_date,
                            private,
                            GROUP_CONCAT(DISTINCT trait_name ORDER BY trait_name SEPARATOR ', ') AS traits,
                            GROUP_CONCAT(DISTINCT geno_name ORDER BY geno_name SEPARATOR ', ') AS genotypes,
                            remaining_tissue
                       FROM datasets
                       INNER JOIN observations o
                        ON o.dataset_id = datasets.dataset_id
                       INNER JOIN measurements m
                        ON m.observation_id = o.observation_id
                       INNER JOIN traits t
                        ON t.trait_id = m.trait_id
                       INNER JOIN genotypes g
                        ON g.geno_id = o.geno_id
                       GROUP BY dataset_id") %>%
  mutate(across(2:ncol(.), ~ifelse(. == 1, "TRUE", 
                                   ifelse(. == 0, "FALSE", .))),
         submission_date = as.Date(as.POSIXct(submission_date, origin="1970-01-01", tz="GMT")),
         DOI = ifelse(!is.na(DOI),
                      paste0("<a href='https://doi.org/",DOI,"' target='_blank'>",DOI,"</a>"),
                      NA),
         ORCID = ifelse(!is.na(ORCID),
                        paste0("<a href='https://orcid.org/",ORCID,"' target='_blank'>",ORCID,"</a>"),
                        NA))

#grab all traits before app even begins----
corals <- dbGetQuery(db, "Select 
                              geno_name AS genotype,
                              mlg_id,
                              CSR_accession,
                              trait_name AS trait,
                            	value,
                            	unit,
                            	method,
                            	datafile_name,
                            	df.filters,
                            	location_name,
                            	location_lat AS lat,
                            	location_long AS 'long'
                              FROM
	                            calculatedmeasurements AS m
                              INNER JOIN
                              genotypes AS g ON m.geno_id = g.geno_id
                              INNER JOIN
                              traits AS t ON m.trait_id = t.trait_id
                              INNER JOIN
                              standards AS s ON m.standard_id = s.standard_id
                              INNER JOIN
                              methods as meth ON m.method_id = meth.method_id
                              INNER JOIN
                              datasets AS dt ON m.dataset_id = dt.dataset_id
                              INNER JOIN
                              locations AS l ON m.location_id = l.location_id
                              INNER JOIN
                              filters AS df ON m.filter_id = df.filter_id") %>%
  mutate(value = as.numeric(value))

#list of geno names
geno_names <- corals %>%
  select(geno_name = genotype, CSR_accession, mlg_id) %>%
  distinct() %>%
  arrange(geno_name)

#reusable functions -----
#create reusable popover button
infoBtn <- function(id) {
  actionButton(id,
               label = "",
               icon = icon("question"),
               style = "info",
               size = "extra-small",
               class='btn action-button btn-info btn-xs shiny-bound-input'
  )
}

#create reusable modal plot
#calculate mean-1sd
lowerSD <- function(x) {
  mean(x)-sd(x)
}
#calculate mean+1sd
upperSD <- function(x) {
  mean(x)+sd(x)
}
modalPlot <- function(pop,geno,z) {
  pop %>%
    ggplot(aes(value,trait)) +
    #create transparent graph, color in black lines below the graph only
    stat_boxplot(geom = "errorbar", width = 0.4, size=3, color = "transparent") +
    stat_summary(fun=mean, geom="segment", aes(xend=..x.., yend=1, y = 0.8), size=3) +
    stat_summary(fun=min, geom="segment", aes(xend=..x.., yend=1, y = 0.8), size=3) +
    stat_summary(fun=max, geom="segment", aes(xend=..x.., yend=1, y = 0.8), size=3) +
    #white rectangle canvas with black border
    geom_rect(aes(xmin=min(value), xmax=max(value), ymin=0.9, ymax=1.1),
              fill="white", color="black", size=3) +
    #annotate line
    stat_summary(data = geno,fun=mean,
                 geom="segment", aes(xend=..x.., y = 1, yend=1.2), size=1,
                 color = "black") +
    stat_summary(fun=lowerSD, geom="segment", aes(xend=..x.., yend=0.95, y = 0.85), size=1) +
    stat_summary(fun=upperSD, geom="segment", aes(xend=..x.., yend=0.95, y = 0.85), size=1) +
    #point + error bar for selected genotype
    stat_summary(data = geno, fun=mean,
                 geom="point",  size=8, shape = 21, color = 'black', stroke = 2,
                 fill = ifelse(z > 1,'green',
                               ifelse(z < -1,'red',
                                      'yellow'))) +
    stat_summary(data = geno,
                 fun.data=mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.1, size =1.5,
                 color = ifelse(z > 1,'green',
                                ifelse(z < -1,'red',
                                       'yellow'))) +
    #text annotations
    stat_summary(data = geno,fun=mean,
                 geom="label", aes(label=round(..x..,2), y=1.2), size=6) +
    geom_text(data = geno,aes(x=round(mean(value),2), y=1.3), size = 5,
              label='Genotype Mean \u00B1 SD') +
    stat_summary(fun=mean, geom="text", aes(label=round(..x..,2), y=0.75), size=6) +
    geom_text(aes(x=round(mean(value),2), y=0.65), label= 'Population \n Mean \u00B1 SD',
              size=5, lineheight = 0.75) +
    stat_summary(fun=min, geom="text", aes(label=round(..x..,2), y=0.75), size=6) +
    geom_text(aes(x=round(min(value),2), y=0.65), label= 'Population \n Min',
              size=5, lineheight = 0.75) +
    stat_summary(fun=max, geom="text", aes(label=round(..x..,2), y=0.75),size=6) +
    geom_text(aes(x=round(max(value),2), y=0.65), label= 'Population \n Max',
              size=5, lineheight = 0.75) +
    theme_minimal() +
    theme(line = element_blank(), text = element_blank(),
          plot.margin=unit(c(0,2,0,2),"cm"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank()) + # get rid of minor grid
    coord_cartesian(clip = 'off')
}

#change genet name based on the chosen genet definition
#calculate composite indicies
genoDefinition <- function(x) {
    if(x == 'name'){
      dat <- corals
    } else if(x == 'snp'){
       dat <- corals %>%
        drop_na(mlg_id) %>%
        mutate(genotype = mlg_id)
    } else if(x == 'csr'){
      dat <- corals %>%
        drop_na(CSR_accession) %>%
        mutate(genotype = CSR_accession)
    }
  
  #create global rankings
  dat <- dat %>%
    #rank non-bleaching traits w/o any stress
    filter(!grepl('stress', filters, ignore.case = T) & !trait %in% bleachingTraits) %>%
    group_by(genotype, trait, method, unit) %>% 
    summarise(mean = mean(value, na.rm = T)) %>%
    mutate(rank = as.numeric(ntile(desc(mean),3))) %>%
    ungroup() %>%
    #rank the bleaching traits
    bind_rows(dat %>%
                filter(trait %in% bleachingTraits) %>%
                group_by(trait, method, unit, genotype) %>% 
                summarise(mean = mean(value, na.rm = T)) %>%
                mutate(rank = as.numeric(ntile(desc(mean),3))) %>%
                ungroup()) %>%
    select(-mean) %>%
    right_join(dat,
               by=c("genotype", "trait", "method", "unit")) %>%
    arrange(genotype, .before="trait") %>%
    #switch ranking for Bayes Rel Risk traits, <1 = low risk; > 1 = high risk
    mutate(rank = case_when(
      method=="Bayes Relative Risk" & rank==1 ~ 3,
      method=="Bayes Relative Risk" & rank==3 ~ 1,
      TRUE ~ rank))

  #calculate composite indicies
  #composite growth index
    compGrowthIndex <- dat  %>%
      #filter to only the traits we want for growth
      filter(trait %in% c('6-month linear growth', 'annual linear growth',
                          '6-month colony volumetric growth',
                          'annual colony volumetric growth',
                          'mass normalized daily calcification',
                          'dark calcification', 'light calcification')) %>%
      mutate(type = case_when(trait %in% c('6-month linear growth', 'annual linear growth',
                                           '6-month colony volumetric growth',
                                           'annual colony volumetric growth') ~ 'linear',
                              TRUE ~ 'calcification')) %>%
      #filter by modified z-score <=3.5
      group_by(trait, method, unit) %>%
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T))) <= 3.5) %>%
      #calculate std z score
      group_by(trait, method, unit) %>%
      mutate(value = scale(value)) %>%
      #calculate avg calcification and linear score
      group_by(genotype,type) %>%
      summarise(score = mean(value, na.rm=T)) %>%
      #compute composite growth index
      group_by(genotype) %>%
      summarise(value = mean(score, na.rm=T),
                trait = "composite growth index",
                unit = "SD",
                method = "composite index",
                datafile_name = "AcDC Composite Index",
                #adding in all filters so always shows
                filters = "ambient temperature; ambient pH; wet season; dry season; bn; ln; lab; reef; Bleaching Stress; OA Stress; Stress Hardening") %>%
      ungroup() %>%
      mutate(rank = ntile(desc(value),3))

  #composite bleaching resistance index
    compBleachingIndex <- dat %>%
      #filter to only the traits we want for growth
      filter(trait %in% c('bleaching photochemical efficiency',
                          'bleaching color score',
                          'bleaching R-score',
                          'ED50')) %>%
      #filter by modified z-score <=3.5
      group_by(trait, method, unit) %>%
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T))) <= 3.5) %>%
      #calculate std z score
      group_by(trait, method, unit) %>%
      mutate(value = scale(value)) %>%
      #compute composite growth index
      group_by(genotype) %>%
      summarise(value = mean(value, na.rm=T),
                trait = "composite bleaching resistance index",
                unit = "SD",
                method = "composite index",
                datafile_name = "AcDC Composite Index",
                #adding in all filters so always shows
                filters = "ambient temperature; ambient pH; wet season; dry season; bn; ln; lab; reef; Bleaching Stress; OA Stress; Stress Hardening") %>%
      ungroup() %>%
      mutate(rank = ntile(desc(value),3))

    #return this df w/ correct genotype definition + composite indicies
    dat <- bind_rows(dat, compGrowthIndex, compBleachingIndex)

    return(dat)
}

#UI ----
ui <- tagList(
  dashboardPage(
    dashboardHeader(title = "Acropora cervicornis Data Coordination Hub", titleWidth = 500,
                    tags$li(a(href = 'https://www.coral.noaa.gov/',
                              target = "_blank",
                              h4("NOAA AOML Coral Program",
                                 style="display:inline; vertical-align: middle;"),
                              img(src='coralProgramLogo_alt.png',
                                  title = "NOAA AOML Coral Program", height="40px"),
                              style = "padding-top:5px; padding-bottom:5px;"),
                            class = "dropdown"),
                    tags$li(a(href="https://www.noaa.gov/",
                              target = "_blank",
                              img(src = 'noaaLogo.png',
                                  title = "National Oceanic and Atmospheric Administration", height = "40px"),
                              style = "padding-top:5px; padding-bottom:5px;"),
                            class = "dropdown")
    ),
    #dashboard sidebar-----
    dashboardSidebar(
      actionLink("goToHomePage", 
                 img(src='large_logo.png', 
                     width='65%',
                     alt='AcDC Logo',
                     style="padding-bottom:5px;")
      ),
      sidebarMenu(
        id = "tabs",
        menuItem("Home", tabName = "home_tab", icon = icon("home", lib = "glyphicon")),
        menuItem("Trait Analysis", tabName = "traits_tab",
                 icon = icon("dashboard", lib = "glyphicon")),
        menuItem("Genet Report", tabName = "genotypes_tab",
                 icon = icon("equalizer", lib = "glyphicon")),
        menuItem("Genet Comparison", tabName = "genotypeCompare_tab",
                 icon = icon("balance-scale", lib = "font-awesome")),
        menuItem("Raw Data", tabName = "raw_tab",
                 icon = icon("database", lib = "font-awesome")),
        menuItem("Data Sources", tabName='dataSources',
                 icon = icon("book", lib = "glyphicon")),
        menuItem("Submit Data", tabName = "uploadData_tab",
                 icon=icon("cloud-upload-alt", lib="font-awesome")),
        menuItem("Glossary", tabName = "glossary_tab",
                 icon = icon("info-sign", lib = "glyphicon"))
      )
    ),
    #dashboard Body ----
    dashboardBody(
      #link to external css
      tags$head(
        tags$link(rel="shortcut icon", href="use.ico"),

        tags$style(type="text/css", '
            
            .btn-info {color: #fff;}
            
            .navbar-custom-menu {padding-right: 10%;}
      
            .skin-blue .main-header .navbar {background-color: #273f78;}
            .skin-blue .main-header .logo {background-color: #273f78;color: #fff; border-bottom: 0 solid transparent;}
            .skin-blue .main-header .logo:hover{background-color: #273f78}
            
            .content-wrapper, .right-side {background-color: #ffffff; z-index: 800; padding-bottom:50px;}

            body {font-family: Georgia, serif;}
            h1, h2, h3, h4, h5, h6, .sidebar-menu, table {font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;}
            
            #sidebarItemExpanded img {display: block; margin-left:auto; margin-right:auto; padding-top: 5%}
            
            .main-header {position: fixed; width:100%; max-height: 100px; z-index:3000;}
            .tab-content {padding-top: 3%;}
        
            ')
      ),
      tags$style(type="text/css",
                 ".shiny-output-error {visibility: hidden;}",
                 ".shiny-output-error:before {visibility: hidden;}",
                 ".fa.fa-angle-left.pull-right {display: none;}",
                 ".sidebar-menu .treeview-menu {padding: 0 0 0 20px;}",
                 ".sidebar {position: fixed; overflow: visible; width: 230px;, z-index: 2000;}",
                 ".left-side, .main-sidebar{z-index:2000; box-shadow: 5px 0px 5px rgba(0, 0, 0, 0.1);}",
                 ".leaflet .legend i{border-radius: 50%; width: 10px; height: 10px;
                   margin-top: 4px;}",
                 "#dataSourceOverview td:first-child {font-weight: bold; text-align: left;}",
                 "#dataSourceOverview th {display:none;}",
                 ".qtip-big {font-size: 15px; line-height: 18px; white-space: nowrap; 
                            word-spacing: 1px;}",
                 "p{font-size: large;}",
      ),
      #modal styles
      tags$style(type="text/css", 
                 ".modal-title {font-size:1.7em;}"),
      tabItems(
        #Home page-----
        tabItem(tabName = "home_tab", fluidPage(
          fluidRow(
            h1("NOAA's", em("Acropora cervicornis"), "Data Coordination hub"),
            p("Welcome to NOAA's", em("Acropora cervicornis"), "Data Coordination hub (AcDC), a decision support tool created by ",
              tags$a(href = 'https://www.coral.noaa.gov/', target = '_blank',
                     "NOAA's Atlantic Oceanographic and Meteorological Laboratory Coral Program"),
              "to assimilate disparate datasets of",
              em("A. cervicornis"), "and identify genets harboring resilient phenotypes."),
            p("The goal of this database is to serve the restoration and scientific community in their efforts to restore populations of ",
              em("A. cervicornis."), "Please explore the sections below and the tabs on the left to learn more about this project."),
            h3("Guide to Using the Database"),
            p("To use this database, please start by exploring one of the three quick
                reference tools in the menu bar:",
              actionLink("goToTraitAnalysis","Trait Analysis,"),
              actionLink("goToGenotypeReport","Genet Report,"),
              actionLink("goToGenetComparison", "or Genet Comparison."),
              "You may choose specific traits or genets to analyze and explore in depth.
                Alternatively, please dive into the data by navigating to the Raw Data tab. To better understand the design of the database and the features of this web tool, please go to the",
              actionLink("goToGlossary","Glossary tab on the left."))),
          fluidRow(
            h3("Database Overview"),
            column(3,dataTableOutput("dbOverview")),
            column(8, style="padding-left:10px;padding-top:25px;",
                   p("The table on the left summarizes the number of traits measured, genets analyzed, and datasets within this database.
                    Below you can find a map of the locations where all the data were
                    sourced from and the source locations of genets used in this database."))),
          fluidRow(
            h3("Map of Data within Database"),
            leafletOutput("dbMap"))
        )),
        #Raw Table page ----
        tabItem(tabName = "raw_tab", fluidPage(
          fluidRow(h2("Raw Table"),
                   p("Select the traits or genotypes you wish to collect all the raw data from. The following table will display all the results and corresponding contextual data from the choices you select. You may download the data by clicking the download button below.")),
          fluidRow(
            column(5,
                   h4("Select the Traits to View Raw Data"),
                   selectizeInput("rawselectTraits", multiple = T,
                                  label = "Select the traits you want to analyze:",
                                  choices = c("biomechanical", "bleaching", "calcification", "disease", "reproduction", "respirometry", "symbiont physiology",
                                              "tissue properties", "total linear extension","volumetric growth", "wound healing"),
                                  selected=NULL)),
            column(5,
                   h4("Select the Genotypes to View Raw Data"),
                   selectizeInput("rawselectGenos", multiple = T,
                                  label = "Select the genotypes you want to analyze:",
                                  choices = geno_names$geno_name, selected=NULL))),
          fluidRow(downloadButton("downloadRaw", "Download")),
          fluidRow(dataTableOutput("rawTable"))
        )),
        #Traits Table page-----
        tabItem(tabName = "traits_tab", fluidPage(
          fluidRow(h2("Analyze by Traits"),
                   p("Select the family of traits you want to explore below. You will then be able to see a suite of summary statistics and a breakdown by genet.")),
          fluidRow(
            column(5,
                   selectInput("traitSelector",
                               label = "Select the traits you want to analyze:",
                               choices = c("biomechanical", "bleaching", "disease", "growth rates",
                                           "host physiology","reproduction", "wound healing")),
                   h3("Data Filters"),
                   radioButtons("locationFilters", 
                                label = tags$b("Observation Setting Filters:",
                                               #tooltip for Location Filters
                                               infoBtn('LocationFiltPop') %>% 
                                                 bsPopover(title = "Observation Setting Filters",
                                                           content = "Filtered data will only be from the selected observation location type.",
                                                           placement = "right",
                                                           trigger = "focus")),
                                c("All Data" = "all",
                                  "Lab Data" = "lab",
                                  "Outplant Monitoring" = "field",
                                  "Nursery Data" = "nursery")),
                   radioButtons("temporalFilters", 
                                label = tags$b("Temporal Filters:",
                                               #tooltip for Temporal Filters
                                               infoBtn('temporalFiltPop') %>%
                                                 bsPopover(title = "Temporal Filters",
                                                           content = "Filter data by season. Wet season is from May 1-October 31 and dry season is from November 1-April 30. All data for a given trait must be collected within a single season for a filter to apply.",
                                                           placement = "right", 
                                                           trigger = "focus")),
                                c("All Data" = "all",
                                  "Wet Season (May-Oct)" = "wet",
                                  "Dry Season (Nov-Apr)" = "dry")),
                   radioButtons("andORFilt", 
                                label = tags$b("Experimental Filters Logic:",
                                               #tooltip for AND/OR Logic
                                               infoBtn("LogicFiltPop") %>%
                                                 bsPopover(title = "Experimental Filters Logic",
                                                           content = "Select whether Experimental Filters below apply 'AND' or 'OR' logic. AND logic filters data that corresponds to the combination of every selected filter. OR logic filters data that corresponds to at least one selected filter.",
                                                           placement = "right", 
                                                           trigger = "focus")),
                                c("OR", "AND"),
                                inline = TRUE),
                   checkboxGroupInput("trait_filters_choices", 
                                      label = tags$b("Experimental Filters:",
                                                     #tooltip for Experimental Filters
                                                     infoBtn("expFiltPop") %>%
                                                       bsPopover(title = "Experimental Filters",
                                                                 content = "Filter data by the experimental conditions. Use the Experimental Filters Logic choices above for greater control. ",
                                                                 placement = "right",
                                                                 trigger = "focus")),
                                      c("Bleaching Stress" = 'Bleaching Stress',
                                        "Ambient Temperatures" = 'ambient temperature',
                                        "OA Stress" = 'OA Stress',
                                        "Ambient pH" = 'ambient pH',
                                        "Stress Hardened" = 'Stress Hardening'),
                                      inline = TRUE,
                                      selected = c('Bleaching Stress',
                                                   'ambient temperature',
                                                   'OA Stress',
                                                   'ambient pH',
                                                   'Stress Hardening')),
                   actionButton("clearAllFilters", "Clear Filters"),
                   radioButtons("TraitAnalysisGenotypeGroupChoice", 
                                label = tags$b("Genet Grouping:",
                                               #tooltip for GenotypeGrouping
                                               infoBtn("GenotypeSelectPop") %>%
                                                 bsPopover(title = "Genet Identifier Selection",
                                                           content = "Select whether the summary statistics below define genets by local nursery name, STAGdb multilocus genotype ID, or the Coral Sample Registry accession number.",
                                                           placement = "right", 
                                                           trigger = "focus")),
                                c("Local Name" ="name","SNP"="snp", "CSR Accession ID"="csr"),
                                inline = TRUE)
                   # create clear filters button here
                   ),
            column(7,
                   h3("Filter by Map"),
                   p("Use the polygon tool to filter data within the selected area."),
                   leafletOutput("TraitsFilterMap"))),
          fluidRow(h2("Trait Summary Statistics"),
                   p("This table lists the summary statistics for the trait you 
                   selected above given the applied filters. Select the trait you 
                   want to explore in depth and the table below will provide a 
                   genet breakdown. For definitions of each trait,",
                     actionLink("goToGlossaryFROMtraits", "please explore the Glossary."))),
          fluidRow(dataTableOutput("traitSummary")),
          fluidRow(h2("Grouped by Genet"),
                   p("This table lists the summary statistics for each genet for the trait you selected above.
                           You can click the mean for each genet to compare it against other genets
                          or you can click the genet to navigate to the Genet Report page and learn more about the selected genet Additionally, you may click the datasets row to view the metadata of the data sources. You may save a copy of this table by clicking the download button directly below the table.")),
          fluidRow(dataTableOutput("traitStats")),
          fluidRow(downloadButton("downloadTraits", "Download")),
        )),
        #Genotypes Tables page-----
        tabItem(tabName = "genotypes_tab", fluidPage(
          useShinyjs(),
          fluidRow(h2("Analyze by Genet"),
                   p("Select a genet of interest. You will be able to view the corresponding metadata, locations of measured traits, and a report card of the chosen genet's performance as measured against the population. 
                     If you would like to compare data from multiple genets,",
                     actionLink("goToGenetComparisonFROMgenotypes", "please go to the 'Genet Comparison' page."))),
          fluidRow(
            h3("Select a Genet to Analyze:")),
          fluidRow(
            column(3,
                   radioButtons("genotypeGroupChoice",
                                label = "Select how you want to choose a genet:",
                                choices = c("Local Name" ="name","SNP"="snp", "CSR Accession ID"="csr"))
            ),
            column(4,
                   selectInput('selectGenotypes', label = 'Select the genet you want to analyze:',
                               choices = c()),
                   style="z-index:1002;"),
          ),
          fluidRow(h2("Genet overview"),
                   column(3,
                          dataTableOutput("genotypeOverivew")),
                   column(9,
                          leafletOutput("genoMap")
                   )),
          fluidRow(h2("Genet Report Card"),
                   p("Below you will find a table that lists the summary statistics for all measured
                       traits available for the genet you selected above. You can click the mean for each genet to compare it against other genets and click the datasets row to view the metadata of the data source. Additionally, you may save a copy of this table by clicking the download button directly below the table.")),
          fluidRow(dataTableOutput("genotypeTable")),
          fluidRow(downloadButton("downloadGenotype", "Download"))
        )),
        #Generate Genotype Report --------
        tabItem(tabName = "genotypeCompare_tab", fluidPage(
          useShinyjs(),
          fluidRow(h2("Genet Comparison"),
                   p("Select the genets you would like to include in your comparison report. You will then be able to view a pooled report card of every trait for each genet you chose.")),
          fluidRow(
            h3("Select the Genets to Build Your Report:"),
            p("If you are interested in comparing a large number of genets, consider using the 'Genet Source Map' tool and drawing a polygon around the geographic area of interest.")),
          fluidRow(
            column(3,
                   radioButtons("genotypeTypeReportChoice",
                                label = "Select how you want to choose a genet:",
                                choices = c("Local Name" ="name","SNP"="snp",
                                            "CSR Accession ID"="csr", "Genet Source Map"="map"))
            ),
            column(9,
                   uiOutput('genotypesReportChoice'),
                   style="z-index:1002;"),
          ),
          fluidRow(h2("Genet Report Card"),
                   p("Below you will find a table that lists the summary statistics for all measured
                       traits available for the genotypes you selected above. You can click the mean for each genet to compare it against other genotypes and click the datasets row to view the metadata of the data source. Additionally, you may save a copy of this table by clicking the download button directly below the table.")),
          fluidRow(dataTableOutput("genotypeReport")),
          fluidRow(downloadButton("downloadGenotypeReport", "Download"))
        )),
        #DataSources Page -----
        tabItem(tabName = "dataSources", fluidPage(
          fluidRow(h2("Data Sources"),
                   p("Below you can select a data source to view its metadata including
                         the title, name and email of the PI, and a short description of the
                         data source, methods and procedures. First select a genet or a trait
                        you want to view datasources of. Then use the data source selector to choose a specific dataset you want to view the metadata of.")),
          fluidRow(
            column(5,selectInput("traits_dataSources",
                                 label = "Select the traits you want to view data sources of:",
                                 choices = c("","biomechanical", "bleaching", "disease", "growth rates", "host physiology",
                                             "reproduction", "wound healing"),
                                 selected = "")),
            column(5,selectInput('genotypes_dataSources',
                                 label = 'Select the genotypes you want to view data sources of:',
                                 choices = c("", geno_names$geno_name), selected = ""))
          ),
          fluidRow(
            column(5,offset = 3,
                   selectInput('selectDataSource',
                               label = 'Select the data source you want to view metadata of:',
                               choices = c("",datasets$title)))
          ),
          fluidRow(
            uiOutput("dataSource_header"),
            column(8,tableOutput("dataSourceOverview"))
          ),
          fluidRow(
            column(4, offset=1,
                   uiOutput("downloadDataSource_output")
            )),
        )),
        #Data Submission Page ----
        tabItem(tabName = "uploadData_tab", fluidPage(
          fluidRow(h2("Submit Data to AcDC")),
          p("To facilitate the submission of data into the database, we have designed data and metadata submission templates.
          The use of these templates is imperative to the QA/QC process of the database and streamlines the data submission process."),
          p("Please click the download button below to download a zip folder containing all of the necessary submission documents as well as a guide to assist you
          with the data submission process. Briefly, the database follows an observation-measurement ontology and data submission
          must follow this outline by placing each unique observation on one row with the corresponding measurements as unique columns
          for that row. The metadata sheet captures the information about the measured traits including methods and standards as well as information about
          the genotypes analyzed in the dataset you are submitting."),
          p("For further questions or to submit datasets, please email",
            a(href="mailto:patrick.kiel@noaa.gov?subject=AcDC%20Data%20Submission",
              "the Coral Program at NOAA-AOML.")),
          fluidRow(
            column(4, offset=1,
                   uiOutput("acdcSubmissionDownload_output")
            ))
        )),
        #Glossary page-----
        tabItem(tabName = "glossary_tab", fluidPage(
          h2("About"),
          includeHTML("www/glossary.html")
        ))
      )),
    
  ), #end dashboardPage
  #dashboard Footer ----
  tags$footer(tags$a(href="https://www.noaa.gov/protecting-your-privacy","Protecting Your Privacy | ", target="_blank"),
              tags$a(href="https://www.noaa.gov/foia-freedom-of-information-act","FOIA | ", target="_blank"),
              tags$a(href="https://www.cio.noaa.gov/services_programs/info_quality.html","Information Quality | ", target="_blank"),
              tags$a(href="https://www.noaa.gov/disclaimer","Disclaimer | ", target="_blank"),
              tags$a(href="https://www.usa.gov/","USA.gov | ", target="_blank"),
              tags$a(href="https://www.ready.gov/","Ready.gov | ", target="_blank"),
              tags$a(href="mailto:webmaster@coral.aoml.noaa.gov?subject=AcDC%20Website","Contact Webmaster", target="_blank"),
              tags$a(href = 'https://www.aoml.noaa.gov/', target = "_blank",
                     img(src='aomlLogo.png',
                         title = "NOAA AOML", height="45px",
                         style="padding-bottom: 5px;"),
                     style = "padding-left: 70px;"),
              tags$a(href = 'https://www.noaa.gov/', target = "_blank",
                     img(src = 'noaaLogo.png',
                         title = "NOAA", height = "35px",
                         style="padding-bottom: 5px;")),
              tags$a(href = 'https://www.commerce.gov/', target = "_blank",
                     img(src = 'docLogo.gif',
                         title = "Department of Comerce", height = "35px",
                         style="padding-bottom: 5px;")),
              align = "center", 
              style = "
              position:fixed;
              bottom:0;
              width:100%;
              height:50px;   /* Height of the footer */
              color: white;
              padding: 10px;
              background-color: black;
              z-index: 1500;")
)#end tagList




#server ----
server <- function(input, output, session) {
  origin <- dbGetQuery(db, "SELECT geno_name AS genotype, origin_nursery AS 'origin nursery'
                              FROM genotypes") 
  
  #global datatable options display NA's in clean format instead of blank spaces----
  rowCallback <- c(
    "function(row, data){",
    "  for(var i=0; i<data.length; i++){",
    "    if(data[i] === null){",
    "      $('td:eq('+i+')', row).html('NA')",
    "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
    "    }",
    "  }",
    "}"  
  )
  
  #global datatable options----
  options(DT.options = list(
    pageLength = 10,
    rowCallback = JS(rowCallback)))
  
  #reactive value to change selection of selectGenotypes
  changingGenotypeSelected <- reactiveVal("")
  
  #reactive value to change selection of selectTraits
  changingTraitSelected <- reactiveVal("")
  
  #reactive value to change selection of selectTraits w/ the correct method
  changingTraitSelectedMethod <- reactiveVal("")
  
  # Direct Link to Genotype Page w/ Query Strings ------
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if(!is.null(query$type) & !is.null(query$id)) {
      #limit query$type to only applicable options
      if(!query$type %in% c('name','snp','csr')) return()
      
      # #if the id entered does not match db, give a warning
      if(!query$id %in% geno_names$geno_name[!is.na(geno_names$geno_name)] & query$type == 'name') {
        showModal(modalDialog(
          title = 'Incorrect Genet Name',
          "The genet you are attempting to access does not exist in our database. Please ensure you have properly entered the genet's name. If the name is correct and this error persists,
              please contact the AcDC team and we will be happy to help ameliorate this issue. In the meantime, please select an available genet to view its physiologal data.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else if(!query$id %in% geno_names$mlg_id[!is.na(geno_names$mlg_id)] & query$type=='snp') {
        showModal(modalDialog(
          title = 'Incorrect STAGdb MLG ID',
          "The genet you are attempting to access does not exist in our database. Please ensure you have properly entered the genet's MLG ID. If the id is correct and this error persists,
              please contact the AcDC team and we will be happy to help ameliorate this issue. In the meantime, please select an available genet to view its physiologal data.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else if(!query$id %in% geno_names$CSR_accession[!is.na(geno_names$CSR_accession)] & query$type=='csr') {
        showModal(modalDialog(
          title = 'Incorrect CSR Accession Record',
          "The genet you are attempting to access does not exist in our database. Please ensure you have properly entered the genet's CSR accession record. If the accession record is correct and this error persists,
              please contact the AcDC team and we will be happy to help ameliorate this issue. In the meantime, please select an available genet to view its physiologal data.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
      
      #give the reactVal the selection
      changingGenotypeSelected(query$id)
      
      #set the correct radio choice
      #the observer will sense the change and set the correct selectChoice as well
      updateRadioButtons(session, 'genotypeGroupChoice', selected = query$type)
      
      updateTabItems(session, 'tabs', selected = 'genotypes_tab')
    }
  })
  
  # Direct Link to Traits Page w/ Query Strings ------
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if(!is.null(query$type) & !is.null(query$id)) {
      #limit query$type to only trait
      if(query$type != 'trait') return()
      
      # #if the id entered does not match db, give a warning
      if(!query$id %in% c("biomechanical", "bleaching", "disease", "growth rates", "reproduction", "respirometry", 
                          "symbiont physiology", "tissue properties", "wound healing")) {
        showModal(modalDialog(
          title = 'Incorrect Trait Family',
          "The trait family you are attempting to access does not exist in our database. Please ensure you have entered the trait correctly. If the name is correct and this error persists,
              please contact the AcDC team and we will be happy to help ameliorate this issue. In the meantime, please select an available trait family to view its physiologal data.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else {
        updateSelectInput(session, 'traitSelector',
                          selected = query$id)
      }
      
      updateTabItems(session, 'tabs', selected = 'traits_tab')
    }
  })
  
  
  # Home Page db Overview ---------------------------------------------------
  #give a warning modal on page load
  showModal(modalDialog(
    title = 'Acropora cervicornis Data Coordination Hub',
    p("This application is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce. All NOAA data are provided on an \'as is\' basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this app will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a Department of Commerce bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the Department of Commerce or the United States Government."),
    footer =  modalButton("Confirm")
  ))
  
  #dbOverview table
  #num genotypes, num traits, num datasets,
  overview <- dbGetQuery(db, "Select DISTINCT 
                                    trait_name,
                                    m.geno_id,
                                    datafile_name,
                                    location_name,
                                    location_lat,
                                    location_long
                                FROM calculatedmeasurements AS m
                                INNER JOIN traits AS t
                                ON t.trait_id = m.trait_id
                                INNER JOIN locations AS l
                                ON l.location_id = m.location_id
                                INNER JOIN datasets AS d
                                ON d.dataset_id = m.dataset_id")
  
  geno_source <- dbGetQuery(db, "Select geno_name,
                                          source_lat,
                                          source_long,
                                          source_date 
                                  FROM genotypes
                                  WHERE source_long IS NOT NULL")
  
  output$dbOverview <- renderDataTable({
    t <- data.frame(Statistic = str_to_title(c('unique genotypes', 'unique traits', 'data sources')),
                    Number = as.character(c(overview %>% select(geno_id) %>% unique() %>% nrow(),
                                            overview %>% select(trait_name) %>% unique() %>% nrow(),
                                            overview %>% select(datafile_name) %>% unique() %>% nrow())))
    
    DT::datatable(t, rownames = F, selection = 'none',
                  options = list(dom = 't', ordering = F), escape = FALSE)
  })
  
  output$dbMap <- renderLeaflet({
    points <- overview %>%
      select(name = location_name, lat = location_lat, long = location_long, geno_id, trait_name) %>%
      drop_na(lat) %>%
      group_by(lat, long) %>%
      summarise(name = name,
                genos = length(unique(geno_id)),
                traits = paste(unique(trait_name), collapse = ", "), .groups = "drop") %>%
      distinct()
    
    leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = points, group = "Observations", lng = ~long, lat = ~lat,
                       popup = ~paste0("<strong>",name,"</strong><br>",
                                       "# of genotypes: ", genos,
                                       "<br>traits: ",traits),
                       label =  ~as.character(name), stroke = FALSE, radius = 6, 
                       color = 'green', fillOpacity = 0.7) %>%
      addCircleMarkers(data = geno_source, group = "Genet Source Locations",
                       lng = ~source_long, lat = ~source_lat, label = ~as.character(geno_name),
                       popup = ~paste0("<strong>",geno_name,":</strong><br>",
                                       "Collection Date: ", ifelse(!is.na(source_date),
                                                                   as.character(as.Date(as.POSIXct(source_date, origin="1970-01-01", tz="GMT"))),
                                                                   "UNKNOWN")),
                       stroke = FALSE, radius = 6, color = 'blue', fillOpacity = 0.7) %>%
      addLayersControl(overlayGroups = c("Observations", "Genet Source Locations"),
                       options = layersControlOptions(collapsed = FALSE)) %>%
      addLegend(colors = c("green", "red"),
                labels = c("Observation", "Genet Source Location"), opacity = 0.5,
                position = "bottomright") %>%
      addScaleBar(position = 'bottomleft',
                  options = scaleBarOptions(imperial = F))
  })
  
  #change legend based on layers
  observe({
    proxy <- leafletProxy("dbMap")
    selectedLayers <- req(input$dbMap_groups)
    
    selectedColors <- character(0)
    selectedLabels <- character(0)
    if('Observations' %in% selectedLayers) {selectedColors = append(selectedColors,'green')
    selectedLabels = append(selectedLabels,'Observation')}
    if('Genet Source Locations' %in% selectedLayers) {selectedColors = append(selectedColors,'blue')
    selectedLabels = append(selectedLabels,'Genet Source Location')}
    
    # Remove any existing legend, create a new one.
    proxy %>% clearControls()
    proxy %>% addLegend(position = "bottomright",
                        colors = selectedColors, labels = selectedLabels,
                        opacity = 0.5)
  })
  
  
  # Raw Data Tab ------------------------------------------------------------
  #user friendly df for readability in wide
  WIDEgeneratedDF <- reactive({
    traits <- character(0)
    if('biomechanical' %in% input$rawselectTraits) {traits = append(traits, c('bulk density'))}
    if('bleaching' %in% input$rawselectTraits) {traits = append(traits, c('bleaching R-score', 'color score', 'photochemical efficiency', 'ED50'))}
    if('calcification' %in% input$rawselectTraits) {traits = append(traits,c('mass', 'light calcification', 'dark calcification'))}
    if('disease' %in% input$rawselectTraits) {traits = append(traits,c('disease control', 'disease exposed', 'disease death','control death', 'low credible interval', 'median relative risk', 'high credible interval'))}
    if('reproduction' %in% input$rawselectTraits) {traits = append(traits,c("colony volume", "colony height", "colony width","colony length","polyp density","oocyte density",                     
                                                                            "oocyte length","oocyte width","oocyte volume","bundle sperm density","bundle egg density"))}
    if('respirometry' %in% input$rawselectTraits) {traits = append(traits, c('P:R','respiration', 'photosynthesis'))}
    if('symbiont physiology' %in% input$rawselectTraits) {traits = append(traits, c('symbiont chlorophyll-a density', 'chlorophyll-a density', 'symbiont density'))}
    if('tissue properties' %in% input$rawselectTraits) {traits = append(traits, c('lipid density', 'tissue dry weight'))}
    if('total linear extension' %in% input$rawselectTraits) {traits = append(traits,c('TLE','branches', 'apical branches'))}
    if('volumetric growth' %in% input$rawselectTraits) {traits = append(traits,c('interstitial space volume', 'colony volume'))}
    if('wound healing' %in% input$rawselectTraits) {traits = append(traits, c('lesion healing', 'lesion area'))}
    
    #set up paramatized query for measurements table
    data_sql <- glue_sql("Select 
                                geno_name AS genotype,
                                m.observation_id AS observation_id,
                                location_name AS location,
                                location_lat AS latitude,
                                location_long AS longitude,
                                depth,
                                type,
                                trait_name AS trait,
                            	  value,
                            	  unit,
                            	  method,
                            	  datafile_name
                              FROM
	                            measurements AS m
                            	INNER JOIN
                            	observations AS o
                            	ON m.observation_id = o.observation_id
                                INNER JOIN
                                standards AS s
                                ON m.standard_id = s.standard_id
                                INNER JOIN
                                genotypes AS g
                                ON m.geno_id = g.geno_id
                                INNER JOIN
                                traits AS t
                                ON m.trait_id = t.trait_id
                                INNER JOIN
                                locations AS l
                                ON o.location_id = l.location_id
                                INNER JOIN
                                methods as meth
                                ON m.method_id = meth.method_id
                                INNER JOIN
                                datasets AS d
                                ON d.dataset_id = o.dataset_id
                                WHERE (trait_name IN ({traits*}) OR trait_class = 'contextual')
                                AND d.private != 1
              									AND m.observation_id IN
              										(SELECT observation_id 
              										FROM measurements m
              										INNER JOIN traits t
              										ON m.trait_id=t.trait_id
              										AND t.trait_name IN ({traits*}))
                                ORDER BY m.observation_id",
                         traits = traits,
                         .con = db)
    data <- dbSendQuery(db, data_sql)
    traitsTemp <- dbFetch(data)
    
    #set up paramatized query for select by genotypes in calculated measurements table
    data_sql3 <- glue_sql("Select 
                               geno_name AS genotype,
                                m.observation_id AS observation_id,
                                location_name AS location,
                                location_lat AS latitude,
                                location_long AS longitude,
                                depth,
                                type,
                                trait_name AS trait,
                            	  value,
                            	  unit,
                            	  method,
                            	  datafile_name
                              FROM
	                            measurements AS m
                            	INNER JOIN
                            	observations AS o
                            	ON m.observation_id = o.observation_id
                                INNER JOIN
                                standards AS s
                                ON m.standard_id = s.standard_id
                                INNER JOIN
                                genotypes AS g
                                ON m.geno_id = g.geno_id
                                INNER JOIN
                                traits AS t
                                ON m.trait_id = t.trait_id
                                INNER JOIN
                                locations AS l
                                ON o.location_id = l.location_id
                                INNER JOIN
                                methods as meth
                                ON m.method_id = meth.method_id
                                INNER JOIN
                                datasets AS d
                                ON d.dataset_id = o.dataset_id
                                WHERE geno_name IN ({genos*})
                                AND d.private != 1
                                ORDER BY m.observation_id",
                          genos = input$rawselectGenos,
                          .con = db)
    data3 <- dbSendQuery(db, data_sql3)
    genosTemp <- dbFetch(data3)
    
    temp <- bind_rows(traitsTemp, genosTemp) %>% distinct()
    
    #capture when no public data are available
    if(nrow(temp)==0 & !is.null(input$rawselectTraits)) {
      showModal(modalDialog(
        title = 'Your selections have no matching public data',
        p("Your applied filters to the raw data returned no matching public records. Data may be in our database, however, the authors have asked for the data to remain private."),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if(nrow(temp)==0 & !is.null(input$rawselectGenos)) {
      showModal(modalDialog(
        title = 'Your selections have no matching public data',
        p("Your applied filters to the raw data returned no matching public records. Data may be in our database, however, the authors have asked for the data to remain private."),
        easyClose = TRUE,
        footer = NULL
      ))
    } else {
      #cast wide
      prefix <- unique(temp$trait) #used to help reorder columns
      temp <- temp %>%
        pivot_wider(names_from = trait,
                    values_from = c(value, unit, method),
                    names_glue = "{trait}_{.value}")
      #reorder columns 
      names_to_order <- map(prefix, ~ names(temp)[grep(paste0(.x,"_"), names(temp))]) %>% unlist
      names_id <- setdiff(names(temp), names_to_order)
      #return reordered columns
      temp %>%
        select(all_of(names_id), all_of(names_to_order)) %>%
        select(-observation_id) %>%
        #remove empty columns
        Filter(function(y)!all(is.na(y)), .)
    }
  })
  
  output$rawTable <- renderDataTable({
    DT::datatable(WIDEgeneratedDF()%>%
                    rename_with(str_to_title),
                  rownames = F,
                  extensions = "FixedColumns",
                  options = list(
                    scrollX = TRUE,
                    fixedColumns = list(leftColumns = 1)
                  ))
  })
  
  #create .csv file with contact + citation information
  citations <- reactive({ 
    df <- WIDEgeneratedDF()
    
    datasets %>%
      filter(datafile_name %in%
               (df %>% pull(datafile_name) %>% unique())) %>%
      select(title, author, email, DOI) %>%
      mutate(DOI = str_extract(DOI, "(?<=>)(.*?)(?=<)"))
  })
  
  output$downloadRaw <- downloadHandler(
    filename = "ACDC_dataExport.zip",
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      #create csv that will be in the zip folder
      fileNames <- c("dataExport.csv","citationExport.csv")
      write.csv(WIDEgeneratedDF(), "dataExport.csv", row.names = F)
      write.csv(citations(), "citationExport.csv", row.names = F)
      
      #create the zip file
      #zip(file,fileNames)
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      showModal(modalDialog(
        title = 'Warning: You must Assume the Responsibility for Downloaded Data',
        p("The data you just downloaded are the sole ownership of the authors who provided the data to us. The data are freely usable, but must be properly cited in any reports, publications, presentations or other relevant products. A csv file containing the author's contact information has also been downloaded for you to assist in proper citation. Further, the National Oceanographic and Atmospheric Administration makes no claims as to the quality assurance of the data presented here. All data are provided as a scientific product on an \'as is\' basis and the user assumes responsibility for its use."),
        footer =  modalButton("Confirm")
      ))
    }
  )
  
  # Traits Analysis --------------------------------------------------------------
  
  #clicked the go to glossary Link
  observeEvent(input$goToGlossary, {
    # do nothing if not clicked yet
    if (input$goToGlossary==0) return()
    updateTabItems(session, 'tabs', selected = 'glossary_tab')
  })
  
  observeEvent(input$goToGlossaryFROMtraits, {
    # do nothing if not clicked yet
    if (input$goToGlossaryFROMtraits==0) return()
    updateTabItems(session, 'tabs', selected = 'glossary_tab')
  })
  
  observeEvent(input$goToHomePage, {
    #do nothing if not clicked yet
    if (input$goToHomePage==0) return()
    updateTabItems(session, 'tabs', selected = 'home_tab')
  })
  
  observeEvent(input$goToGenetComparison, {
    #do nothing if not clicked yet
    if(input$goToGenetComparison==0) return()
    updateTabItems(session,'tabs', selected = 'genotypeCompare_tab')
  })
  
  observeEvent(input$goToGenetComparisonFROMgenotypes, {
    #do nothing if not clicked yet
    if(input$goToGenetComparisonFROMgenotypes==0) return()
    updateTabItems(session,'tabs', selected = 'genotypeCompare_tab')
  })
  
  observeEvent(input$goToTraitAnalysis, {
    #do nothing if not clicked yet
    if(input$goToTraitAnalysis==0) return()
    updateTabItems(session,'tabs', selected = 'traits_tab')
  })
  
  observeEvent(input$goToGenotypeReport, {
    #do nothing if not clicked yet
    if(input$goToGenotypeReport==0) return()
    updateTabItems(session,'tabs', selected = 'genotypes_tab')
  })
  
  #adjust trait filters
  #reactive df for traits
  totalTraitsDF <- reactive({ 
    traits <- character(0)
    if('biomechanical' %in% input$traitSelector) {traits = append(traits, biomechanicalTraits)}
    if('bleaching' %in% input$traitSelector) {traits = append(traits, bleachingTraits)}
    if('disease' %in% input$traitSelector) {traits = append(traits, diseaseTraits)}
    if('growth rates' %in% input$traitSelector) {traits = append(traits,growthTraits)}
    if('host physiology' %in% input$traitSelector) {traits = append(traits,hostPhysiologyTraits)}
    if('reproduction' %in% input$traitSelector) {traits = append(traits, reproductionTraits)}
    if('wound healing' %in% input$traitSelector) {traits = append(traits, woundHealingTraits)}

    df <- genoDefinition(input$TraitAnalysisGenotypeGroupChoice) %>%
          filter(trait %in% traits)
    
    #location filters
    if(input$locationFilters=='field'){
      df <- df %>% filter(grepl('reef', filters))
    }
    if(input$locationFilters=='lab'){
      df <- df %>% filter(grepl('lab', filters))
    }
    if(input$locationFilters=='nursery'){
      df <- df %>% filter(grepl('ln', filters) | grepl('bn', filters))
    }
    #temporaral filters
    if(input$temporalFilters=='wet'){
      df <- df %>% filter(grepl('wet season', filters))
    }
    if(input$temporalFilters=='dry'){
      df <- df %>% filter(grepl('dry season', filters))
    }
    
    #filtering logic
    #all filters unchecked empty df
    if(length(input$trait_filters_choices)==0){
      df <- df[0,]
      #OR logic
    } else if(input$andORFilt == 'OR') {
      df <- df[grepl(paste(input$trait_filters_choices, collapse = "|"), df$filters, ignore.case = T),]
      #AND logic
    } else if(input$andORFilt == 'AND') {
      df <- df %>%
        filter(grepl(paste(sprintf("(?=.*%s)", input$trait_filters_choices), collapse=""), df$filters, perl=TRUE, ignore.case = T))
    }
    
    #recalculate composite indicies w/ applied filters
    #composite growth index
    if('growth rates' %in% input$traitSelector) {
      compGrowthIndex <- df  %>%
        #filter to only the traits we want for growth
        filter(trait %in% c('6-month linear growth', 'annual linear growth',
                            '6-month colony volumetric growth',
                            'annual colony volumetric growth',
                            'mass normalized daily calcification',
                            'dark calcification', 'light calcification')) %>%
        mutate(type = case_when(trait %in% c('6-month linear growth', 'annual linear growth',
                                             '6-month colony volumetric growth',
                                             'annual colony volumetric growth') ~ 'linear',
                                TRUE ~ 'calcification')) %>%
        #filter by modified z-score <=3.5
        group_by(trait, method, unit) %>%
        filter(abs(scale(value,
                         center = median(value, na.rm = T),
                         scale = mad(value, na.rm=T))) <= 3.5) %>%
        #calculate std z score
        group_by(trait, method, unit) %>%
        mutate(value = scale(value)) %>%
        #calculate avg calcification and linear score
        group_by(genotype,type) %>%
        summarise(score = mean(value, na.rm=T)) %>%
        #compute composite growth index
        group_by(genotype) %>%
        summarise(value = mean(score, na.rm=T),
                  trait = "composite growth index") %>%
        ungroup()
      
      compGrowthIndex <- df %>%
        #select only the composite growth index, replace it's value w/ recalculated
        filter(trait == "composite growth index") %>%
        select(-value) %>%
        left_join(compGrowthIndex,
                  by=c("genotype","trait"))
      
      #add back in all the other filtered data
      df <- df %>%
        filter(trait != "composite growth index") %>%
        bind_rows(compGrowthIndex)
    }
    #composite bleaching resistance index
    if('bleaching' %in% input$traitSelector) {
      compBleachingIndex <- df %>%
        #filter to only the traits we want for growth
        filter(trait %in% c('bleaching photochemical efficiency',
                            'bleaching color score',
                            'bleaching R-score',
                            'ED50')) %>%
        #filter by modified z-score <=3.5
        group_by(trait, method, unit) %>%
        filter(abs(scale(value,
                         center = median(value, na.rm = T),
                         scale = mad(value, na.rm=T))) <= 3.5) %>%
        #calculate std z score
        group_by(trait, method, unit) %>%
        mutate(value = scale(value)) %>%
        #compute composite growth index
        group_by(genotype) %>%
        summarise(value = mean(value, na.rm=T),
                  trait = "composite bleaching resistance index") %>%
        ungroup()
      
      compBleachingIndex <- df %>%
        #select only the composite bleaching index, replace it's value w/ recalculated
        filter(trait == "composite bleaching resistance index") %>%
        select(-value) %>%
        left_join(compBleachingIndex,
                  by=c("genotype","trait"))
      #add back in all the other filtered data
      df <- df %>% 
        filter(trait != "composite bleaching resistance index") %>%
        bind_rows(compBleachingIndex)
    }
    
    df %>% drop_na(value)
  })
  
  #reset trait filters
  observeEvent(input$clearAllFilters, {
    #do nothing if not clicked yet
    if(input$clearAllFilters==0) return()
    updateRadioButtons(session, "locationFilters",
                       selected = "all")
    updateRadioButtons(session, "temporalFilters",
                       selected = "all")
    updateRadioButtons(session, "andORFilt",
                       selected = "OR")
    updateCheckboxGroupInput(session, "trait_filters_choices",
                             selected = c('Bleaching Stress',
                                       'ambient temperature',
                                       'OA Stress',
                                       'ambient pH',
                                       'Stress Hardening'))
    #clear drawn polygon on map
    leafletProxy('TraitsFilterMap') %>%
      removeDrawToolbar(clearFeatures=TRUE) %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = drawPolygonOptions(shapeOptions=drawShapeOptions(fillOpacity = .4,
                                                                                       color = 'black',
                                                                                       fillColor = '#fff',
                                                                                       weight = 3)),
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE)
    
    drawnValues$activePolygon <- FALSE
  })
  
  #filter traits map
  output$TraitsFilterMap <- renderLeaflet({
    points <- totalTraitsDF() %>%
      select(name = location_name, lat, long) %>%
      distinct() %>%
      drop_na(lat)
    
    leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = points, group = "Observations", lng = ~long, lat = ~lat,
                       stroke = FALSE, radius = 6, 
                       color = 'green', fillOpacity = 0.7) %>%
      addLegend(colors = c("green"),
                labels = c("Observation"), opacity = 0.5,
                position = "bottomright") %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = drawPolygonOptions(shapeOptions=drawShapeOptions(fillOpacity = .4,
                                                                                       color = 'black',
                                                                                       fillColor = '#fff',
                                                                                       weight = 3)),
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE) %>%
      addScaleBar(position = 'bottomleft',
                  options = scaleBarOptions(imperial = F))
  })
  
  drawnValues <- reactiveValues(activePolygon=FALSE)
  
  #prevent adding polygons after the first is created
  observeEvent(input$TraitsFilterMap_draw_new_feature,{
    drawnValues$activePolygon <- TRUE
    leafletProxy('TraitsFilterMap') %>%
      removeDrawToolbar() %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = FALSE,
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE)
  })
  
  #add original toolbar back in when polygon is deleted
  observeEvent(input$TraitsFilterMap_draw_deleted_features,{
    drawnValues$activePolygon <- FALSE
    leafletProxy('TraitsFilterMap') %>%
      removeDrawToolbar() %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = drawPolygonOptions(shapeOptions=drawShapeOptions(fillOpacity = .4,
                                                                                       color = 'black',
                                                                                       fillColor = '#fff',
                                                                                       weight = 3)),
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE)
  })
  
  #filter the traitsDF by the polygon
  PolyFilteredDF <- eventReactive(input$TraitsFilterMap_draw_new_feature,{
    poly <- input$TraitsFilterMap_draw_new_feature$geometry$coordinates[[1]]
    bounds <- data.frame(long = map_dbl(poly, `[[`,1),
                         lat = map_dbl(poly, `[[`, 2))
    df <- totalTraitsDF() %>%
      mutate(inside = point.in.polygon(long, lat, bounds$long, bounds$lat)) %>%
      filter(inside >0) %>%
      select(-inside)
    
    #recalculate composite indicies based on geographic filter
    #composite growth index
    if('growth rates' %in% input$traitSelector) {
      compGrowthIndex <- df  %>%
        #filter to only the traits we want for growth
        filter(trait %in% c('6-month linear growth', 'annual linear growth',
                            '6-month colony volumetric growth',
                            'annual colony volumetric growth',
                            'mass normalized daily calcification',
                            'dark calcification', 'light calcification')) %>%
        mutate(type = case_when(trait %in% c('6-month linear growth', 'annual linear growth',
                                             '6-month colony volumetric growth',
                                             'annual colony volumetric growth') ~ 'linear',
                                TRUE ~ 'calcification')) %>%
        #filter by modified z-score <=3.5
        group_by(trait, method, unit) %>%
        filter(abs(scale(value,
                         center = median(value, na.rm = T),
                         scale = mad(value, na.rm=T))) <= 3.5) %>%
        #calculate std z score
        group_by(trait, method, unit) %>%
        mutate(value = scale(value)) %>%
        #calculate avg calcification and linear score
        group_by(genotype,type) %>%
        summarise(score = mean(value, na.rm=T)) %>%
        #compute composite growth index
        group_by(genotype) %>%
        summarise(value = mean(score, na.rm=T),
                  trait = "composite growth index") %>%
        ungroup()
      
      compGrowthIndex <- df %>%
        #select only the composite growth index, replace it's value w/ recalculated
        filter(trait == "composite growth index") %>%
        select(-value) %>%
        left_join(compGrowthIndex,
                  by=c("genotype","trait"))
      
        #add back in all the other filtered data
      df <- df %>%
        filter(trait != "composite growth index") %>%
        bind_rows(compGrowthIndex)
    }
    #composite bleaching resistance index
    if('bleaching' %in% input$traitSelector) {
      compBleachingIndex <- df %>%
        #filter to only the traits we want for growth
        filter(trait %in% c('bleaching photochemical efficiency',
                            'bleaching color score',
                            'bleaching R-score',
                            'ED50')) %>%
        #filter by modified z-score <=3.5
        group_by(trait, method, unit) %>%
        filter(abs(scale(value,
                         center = median(value, na.rm = T),
                         scale = mad(value, na.rm=T))) <= 3.5) %>%
        #calculate std z score
        group_by(trait, method, unit) %>%
        mutate(value = scale(value)) %>%
        #compute composite growth index
        group_by(genotype) %>%
        summarise(value = mean(value, na.rm=T),
                  trait = "composite bleaching resistance index") %>%
        ungroup()
      
      compBleachingIndex <- df %>%
        #select only the composite bleaching index, replace it's value w/ recalculated
        filter(trait == "composite bleaching resistance index") %>%
        select(-value) %>%
        left_join(compBleachingIndex,
                  by=c("genotype","trait"))
        #add back in all the other filtered data
        df <- df %>% 
          filter(trait != "composite bleaching resistance index") %>%
          bind_rows(compBleachingIndex)
    }
    
    df
    
  })
  
  #reset the activePolygon when trait selector is change
  observeEvent(input$traitSelector,{
    drawnValues$activePolygon <- FALSE
  }, ignoreInit = TRUE)
  
  #grab either the traitsDF or PolyFilteredDF depending on user input
  traitsDF <- reactive({
    if(drawnValues$activePolygon) {
      PolyFilteredDF()
    }else{totalTraitsDF()}
  })
  
  #create filtered choices reactive values
  filteredVals <- reactiveValues()
  
  traitSummaryDF <- reactive({
    temp <- traitsDF() %>%
      group_by(trait, method, unit) %>%
      summarise(genotypes = length(unique(genotype)),
                mean = mean(value, na.rm =T),
                sd = sd(value, na.rm = T),
                n = n(),
                min = min(value, na.rm = T),
                max = max(value, na.rm = T),
                datasets = length(unique(datafile_name)), .groups = "drop") %>%
      mutate(across(c(mean,sd,min,max), ~format(round(.x, 3), nsmall=3)),
             across(c(genotypes, mean, sd, n, min, max, datasets),~as.character(.))) %>%
      distinct() %>%
      arrange(trait, method, unit) %>%
      relocate(unit, .after = mean) %>%
      mutate(Select = if_else(row_number()==1,
                              '<input type="radio" name="%s" value="%s" checked="checked"/>',
                              '<input type="radio" name="%s" value="%s"/>'),
             .before="trait")
    
    if(changingTraitSelected()!="" & changingTraitSelected() %in% temp$trait) {
      temp <- temp %>%
        #make the selection based on selected Trait
        mutate(Select = if_else(row_number()==which(temp$trait==changingTraitSelected()),
                                '<input type="radio" name="%s" value="%s" checked="checked"/>',
                                '<input type="radio" name="%s" value="%s"/>'),
               .before="trait")
    }
    
    if(changingTraitSelected()=="" | !changingTraitSelected() %in% temp$trait) {
      filteredVals$filteredTrait <- temp$trait[1]
      filteredVals$filteredMethod  <- temp$method[1]
      filteredVals$filteredUnit <- temp$unit[1]
    } else {
      filteredVals$filteredTrait <- changingTraitSelected()
      filteredVals$filteredMethod  <- changingTraitSelectedMethod()
      filteredVals$filteredUnit <- temp$unit[which(temp$trait==changingTraitSelected())]
    }
    
    temp
  })
  
  #trait summary for selected trait
  output$traitSummary <- renderDataTable({
    
    if(changingTraitSelected()=="" | !changingTraitSelected() %in% traitSummaryDF()$trait) {
      DT::datatable(traitSummaryDF() %>%
                      rename_with(str_to_title) %>%
                      mutate(Select = if_else(row_number()==1,
                                              '<input type="radio" name="%s" value="%s" checked="checked"/>',
                                              '<input type="radio" name="%s" value="%s"/>')),
                    rownames = F, selection =  list(mode = 'single', selected = c(1)),
                    escape = FALSE,
                    options = list(dom = 't', ordering=F, pageLength = 100)) %>%
        formatStyle(11, cursor = "pointer") %>%
        formatStyle(1, textAlign="center")
    } else {
      DT::datatable(traitSummaryDF() %>%
                      rename_with(str_to_title),
                    rownames = F, selection =  list(mode = 'single', selected = c(which(traitSummaryDF()$trait==changingTraitSelected() & traitSummaryDF()$method==changingTraitSelectedMethod()))),
                    escape = FALSE,
                    options = list(dom = 't', ordering=F, pageLength = 100)) %>%
        formatStyle(11, cursor = "pointer") %>%
        formatStyle(1, textAlign="center")
    }
  })
  
  #Popup list of datasets when num datasets is clicked in trait summary tbl
  observeEvent(input$traitSummary_cell_clicked, {
    info = input$traitSummary_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 10th column
    if (is.null(info$value) || info$col != 10) return()
    showModal(modalDialog(
      title = paste("Datasets that contain the", traitSummaryDF()[info$row,2], "trait."),
      dataTableOutput("traitsModalTable"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  traitsModaldf <- reactive({
    list = traitsDF() %>%
      filter(trait==traitSummaryDF()$trait[input$traitSummary_cell_clicked$row] & method==traitSummaryDF()$method[input$traitSummary_cell_clicked$row]& unit==traitSummaryDF()$unit[input$traitSummary_cell_clicked$row]) %>%
      select(datafile_name) %>% distinct() %>%
      left_join(datasets, by = "datafile_name") %>%
      select(`Data Source Title` = title) %>%
      rename_with(str_to_title)
  })
  
  output$traitsModalTable <- renderDataTable({
    DT::datatable(traitsModaldf(),
                  rownames = T, selection = 'none',
                  options = list(dom = 'tf', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click datasource brings you to that datasource page
  observeEvent(input$traitsModalTable_cell_clicked, {
    info = input$traitsModalTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    updateTabItems(session, 'tabs', selected = 'dataSources')
    
    newChoices <- corals %>%
      select(datafile_name) %>% distinct() %>%
      left_join(datasets, by = "datafile_name") %>%
      pull(title)
    updateSelectInput(session, 'selectDataSource', choices=newChoices, selected = info$value)
    updateSelectInput(session, 'traits_dataSources', selected = "")
    updateSelectInput(session, 'genotypes_dataSources', selected = "")
    removeModal()
  })
  
  #give radio choices to the grouped by genotype table
  observeEvent(input$traitSummary_cell_clicked, {
    info = input$traitSummary_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is the 10th column
    if (is.null(info$value) || info$col==10) return()
    
    filteredVals$filteredTrait <- traitSummaryDF()$trait[info$row]
    filteredVals$filteredMethod  <- traitSummaryDF()$method[info$row]
    filteredVals$filteredUnit <- traitSummaryDF()$unit[info$row]
    
    #change the radio button
    output$traitSummary <- renderDataTable({
      DT::datatable(traitSummaryDF() %>%
                      mutate(Select = if_else(row_number()==info$row,
                                              '<input type="radio" name="%s" value="%s" checked="checked"/>',
                                              '<input type="radio" name="%s" value="%s"/>')) %>%
                      rename_with(str_to_title),
                    rownames = F, selection =  list(mode = 'single', selected = c(info$row)),
                    escape = FALSE,
                    options = list(dom = 't', ordering=F, pageLength = 100)) %>%
        formatStyle(11, cursor = "pointer") %>%
        formatStyle(1, textAlign="center")
    })
    
  })
  
  #trait stats grouped by genotype,method, unit
  traitStatsDF <- reactive({
    traitsDF() %>%
      filter(trait == filteredVals$filteredTrait & method == filteredVals$filteredMethod & unit == filteredVals$filteredUnit) %>%
      group_by(genotype, trait, method, unit, rank) %>% 
      summarise(mean = mean(value, na.rm = T),
                sd = sd(value, na.rm = T),
                n = n(),
                min = min(value, na.rm = T),
                max = max(value, na.rm = T), 
                datasets = length(unique(datafile_name)),
                .groups = "drop") %>%
      mutate(across(c(mean,sd,min,max), ~format(round(.x, 3), nsmall=3)),
             across(c(mean, sd, n, min, max, datasets), ~as.character(.))) %>%
      distinct() %>%
      arrange(desc(mean)) %>%
      relocate(rank, unit, .after = mean) %>%
      mutate(rank = ifelse(rank==1,
                           '<img src="1greenLight_small.png" height="25" />',
                           ifelse(rank==2,
                                  '<img src="2yellowLight_small.png" height="25" />',
                                  '<img src="3redLight_small.png" height="25" />')), .after='mean')
  })
  #trait summary stats table
  output$traitStats <- renderDataTable({
    globalRank_text <- tags$span(
      "Global Rank", 
      infoBtn('traitGlobalRank') %>% 
        bsPopover(title = "Global Rank",
                  content = "This column is a quick reference tool to indicate above average, average, or below average performance for a given genet for the applicable trait, method, and unit based off the total data and not the applied data filters. Green indicates the genet falls in the first tercile, yellow indicates the genet falls in the middle tercile, and red indicates the genet falls in the bottom tercile.",
                  placement = "top",
                  trigger = "hover")
    ) %>% 
      as.character()
    
    mean_text <- tags$span(
      "Mean", 
      infoBtn('traitMean') %>% 
        spsComps::bsPopover(title = "Mean",
                            content = "This mean is calculated from all data. Click the mean to view a modified box-plot to understand how the calculated mean compares to the population given the applied filters.",
                            placement = "top",
                            trigger = "hover")
    ) %>% 
      as.character()
    
    
    
    DT::datatable(traitStatsDF() %>% rename_with(str_to_title) %>%
                    #add in popovers
                    rename(!!mean_text:=Mean,
                           !!globalRank_text:=Rank),
                  rownames = F, selection = 'none', escape = F,
                  options = list(dom = 'tipl', pageLength = 1000)) %>%
      formatStyle(c(1,4,11), cursor = 'pointer')
  })
  
  #Change tab to selected genotype
  observeEvent(input$traitStats_cell_clicked, {
    info = input$traitStats_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 0) return()
    changingGenotypeSelected(info$value)
    
    updateRadioButtons(session, "genotypeGroupChoice", selected = input$TraitAnalysisGenotypeGroupChoice)
    
    updateSelectInput(session, 'selectGenotypes', selected = info$value)
    updateTabItems(session, 'tabs', selected = 'genotypes_tab')
  })
  
  #Popup list of datasets when num datasets is clicked in by genos tbl
  observeEvent(input$traitStats_cell_clicked, {
    info = input$traitStats_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 10) return()
    showModal(modalDialog(
      title = paste("Datasets that contain the", traitStatsDF()$genotype[info$row],
                    "genet and", traitStatsDF()$trait[info$row], "trait."),
      dataTableOutput("traitsModalTable2"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  traitsModalDF2 <- reactive({
    list = traitsDF() %>%
      filter(genotype == traitStatsDF()$genotype[input$traitStats_cell_clicked$row] & trait == traitStatsDF()$trait[input$traitStats_cell_clicked$row]) %>%
      select(datafile_name) %>% distinct() %>%
      left_join(datasets, by = "datafile_name") %>%
      select(`Data Source Title` = title) %>%
      rename_with(str_to_title)
  })
  
  output$traitsModalTable2 <- renderDataTable({
    DT::datatable(traitsModalDF2(), rownames = T, selection = 'none',
                  options = list(dom = 'tf', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click datasource brings you to that datasource page
  observeEvent(input$traitsModalTable2_cell_clicked, {
    info = input$traitsModalTable2_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    updateTabItems(session, 'tabs', selected = 'dataSources')
    updateSelectInput(session, 'selectDataSource', selected = info$value)
    removeModal()
  })
  
  #Modal for plot
  observeEvent(input$traitStats_cell_clicked, {
    info = input$traitStats_cell_clicked
    
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 3) return()
    showModal(modalDialog(
      id = 'plotModal',
      title = paste("Comparison of", traitStatsDF()$genotype[info$row],
                    "to the population mean of the", filteredVals$filteredTrait,
                    "trait, measured in", filteredVals$filteredUnit, "units by the", filteredVals$filteredMethod, "method."),
      span("Loading...", id="UpdateAnimate",
           style="font-size:1.5em; margin-left: auto;
                     margin-right: auto;"),
      imageOutput("trait_genoPlot"),
      easyClose = TRUE,
      footer = NULL,
    ))
  })
  
  output$trait_genoPlot <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    
    info <- input$traitStats_cell_clicked
    
    temp <- traitsDF() %>% 
      filter(trait == filteredVals$filteredTrait &
               method == filteredVals$filteredMethod &
               unit == filteredVals$filteredUnit) %>%
      #remove extreme outliers for z-score calculation
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T)))<=5) %>%
      #calculate genotype z-scores
      group_by(genotype) %>%
      summarise(mean = mean(value,na.rm=T)) %>%
      mutate(z = scale(mean)) %>%
      select(-mean) %>%
      #add data back in
      right_join(traitsDF() %>% 
                   filter(trait == filteredVals$filteredTrait &
                            method == filteredVals$filteredMethod &
                            unit == filteredVals$filteredUnit),
                 by="genotype") %>%
      #remove extreme outliers for data display
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T)))<=5)
    
    temp2 <- temp %>%
      filter(genotype == traitStatsDF()$genotype[info$row])
    z <- temp2 %>% pull(z) %>% unique()
    
    png(outfile, 
        width = 500*4, 
        height = 250*4,
        res = 72*4)
    print(modalPlot(temp,temp2,z))
    dev.off()
    shinyjs::hideElement(id='UpdateAnimate')
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 500,
         height = 250,
         alt = "This is alternate text")
  }, deleteFile = T)
  
  
  #create .csv file with contact + citation information
  traitCitations <- reactive({ 
    df <- traitsDF() %>%
      filter(trait == filteredVals$filteredTrait)
    
    datasets %>%
      filter(datafile_name %in%
               (df %>% pull(datafile_name) %>% unique())) %>%
      select(title, author, email, DOI) %>%
      mutate(DOI = str_extract(DOI, "(?<=>)(.*?)(?=<)"))
  })
  
  #download traits genotype ranking
  output$downloadTraits <- downloadHandler(
    #file name is not working for some reason
    filename = paste0("AcDC_",
                      gsub("(^|[^[:alnum:]])([[:alnum:]])",
                           "\\U\\2",
                           filteredVals$filteredTrait,
                           perl = TRUE),
                      "_genotypeSumStats.zip"),
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      #create csv that will be in the zip folder
      fileNames <- c("dataExport.csv","citationExport.csv")
      
      write.csv(traitStatsDF() %>%
                  #this is not changing for some reason
                  #mutate(unit = gsub("<\\U\\+0394>", "change", unit)) %>%
                  select(-rank),
                "dataExport.csv", row.names = F)
      
      write.csv(traitCitations(), "citationExport.csv", row.names = F)
      
      #create the zip file
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      showModal(modalDialog(
        title = 'Warning: You must Assume the Responsibility for Downloaded Data',
        p("The data you just downloaded are the summary statistics of data provided by the authors. The National Oceanographic and Atmospheric Administration makes no claims as to the quality assurance of the data presented here. All data are provided as a scientific product on an \'as is\' basis and the user assumes responsibility for its use. Please cite the authors of the original data and the AcDC website in all publications and presentations that use the downloaded data."),
        footer =  modalButton("Confirm")
      ))
    }
  )
  
  # Genotypes Report -----------------------------------------------------------
  #filter how user selects which genotype
  observeEvent(input$genotypeGroupChoice, {
    #no value to reactive value + by local name selected
    if(input$genotypeGroupChoice == 'name' & (changingGenotypeSelected()=="" || !changingGenotypeSelected() %in% geno_names$geno_name[!is.na(geno_names$geno_name)])) {
      updateSelectInput(session, 'selectGenotypes', choices=geno_names$geno_name[!is.na(geno_names$geno_name)])
      #reactive value provided + by local name selected
    } else if(input$genotypeGroupChoice == 'name' & changingGenotypeSelected()!="") {
      updateSelectInput(session, 'selectGenotypes', choices=geno_names$geno_name[!is.na(geno_names$geno_name)],
                        selected = changingGenotypeSelected())
      #no value to reactive value + by snp selected
    } else if(input$genotypeGroupChoice == 'snp' & (changingGenotypeSelected()=="" || !changingGenotypeSelected() %in% geno_names$mlg_id[!is.na(geno_names$mlg_id)])) {
      updateSelectInput(session, 'selectGenotypes',
                        choices=geno_names$mlg_id[!is.na(geno_names$mlg_id)])
      #reactive value provided + by snp selected
    } else if(input$genotypeGroupChoice == 'snp' & changingGenotypeSelected()!="") {
      updateSelectInput(session, 'selectGenotypes',
                        choices=geno_names$mlg_id[!is.na(geno_names$mlg_id)],
                        selected = changingGenotypeSelected())
      #no value to reactive value + by csr selected
    } else if(input$genotypeGroupChoice == 'csr' & (changingGenotypeSelected()=="" || !changingGenotypeSelected() %in% geno_names$CSR_accession[!is.na(geno_names$CSR_accession)])) {
      updateSelectInput(session, 'selectGenotypes',
                        choices=geno_names$CSR_accession[!is.na(geno_names$CSR_accession)])
      #reactive value provided + by csr selected
    } else if(input$genotypeGroupChoice == 'csr' & changingGenotypeSelected()!="") {
      updateSelectInput(session, 'selectGenotypes',
                        choices=geno_names$CSR_accession[!is.na(geno_names$CSR_accession)],
                        selected = changingGenotypeSelected())
    }
  })
  
  #overview table
  genotypeOverviewdf <- reactive({
    if(input$genotypeGroupChoice == 'name') {
      genoOverview_sql <- glue_sql("Select 
                                      geno_name AS genotype,
                                      alt_names,
                                      origin_nursery,
                                      source_lat,
                                      source_long,
                                      sexual_recruit,
                                      source_date,
                                      CSR_accession,
                                      mlg_id,
                                      mlg_geno_id,
                                      microsat_id
                                    FROM genotypes
                                    WHERE geno_name IN ({geno*})",
                                   geno = input$selectGenotypes,
                                   .con = db)
    } else if(input$genotypeGroupChoice == 'snp') {
      genoOverview_sql <- glue_sql("Select 
                                      geno_name AS genotype,
                                      alt_names,
                                      origin_nursery,
                                      source_lat,
                                      source_long,
                                      sexual_recruit,
                                      source_date,
                                      CSR_accession,
                                      mlg_id,
                                      mlg_geno_id,
                                      microsat_id
                                    FROM genotypes
                                    WHERE mlg_id IN ({snp*})",
                                   snp = input$selectGenotypes,
                                   .con = db)
    } else if(input$genotypeGroupChoice == 'csr') {
      genoOverview_sql <- glue_sql("Select 
                                      geno_name AS genotype,
                                      alt_names,
                                      origin_nursery,
                                      source_lat,
                                      source_long,
                                      sexual_recruit,
                                      source_date,
                                      CSR_accession,
                                      mlg_id,
                                      mlg_geno_id,
                                      microsat_id
                                    FROM genotypes
                                    WHERE CSR_accession IN ({csr*})",
                                   csr = input$selectGenotypes,
                                   .con = db)
    }
    
    genoOverview <- dbSendQuery(db, genoOverview_sql)
    t <- dbFetch(genoOverview) %>%
      mutate(sexual_recruit = ifelse(sexual_recruit == 0, 'No', 'Yes'),
             across(everything(), ~as.character(.)),
             source_location = paste(source_lat,source_long, sep = ", ")) %>%
      select(-c(source_lat, source_long)) %>%
      mutate(mlg_id = ifelse(startsWith(mlg_id, "HG") & !is.null(mlg_geno_id),       
                             paste0(mlg_id,"<br / >(<a href='https://coralsnp.science.psu.edu/reports/samples/with_genotype?sort_id=default&order=default&genotype_id=",mlg_geno_id,"&coral_mlg_clonal_id=", mlg_id, "' target='_blank'>View STAGdb record</a>)"),
                             NA)) %>%
      mutate(CSR_accession = ifelse(!is.null(CSR_accession) & !is.na(CSR_accession),
                                    paste0(CSR_accession,"<br />(<a href='http://crfcoralregistry.com/#main/0/registry/edit?id=",CSR_accession,"' target='_blank'>View CSR record</a>)"),
                                    NA)) %>%
      select(-mlg_geno_id) %>% #drop the mlg_geno_id col
      select(where(~!all(is.na(.)))) %>% # select only cols w/ data
      pivot_longer(cols = everything(),
                   names_to = 'Info',
                   values_to= 'Value') %>%
      group_by(Info) %>%
      summarise(Value = paste(unique(Value[!is.na(Value)]), collapse = "; "), .groups = "drop") %>%
      filter(Value != '') %>%
      mutate(Value = ifelse(Info %in% c('source_location', 'source_date'),
                            str_replace_all(Value, ';', ' <br /> '),
                            Value)) %>%
      mutate(Info = str_replace(Info, '_', ' ')) %>%
      drop_na() %>%
      mutate(Info = paste0(str_to_title(Info),":")) %>%
      mutate(Info = ifelse(Info == 'Mlg Id:', 'STAGdb MLG ID:',
                           ifelse(Info == 'Csr Accession:', 'CSR Accession:',
                                  ifelse(Info=='Microsat Id:','Microsatelite ID:',Info)))) %>%
      arrange(factor(Info, levels = c("Genotype:", "Alt Names:", "Origin Nursery:", "Source Location:",
                                      "Sexual Recruit:", "Source Date:", "Microsatelite ID:", "STAGdb MLG ID:",
                                      "CSR Accession:")))
  })
  
  output$genotypeOverivew <- renderDataTable({
    DT::datatable(genotypeOverviewdf(), rownames = F, selection = 'none',
                  options = list(dom = 't', ordering = F), escape = FALSE) %>%
      formatStyle(columns = c(2),
                  valueColumns = c(1),
                  cursor = styleEqual(c("Genotype:","STAGdb MLG ID:","CSR Accession:"),
                                      c('pointer', 'pointer', 'pointer'))) 
  })
  
  #observe table for clicking, change genotypeGroupingChoice & selectGenotypes
  observeEvent(input$genotypeOverivew_cell_clicked, {
    info = input$genotypeOverivew_cell_clicked
    
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    
    #collect desired genotype classification
    if(grepl('STAGdb', info$value, fixed = T)) {
      group_select <- 'snp'
    } else if(grepl('CSR', info$value, fixed = T)) {
      group_select <- 'csr'
    } else if(info$row == 1) {
      group_select <- 'name'
    } else{return()}
    
    #capture when multiple SNP matches -> give option to select specific CSR or local name
    if(nchar(str_extract(info$value,".*(?=<br )"))>40 & group_select=='csr'){
      showModal(modalDialog(
        title = "Multiple matching CSR Accession Records",
        h4("Select the CSR Accession Record you want to view:"),
        dataTableOutput("multipleCSRtable"),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if(group_select =='name' & dim(str_split(info$value, "; ", simplify = T))[2]>1) {
      showModal(modalDialog(
        title = "Multiple matching Local Nursery Names",
        h4("Select the Local Nursery Name you want to view:"),
        dataTableOutput("multipleNurseryNametable"),
        easyClose = TRUE,
        footer = NULL
      ))
    } else if(group_select == 'name') {
      changingGenotypeSelected(info$value)
      updateRadioButtons(session, "genotypeGroupChoice", selected = 'name')
    } else{
      changingGenotypeSelected(str_extract(info$value,'.*(?=<br )'))
      updateRadioButtons(session, "genotypeGroupChoice", selected = group_select)
    }
  })
  
  output$multipleCSRtable <- renderDataTable({
    t = geno_names %>%
      filter(mlg_id == input$selectGenotypes) %>%
      select(CSR_accession) %>%
      rename(`CSR Accession:` = CSR_accession)
    
    DT::datatable(t,
                  rownames=F, selection = 'none',
                  options = list(dom = 't', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click modal table with multiple csr and go to that specific csr
  observeEvent(input$multipleCSRtable_cell_clicked, {
    info = input$multipleCSRtable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 0) return()
    changingGenotypeSelected(info$value)
    updateRadioButtons(session, "genotypeGroupChoice", selected = 'csr') 
    removeModal()
  })
  
  output$multipleNurseryNametable <- renderDataTable({
    t = geno_names %>%
      filter(mlg_id == input$selectGenotypes) %>%
      select(geno_name) %>%
      rename(Genotype = geno_name)
    
    DT::datatable(t,
                  rownames =F, selection = 'none',
                  options = list(dom = 't', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click modal table with multiple nursery names and go to that specific nursery named genet
  observeEvent(input$multipleNurseryNametable_cell_clicked, {
    info = input$multipleNurseryNametable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 0) return()
    changingGenotypeSelected(info$value)
    updateRadioButtons(session, "genotypeGroupChoice", selected = 'name') 
    removeModal()
  })
  
  #observations and source of the selected genotype
  output$genoMap <- renderLeaflet({
    if(input$genotypeGroupChoice == 'name') {
      genoMap_sql <- glue_sql("Select
                                location_name as 'name',
                                location_lat AS 'lat',
                                location_long AS 'lng',
                                type,
                                trait_name AS 'trait',
                                value
                              FROM calculatedmeasurements AS m
                              INNER JOIN locations
                              ON locations.location_id = m.location_id
                              INNER JOIN genotypes
                              ON genotypes.geno_id = m.geno_id
                              INNER JOIN traits
                              ON traits.trait_id = m.trait_id
                              WHERE geno_name IN ({geno*})",
                              geno = input$selectGenotypes,
                              .con = db)
      genoMap <- dbSendQuery(db, genoMap_sql)
      genoMapPoints <- dbFetch(genoMap) %>%
        distinct() %>%
        group_by(name, lat, lng, type) %>%
        summarise(n = length(value),
                  traits = paste(unique(trait), collapse= ", "),
                  .groups = "drop")
      
      genoSource_sql <- glue_sql("Select
                                  geno_name AS 'name',
                                  source_lat AS 'lat',
                                  source_long AS 'lng',
                                  source_date AS 'date'
                              FROM genotypes
                              WHERE geno_name IN ({geno*})",
                                 geno = input$selectGenotypes,
                                 .con = db)
      geno_sourceS <- dbSendQuery(db, genoSource_sql)
      geno_source <- dbFetch(geno_sourceS)
    } else if(input$genotypeGroupChoice == 'snp') {
      genoMap_sql <- glue_sql("Select
                                location_name as 'name',
                                location_lat AS 'lat',
                                location_long AS 'lng',
                                type,
                                trait_name AS 'trait',
                                value
                              FROM calculatedmeasurements AS m
                              INNER JOIN locations
                              ON locations.location_id = m.location_id
                              INNER JOIN genotypes
                              ON genotypes.geno_id = m.geno_id
                              INNER JOIN traits
                              ON traits.trait_id = m.trait_id
                              WHERE mlg_id IN ({snp*})",
                              snp = input$selectGenotypes,
                              .con = db)
      genoMap <- dbSendQuery(db, genoMap_sql)
      genoMapPoints <- dbFetch(genoMap) %>%
        distinct() %>%
        group_by(name, lat, lng, type) %>%
        summarise(n = length(value),
                  traits = paste(unique(trait), collapse= ", "),
                  .groups = "drop")
      
      genoSource_sql <- glue_sql("Select
                                  geno_name AS 'name',
                                  source_lat AS 'lat',
                                  source_long AS 'lng',
                                  source_date AS 'date'
                              FROM genotypes
                              WHERE mlg_id IN ({snp*})",
                                 snp = input$selectGenotypes,
                                 .con = db)
      geno_sourceS <- dbSendQuery(db, genoSource_sql)
      geno_source <- dbFetch(geno_sourceS)
    } else if(input$genotypeGroupChoice == 'csr') {
      genoMap_sql <- glue_sql("Select
                                location_name as 'name',
                                location_lat AS 'lat',
                                location_long AS 'lng',
                                type,
                                trait_name AS 'trait',
                                value
                              FROM calculatedmeasurements AS m
                              INNER JOIN locations
                              ON locations.location_id = m.location_id
                              INNER JOIN genotypes
                              ON genotypes.geno_id = m.geno_id
                              INNER JOIN traits
                              ON traits.trait_id = m.trait_id
                              WHERE CSR_accession IN ({csr*})",
                              csr = input$selectGenotypes,
                              .con = db)
      genoMap <- dbSendQuery(db, genoMap_sql)
      genoMapPoints <- dbFetch(genoMap) %>%
        distinct() %>%
        group_by(name, lat, lng, type) %>%
        summarise(n = length(value),
                  traits = paste(unique(trait), collapse= ", "),
                  .groups = "drop")
      
      genoSource_sql <- glue_sql("Select
                                  geno_name AS 'name',
                                  source_lat AS 'lat',
                                  source_long AS 'lng',
                                  source_date AS 'date'
                              FROM genotypes
                              WHERE CSR_accession IN ({csr*})",
                                 csr = input$selectGenotypes,
                                 .con = db)
      geno_sourceS <- dbSendQuery(db, genoSource_sql)
      geno_source <- dbFetch(geno_sourceS)
    }
    
    #set factors
    geno_source$type <- as.factor("source")
    #change nursery type to nursery
    genoMapPoints$type[which(genoMapPoints$type %in% c('ln','bn'))] <- 'nursery'
    genoMapPoints$type <- as.factor(genoMapPoints$type)
    
    factpal <- colorFactor(c("darkviolet", "red", "green", "blue"),
                           c("lab", "nursery", "reef", "source"))
    if(nrow(genoMapPoints)>0 | nrow(geno_source) > 0){
      leaflet() %>%
        addTiles() %>%
        addCircleMarkers(data = genoMapPoints, group = "Observations", lng = ~lng, lat = ~lat,
                         popup = ~paste0("<strong>",name,"</strong><br>",
                                         "# of fragments: ", n, "<br>",
                                         "Traits: ", traits),
                         label =  ~as.character(name),
                         stroke = FALSE, radius = 6, fillOpacity = 0.7,
                         color = ~factpal(type)) %>%
        addCircleMarkers(data = geno_source, group = "Genet Source Locations",
                         lng = ~lng, lat = ~lat, label = ~as.character(name),
                         popup = ~paste0("<strong>",name,":</strong><br>",
                                         "Collection Date: ", date),
                         stroke = FALSE, radius = 6, color = ~factpal(type), fillOpacity = 0.7) %>%
        addLegend(colors = c("darkviolet", "red", "green", "blue"),
                  labels = c("lab", "nursery", "reef", "source location"), opacity = 0.5,
                  position = "bottomright") %>%
        addScaleBar(position = 'bottomleft',
                    options = scaleBarOptions(imperial = F))}
  })
  
  genotypesDF <- reactive({
    data <- genoDefinition(input$genotypeGroupChoice) %>%
            filter(genotype == input$selectGenotypes)
    
  })
  
  genotypeStatsDF <-reactive({
    genotypesDF() %>%
      group_by(trait, method, unit, rank) %>% 
      summarise(mean = mean(value, na.rm = T),
                sd = sd(value, na.rm = T),
                n = n(),
                min = min(value),
                max = max(value), 
                datasets = length(unique(datafile_name)),
                .groups = "drop") %>%
      mutate(across(c(mean,sd,min,max), ~format(round(.x, 3), nsmall=3)),
             across(everything(), ~as.character(.))) %>%
      distinct() %>%
      arrange(trait, method, unit) %>%
      mutate(rank = ifelse(rank==1,
                           '<img src="1greenLight_small.png" height="25" />',
                           ifelse(rank==2,
                                  '<img src="2yellowLight_small.png" height="25" />',
                                  '<img src="3redLight_small.png" height="25" />'))) %>%
      relocate(rank, .after='mean')
  })
  
  output$genotypeTable <- renderDataTable({
    globalRank_text <- tags$span(
      "Global Rank", 
      infoBtn('traitGlobalRank') %>% 
        bsPopover(title = "Global Rank",
                  content = "This column is a quick reference tool to indicate above average, average, or below average performance for a given genet for the applicable trait, method, and unit based off the total data and not the applied data filters. Green indicates the genet falls in the first tercile, yellow indicates the genet falls in the middle tercile, and red indicates the genet falls in the bottom tercile.",
                  placement = "top",
                  trigger = "hover")
    ) %>% 
      as.character()
    
    mean_text <- tags$span(
      "Mean", 
      infoBtn('traitMean') %>% 
        spsComps::bsPopover(title = "Mean",
                            content = "This mean is calculated from all data. Click the mean to view a modified box-plot to understand how the calculated mean compares to the population.",
                            placement = "top",
                            trigger = "hover")
    ) %>% 
      as.character()
    
    DT::datatable(genotypeStatsDF() %>%
                    rename_with(str_to_title) %>%
                    rename(!!mean_text:=Mean,
                           !!globalRank_text:=Rank),
                  rownames = F, selection = 'none', escape = F,
                  options = list(dom = 'tip', pageLength = 1000)) %>%
      formatStyle(c(1,4,10), cursor = "pointer")
  })
  
  #Change to selected traits tab
  observeEvent(input$genotypeTable_cell_clicked, {
    info = input$genotypeTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 0) return()
    updateTabItems(session, 'tabs', selected = 'traits_tab')
    
    changingTraitSelected(info$value)
    changingTraitSelectedMethod(genotypeStatsDF()$method[info$row])
    
    updateRadioButtons(session, "TraitAnalysisGenotypeGroupChoice", selected = input$genotypeGroupChoice)
    
    newSelectedTrait <- case_when(info$value %in% biomechanicalTraits ~  'biomechanical',
                                  info$value %in% bleachingTraits ~ 'bleaching',
                                  info$value %in% diseaseTraits ~ 'disease',
                                  info$value %in% growthTraits ~ 'growth rates',
                                  info$value %in% hostPhysiologyTraits  ~ 'host physiology',
                                  info$value %in% reproductionTraits ~ 'reproduction',
                                  info$value %in% woundHealingTraits ~'wound healing',
                                  TRUE ~ 'growth rates')
    
    updateSelectInput(session, 'traitSelector',
                      selected = newSelectedTrait)
  })
  
  #Popup list of datasets when num datasets is clicked in genotype report card
  observeEvent(input$genotypeTable_cell_clicked, {
    info = input$genotypeTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 9) return()
    showModal(modalDialog(
      title = paste("Datasets that contain the", input$selectGenotypes,
                    "genet and the", genotypeStatsDF()$trait[info$row],
                    "trait."),
      dataTableOutput("genotypesModalTable"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  output$genotypesModalTable <- renderDataTable({
    list = genotypesDF() %>%
      filter(trait == genotypeStatsDF()$trait[input$genotypeTable_cell_clicked$row]) %>%
      select(datafile_name) %>% distinct() %>%
      left_join(datasets, by = "datafile_name") %>%
      select(`Data Source Title` = title)
    
    DT::datatable(list %>% rename_with(str_to_title),
                  rownames = T, selection = 'none',
                  options = list(dom = 'tf', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click datasource brings you to that datasource page
  observeEvent(input$genotypesModalTable_cell_clicked, {
    info = input$genotypesModalTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    updateTabItems(session, 'tabs', selected = 'dataSources')
    updateSelectInput(session, 'selectDataSource', selected = info$value)
    removeModal()
  })
  
  #Modal for plot
  observeEvent(input$genotypeTable_cell_clicked, {
    info = input$genotypeTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 3rd column
    if (is.null(info$value) || info$col != 3) return()
    showModal(modalDialog(
      id = 'plotModal',
      title = paste("Comparison of", input$selectGenotypes,
                    "to the population mean of the", genotypeStatsDF()$trait[info$row],
                    "trait."),
      span("Loading...", id="UpdateAnimate",
           style="font-size:1.5em; margin-left: auto;
                     margin-right: auto;"),
      imageOutput("genotypesPlot", height="100%", width = "100%"),
      easyClose = TRUE,
      footer = NULL,
    ))
  })
  
  output$genotypesPlot <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    info <- input$genotypeTable_cell_clicked
    selectedTrait <- genotypeStatsDF()$trait[info$row]
    selectedMethod <- genotypeStatsDF()$method[info$row]
    selectedUnit <- genotypeStatsDF()$unit[info$row]
    
    temp <- genoDefinition(input$genotypeGroupChoice) %>%
      filter(trait == selectedTrait &
               method == selectedMethod &
               unit == selectedUnit) %>%
    #remove extreme outliers for z-score calculation
    filter(abs(scale(value,
                     center = median(value, na.rm = T),
                     scale = mad(value, na.rm=T)))<=5) %>%
    #calculate genotype z-scores
    group_by(genotype) %>%
      summarise(mean = mean(value,na.rm=T)) %>%
      mutate(z = scale(mean)) %>%
      select(-mean) %>%
      #add data back in
      right_join(genoDefinition(input$genotypeGroupChoice) %>%
                   filter(trait == selectedTrait &
                            method == selectedMethod &
                            unit == selectedUnit),
                 by="genotype") %>%
      #remove extreme outliers for data display
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T)))<=5)
    
    temp2 <- genotypesDF() %>%
      filter(trait == selectedTrait & method == selectedMethod &
               unit == selectedUnit)
    z <- temp %>% filter(genotype == unique(temp2$genotype)) %>%
      pull(z) %>% unique()
    
    
    png(outfile, 
        width = 500*4, 
        height = 250*4,
        res = 72*4)
    print(modalPlot(temp,temp2,z))
    dev.off()
    shinyjs::hideElement(id='UpdateAnimate')
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 500,
         height = 250,
         alt = "This is alternate text")
  }, deleteFile = T)
  
  #create .csv file with contact + citation information
  genotypeCitations <- reactive({ 
    df <- genotypesDF()
    
    datasets %>%
      filter(datafile_name %in%
               (df %>% pull(datafile_name) %>% unique())) %>%
      select(title, author, email, DOI) %>%
      mutate(DOI = str_extract(DOI, "(?<=>)(.*?)(?=<)"))
  })
  
  output$downloadGenotype <- downloadHandler(
    filename = paste0("AcDC_",
                      gsub("(^|[^[:alnum:]])([[:alnum:]])",
                           "\\U\\2",
                           input$selectGenotypes,
                           perl = TRUE),
                      "_genotypeReportCard.zip"),
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      #create csv that will be in the zip folder
      fileNames <- c("dataExport.csv","citationExport.csv")
      write.csv(genotypeStatsDF() %>%
                  select(-rank) %>%
                  relocate(genotype, .before="trait"),
                #this is not changing for some reason
                #mutate(unit = gsub("<\\U\\+0394>", "change", unit)),
                "dataExport.csv", row.names = F)
      
      write.csv(genotypeCitations(), "citationExport.csv", row.names = F)
      
      #create the zip file
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      showModal(modalDialog(
        title = 'Warning: You must Assume the Responsibility for Downloaded Data',
        p("The data you just downloaded are the summary statistics of data provided by the authors. The National Oceanographic and Atmospheric Administration makes no claims as to the quality assurance of the data presented here. All data are provided as a scientific product on an \'as is\' basis and the user assumes responsibility for its use. Please cite the authors of the original data and the AcDC website in all publications and presentations that use the downloaded data."),
        footer =  modalButton("Confirm")
      ))
    }
  )
  
  #Genet Comparison ---------------------------------------------------
  
  genotypeReportSelectorChoices <- reactive({
    #local name selected
    if(input$genotypeTypeReportChoice == 'name') {
      geno_names$geno_name
      #snp selected
    } else if(input$genotypeTypeReportChoice == 'snp') {
      geno_names$mlg_id[!is.na(geno_names$mlg_id)]
    } else if(input$genotypeTypeReportChoice == 'csr') {
      geno_names$CSR_accession[!is.na(geno_names$CSR_accession)]
    } 
  })
  
  output$genotypesReportChoice <-renderUI({
    if(input$genotypeTypeReportChoice == 'map') {
      tagList(
        p("Use the polygon tool to select genotypes by their source locations."),
        leafletOutput("genotypesReportMapSelector")
      )
    } else {
      selectizeInput('genotypesReportSelector', label = 'Select the genotypes you want to analyze:',
                     choices = genotypeReportSelectorChoices(),
                     multiple = T)
    }
  })
  
  #observations and source of the selected genotype
  output$genotypesReportMapSelector <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = geno_source, group = "Genet Source Locations", lng = ~source_long, lat = ~source_lat,
                       popup = ~paste0("<strong>",geno_name,"</strong>"),
                       label =  ~as.character(geno_name),
                       stroke = FALSE, radius = 6, fillOpacity = 0.7,
                       color = 'red') %>%
      addLegend(colors = c("red"),
                labels = c("source location"), opacity = 0.5,
                position = "bottomright") %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = drawPolygonOptions(shapeOptions=drawShapeOptions(fillOpacity = .4,
                                                                                       color = 'black',
                                                                                       fillColor = '#fff',
                                                                                       weight = 3)),
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE) %>%
      addScaleBar(position = 'bottomleft',
                  options = scaleBarOptions(imperial = F))
  })
  
  genotypeReportDrawnValues <- reactiveValues(activePolygon=FALSE)
  
  #prevent adding polygons after the first is created
  observeEvent(input$genotypesReportMapSelector_draw_new_feature,{
    genotypeReportDrawnValues$activePolygon <- TRUE
    leafletProxy('genotypesReportMapSelector') %>%
      removeDrawToolbar() %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = FALSE,
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE)
  })
  
  #add original toolbar back in when polygon is deleted
  observeEvent(input$genotypesReportMapSelector_draw_deleted_features,{
    genotypeReportDrawnValues$activePolygon <- FALSE
    leafletProxy('genotypesReportMapSelector') %>%
      removeDrawToolbar() %>%
      addDrawToolbar(rectangleOptions = FALSE,
                     editOptions = editToolbarOptions(edit = FALSE,
                                                      remove = TRUE),
                     polylineOptions = FALSE,
                     polygonOptions = drawPolygonOptions(shapeOptions=drawShapeOptions(fillOpacity = .4,
                                                                                       color = 'black',
                                                                                       fillColor = '#fff',
                                                                                       weight = 3)),
                     circleOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     markerOptions = FALSE)
  })
  
  #if raidobutton is clicked after map has already been open
  observeEvent(input$genotypeTypeReportChoice, {
    if(input$genotypeTypeReportChoice != 'map') {
      genotypeReportDrawnValues$activePolygon <- FALSE
      removeUI('genotypesReportMapSelector')
    }
  })
  
  #filter the geno_source DF by the polygon
  PolyFilteredGenotypeSelection <- eventReactive(input$genotypesReportMapSelector_draw_new_feature,{
    poly <- input$genotypesReportMapSelector_draw_new_feature$geometry$coordinates[[1]]
    bounds <- data.frame(long = map_dbl(poly, `[[`,1),
                         lat = map_dbl(poly, `[[`, 2))
    geno_source %>%
      mutate(inside = point.in.polygon(source_long, source_lat, bounds$long, bounds$lat)) %>%
      filter(inside >0) %>%
      select(-inside)
  })
  
  #filtering of gentypeReportsdf
  genotypeReportDF <- reactive({
    if(input$genotypeTypeReportChoice == 'map' & !genotypeReportDrawnValues$activePolygon) {
      #return blank data to the table
      temp <- corals[0,]
    } else if(input$genotypeTypeReportChoice == 'map' & genotypeReportDrawnValues$activePolygon) {
      temp <- genoDefinition('name') %>%
        filter(genotype %in% PolyFilteredGenotypeSelection()$geno_name)
    } else {
      temp <- genoDefinition(input$genotypeTypeReportChoice) %>%
        filter(genotype %in% input$genotypesReportSelector)
    }
    temp
  })
  
  #output for genotypeReport
  genotypeReportStatsDF <-reactive({
    genotypeReportDF() %>%
      group_by(trait, method, unit, genotype, rank) %>% 
      summarise(mean = mean(value, na.rm = T),
                sd = sd(value, na.rm = T),
                n = n(),
                min = min(value),
                max = max(value), 
                datasets = length(unique(datafile_name)),
                .groups = "drop") %>%
      mutate(across(c(mean,sd,min,max), ~format(round(.x, 3), nsmall=3)),
             across(everything(), ~as.character(.))) %>%
      distinct() %>%
      arrange(genotype) %>%
      mutate(rank = ifelse(rank==1,
                           '<img src="1greenLight_small.png" height="25" />',
                           ifelse(rank==2,
                                  '<img src="2yellowLight_small.png" height="25" />',
                                  '<img src="3redLight_small.png" height="25" />'))) %>%
      relocate(rank, .after='mean') %>%
      relocate(genotype, .before='trait')
  })
  
  output$genotypeReport <- renderDataTable({
    globalRank_text <- tags$span(
      "Global Rank", 
      infoBtn('traitGlobalRank') %>% 
        bsPopover(title = "Global Rank",
                  content = "This column is a quick reference tool to indicate above average, average, or below average performance for a given genet for the applicable trait, method, and unit based off the total data and not the applied data filters. Green indicates the genet falls in the first tercile, yellow indicates the genet falls in the middle tercile, and red indicates the genet falls in the bottom tercile.",
                  placement = "top",
                  trigger = "hover")
    ) %>% 
      as.character()
    
    mean_text <- tags$span(
      "Mean", 
      infoBtn('traitMean') %>% 
        spsComps::bsPopover(title = "Mean",
                            content = "This mean is calculated from all data. Click the mean to view a modified box-plot to understand how the calculated mean compares to the population.",
                            placement = "top",
                            trigger = "hover")
    ) %>% 
      as.character()
    
    
    DT::datatable(genotypeReportStatsDF() %>%
                    rename_with(str_to_title) %>%
                    rename(!!mean_text:=Mean,
                           !!globalRank_text:=Rank),
                  rownames = F, selection = 'none', escape = F, 
                  options = list(dom = 'tip', pageLength = 1000)) %>%
      formatStyle(c(1,2,5,11), cursor = "pointer")
  })
  
  #Change to selected Genotypes tab
  observeEvent(input$genotypeReport_cell_clicked, {
    info = input$genotypeReport_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 0) return()
    
    changingGenotypeSelected(info$value)
    
    #make genotype report radio selection the same as the genet comparison
    if(input$genotypeTypeReportChoice == 'name') {
      updateRadioButtons(session, "genotypeGroupChoice", selected = 'name')
    } else if(input$genotypeTypeReportChoice == 'snp') {
      updateRadioButtons(session, "genotypeGroupChoice", selected = 'snp')
    } else if(input$genotypeTypeReportChoice == 'csr') {
      updateRadioButtons(session, "genotypeGroupChoice", selected = 'csr')
    } else if(input$genotypeTypeReportChoice == 'map') {
      updateRadioButtons(session, "genotypeGroupChoice", selected = 'name')
    }
    
    updateSelectInput(session, 'selectGenotypes', selected = info$value)
    
    updateTabItems(session, 'tabs', selected = 'genotypes_tab')
  })
  
  #Change to selected Traits tab
  observeEvent(input$genotypeReport_cell_clicked, {
    info = input$genotypeReport_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    updateTabItems(session, 'tabs', selected = 'traits_tab')
    
    changingTraitSelected(info$value)
    changingTraitSelectedMethod(genotypeReportStatsDF()$method[info$row])
    
    #make trait analysis genotype radio selection the same as the genet comparison
    if(input$genotypeTypeReportChoice == 'name') {
      updateRadioButtons(session, "TraitAnalysisGenotypeGroupChoice", selected = 'name')
    } else if(input$genotypeTypeReportChoice == 'snp') {
      updateRadioButtons(session, "TraitAnalysisGenotypeGroupChoice", selected = 'snp')
    } else if(input$genotypeTypeReportChoice == 'csr') {
      updateRadioButtons(session, "TraitAnalysisGenotypeGroupChoice", selected = 'csr')
    } else if(input$genotypeTypeReportChoice == 'map') {
      updateRadioButtons(session, "TraitAnalysisGenotypeGroupChoice", selected = 'name')
    }
    
    newSelectedTrait <- case_when(info$value %in% biomechanicalTraits ~  'biomechanical',
                                  info$value %in% bleachingTraits ~ 'bleaching',
                                  info$value %in% diseaseTraits ~ 'disease',
                                  info$value %in% growthTraits ~ 'growth rates',
                                  info$value %in% hostPhysiologyTraits  ~ 'host physiology',
                                  info$value %in% reproductionTraits ~ 'reproduction',
                                  info$value %in% woundHealingTraits ~'wound healing',
                                  TRUE ~ 'growth rates')
    
    updateSelectInput(session, 'traitSelector',
                      selected = newSelectedTrait)
  })
  
  #Popup list of datasets when num datasets is clicked
  observeEvent(input$genotypeReport_cell_clicked, {
    info = input$genotypeReport_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 11th column
    if (is.null(info$value) || info$col != 10) return()
    showModal(modalDialog(
      title = paste("Datasets that contain the", genotypeReportStatsDF()$genotype[info$row],
                    "genet and the", genotypeReportStatsDF()$trait[info$row],
                    "trait."),
      dataTableOutput("genotypeReportModalTable"),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  output$genotypeReportModalTable <- renderDataTable({
    list = genotypeReportDF() %>%
      filter(trait == genotypeReportStatsDF()$trait[input$genotypeReport_cell_clicked$row] &
               genotype == genotypeReportStatsDF()$genotype[input$genotypeReport_cell_clicked$row]) %>%
      select(datafile_name) %>% distinct() %>%
      left_join(datasets, by = "datafile_name") %>%
      select(`Data Source Title` = title)
    
    DT::datatable(list %>% rename_with(str_to_title),
                  rownames = T, selection = 'none',
                  options = list(dom = 'tf', pageLength = 10)) %>%
      formatStyle(1, cursor = 'pointer')
  })
  
  #click datasource brings you to that datasource page
  observeEvent(input$genotypeReportModalTable_cell_clicked, {
    info = input$genotypeReportModalTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st column
    if (is.null(info$value) || info$col != 1) return()
    updateTabItems(session, 'tabs', selected = 'dataSources')
    updateSelectInput(session, 'selectDataSource', selected = info$value)
    removeModal()
  })
  
  #Modal for plot
  observeEvent(input$genotypeReport_cell_clicked, {
    info = input$genotypeReport_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 5th column
    if (is.null(info$value) || info$col != 4) return()
    showModal(modalDialog(
      id = 'plotModal',
      title = paste("Comparison of", genotypeReportStatsDF()$genotype[info$row],
                    "to the population mean of the", genotypeReportStatsDF()$trait[info$row],
                    "trait."),
      span("Loading...", id="UpdateAnimate",
           style="font-size:1.5em; margin-left: auto;
                     margin-right: auto;"),
      imageOutput("genotypesReportPlot"),
      easyClose = TRUE,
      footer = NULL,
    ))
  })
  
  output$genotypesReportPlot <- renderImage({
    outfile <- tempfile(fileext = '.png')
    
    info <- input$genotypeReport_cell_clicked
    selectedTrait <- genotypeReportStatsDF()$trait[info$row]
    selectedUnit <- genotypeReportStatsDF()$unit[info$row]
    selectedMethod <- genotypeReportStatsDF()$method[info$row]
    selectedGenotype <- genotypeReportStatsDF()$genotype[info$row]
    
    temp <- genoDefinition(input$genotypeTypeReportChoice) %>%
      filter(trait == selectedTrait &
               unit == selectedUnit &
               method ==  selectedMethod) %>%
      #remove extreme outliers for z-score calculation
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T)))<=5) %>%
    #calculate genotype z-scores
    group_by(genotype) %>%
      summarise(mean = mean(value,na.rm=T)) %>%
      mutate(z = scale(mean)) %>%
      select(-mean) %>%
      #add data back in
      right_join(genoDefinition(input$genotypeTypeReportChoice) %>%
                   filter(trait == selectedTrait &
                            unit == selectedUnit &
                            method ==  selectedMethod),
                 by="genotype") %>%
      #remove extreme outliers for data display
      filter(abs(scale(value,
                       center = median(value, na.rm = T),
                       scale = mad(value, na.rm=T)))<=5)
    
    temp2 <- genotypeReportDF() %>%
      filter(trait == selectedTrait & 
               unit == selectedUnit &
               method == selectedMethod &
               genotype == selectedGenotype)
    
    z <- temp %>% filter(genotype == selectedGenotype) %>%
      pull(z) %>% unique()
    
    png(outfile, 
        width = 500*4, 
        height = 250*4,
        res = 72*4)
    print(modalPlot(temp,temp2,z))
    dev.off()
    shinyjs::hideElement(id='UpdateAnimate')
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 500,
         height = 250,
         alt = "This is alternate text")
  }, deleteFile = T)
  
  #create .csv file with contact + citation information
  genetCompareCitations <- reactive({ 
    df <- genotypeReportDF()
    
    datasets %>%
      filter(datafile_name %in%
               (df %>% pull(datafile_name) %>% unique())) %>%
      select(title, author, email, DOI) %>%
      mutate(DOI = str_extract(DOI, "(?<=>)(.*?)(?=<)"))
  })
  
  output$downloadGenotypeReport <- downloadHandler(
    filename = "AcDC_FullReportCards.zip",
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      #create csv that will be in the zip folder
      fileNames <- c("dataExport.csv","citationExport.csv")
      write.csv(genotypeReportStatsDF() %>%
                  #this mutate is not working for some reason??
                  #mutate(unit = gsub("<\\U\\+0394>", "change", unit)) %>%
                  select(-rank),
                "dataExport.csv", row.names = F)
      
      write.csv(genetCompareCitations(), "citationExport.csv", row.names = F)
      
      #create the zip file
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      showModal(modalDialog(
        title = 'Warning: You Must Assume the Responsibility for Downloaded Data',
        p("The data you just downloaded are the summary statistics of data provided by the authors. The National Oceanographic and Atmospheric Administration makes no claims as to the quality assurance of the data presented here. All data are provided as a scientific product on an \'as is\' basis and the user assumes responsibility for its use. Please cite the authors of the original data and the AcDC website in all publications and presentations that use the downloaded data."),
        footer =  modalButton("Confirm")
      ))
    }
  )
  
  #DataSources Tab----------------------------------------------------------------------
  #selectDataSource Filter
  observe({
    if('biomechanical' %in% input$traits_dataSources) {traits = biomechanicalTraits}
    if('bleaching' %in% input$traits_dataSources) {traits = bleachingTraits}
    if('disease' %in% input$traits_dataSources) {traits = diseaseTraits}
    if('growth rates' %in% input$traits_dataSources) {traits = growthTraits}
    if('host physiology' %in% input$traits_dataSources) {traits = hostPhysiologyTraits}
    if('reproduction' %in% input$traits_dataSources) {traits = reproductionTraits}
    if('wound healing' %in% input$traits_dataSources) {traits = woundHealingTraits}

    if(input$traits_dataSources != "" & input$genotypes_dataSources != "") {
      newChoices <- corals %>%
        filter(trait %in% traits & genotype %in% input$genotypes_dataSources) %>%
        select(datafile_name) %>% distinct() %>% 
        left_join(datasets, by = "datafile_name") %>%
        filter(!is.na(title)) %>%
        pull(title)
      
      updateSelectInput(session, 'selectDataSource', choices=newChoices)}
    
    if(input$traits_dataSources != "" & input$genotypes_dataSources == "") {
      newChoices <- corals %>%
        filter(trait %in% traits) %>%
        select(datafile_name) %>% distinct() %>%
        left_join(datasets, by = "datafile_name") %>%
        filter(!is.na(title)) %>%
        pull(title)
      
      updateSelectInput(session, 'selectDataSource', choices=newChoices)}
    
    if(input$traits_dataSources == "" & input$genotypes_dataSources != "") {
      newChoices <- corals %>%
        filter(genotype %in% input$genotypes_dataSources) %>%
        select(datafile_name) %>% distinct() %>%
        left_join(datasets, by = "datafile_name") %>%
        filter(!is.na(title)) %>%
        pull(title)
      
      updateSelectInput(session, 'selectDataSource', choices=newChoices)}
    if(input$traits_dataSources == "" & input$genotypes_dataSources == "") {
      newChoices <- corals %>%
        select(datafile_name) %>% distinct() %>%
        left_join(datasets, by = "datafile_name") %>%
        filter(!is.na(title)) %>%
        pull(title)
      
      updateSelectInput(session, 'selectDataSource', choices=newChoices)}
  })
  
  output$dataSource_header <- renderUI({
    HTML(paste('<h3>Dataset:', input$selectDataSource, '</h3>'))
  })
  
  dataSourceDF <- reactive({
    datasets %>%
      filter(title == input$selectDataSource) %>%
      mutate(across(everything(),~as.character(.)),
             remaining_tissue = ifelse(remaining_tissue=='TRUE',
                                       'Please contact the author to discuss if samples are still preserved and available for additional analysis.',
                                       NA)) %>%
      rename(`Methods & Procedures` = meth_proc) %>%
      pivot_longer(cols = everything(),
                   names_to = 'Info',
                   values_to= 'Value') %>%
      filter(!Info %in% c('title','dataset_id', 'datafile_name', 'complete')) %>%
      mutate(Info = paste0(str_to_title(str_replace_all(Info,"_"," ")),":")) %>%
      mutate(Info = ifelse(Info %in% c('Orcid:', 'Doi:'), toupper(Info), Info)) %>%
      drop_na()
  })
  
  output$dataSourceOverview <- renderTable({
    dataSourceDF()
  }, sanitize.text.function = function(x) x)
  
  dataSourceRaw <- reactive({
    selectedTitle <- input$selectDataSource
    selectedDatafile <- datasets %>%
      filter(title == selectedTitle) %>% pull(datafile_name)
    #set up paramatized query for measurements table
    data_sql <- glue_sql("Select 
                                geno_name AS genotype,
                                m.observation_id AS observation_id,
                                location_name AS location,
                                location_lat AS latitude,
                                location_long AS longitude,
                                depth,
                                type,
                                trait_name AS trait,
                              	value,
                              	unit,
                              	method,
                              	datafile_name
                              FROM
	                            measurements AS m
                            	INNER JOIN
                            	observations AS o
                            	ON m.observation_id = o.observation_id
                                INNER JOIN
                                standards AS s
                                ON m.standard_id = s.standard_id
                                INNER JOIN
                                genotypes AS g
                                ON m.geno_id = g.geno_id
                                INNER JOIN
                                traits AS t
                                ON m.trait_id = t.trait_id
                                INNER JOIN
                                locations AS l
                                ON o.location_id = l.location_id
                                INNER JOIN
                                methods as meth
                                ON m.method_id = meth.method_id
                                INNER JOIN
                                datasets AS d
                                ON d.dataset_id = o.dataset_id
                                WHERE datafile_name IN ({dataset*})
                                AND d.private != 1
                                ORDER BY m.observation_id",
                         dataset = selectedDatafile,
                         .con = db)
    data <- dbSendQuery(db, data_sql)
    dbFetch(data) %>% distinct()
  })
  
  #create .csv file with contact + citation information
  dataSourceCitations <- reactive({ 
    datasets %>%
      filter(title == input$selectDataSource) %>%
      select(title, author, email, DOI) %>%
      mutate(DOI = str_extract(DOI, "(?<=>)(.*?)(?=<)"))
  })
  
  is_private <- reactive({
    datasets %>%
      filter(title == input$selectDataSource) %>%
      pull(private)
  })
  
  #display download button if not private
  output$downloadDataSource_output <- renderUI({
    if(is_private() == 'FALSE') {
      downloadButton("downloadDataSource", "Download Data")
    }
  })
  
  output$downloadDataSource <- downloadHandler(
    filename = "ACDC_dataExport.zip",
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      
      #create csv that will be in the zip folder
      fileNames <- c("dataExport.csv","citationExport.csv")
      write.csv(dataSourceRaw(), "dataExport.csv", row.names = F)
      write.csv(dataSourceCitations(), "citationExport.csv", row.names = F)
      
      #create the zip file
      #zip(file,fileNames)
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      showModal(modalDialog(
        title = 'Warning: You Must Assume the Responsibility for Downloaded Data',
        p("The data you just downloaded is the sole ownership of the authors who provided the data to us. The data are freely usable, but must be properly cited in any reports, publications, presentations or other relevant products. A csv file containing the author's contact information has also been downloaded for you to assist in proper citation. Further, the National Oceanographic and Atmospheric Administration makes no claims as to the quality assurance of the data presented here. All data are provided as a scientific product on an \'as is\' basis and the user assumes responsibility for its use."),
        footer =  modalButton("Confirm")
      ))
    }
  )
  
  # Data Submission Tab -----
  
  
  #display download button if not private
  output$acdcSubmissionDownload_output <- renderUI({
    downloadButton("acdcSubmissionDownload", "Download Submission Forms")
  })
  
  output$acdcSubmissionDownload <- downloadHandler(
    filename = "acdcSubmission.zip",
    content = function(file){
      #go to a temp dir to avoid permission issues
      #owd <- setwd(tempdir())
      #on.exit(setwd(owd))
      
      #create file names that will be in the zip folder
      fileNames <- c("dataSubmissionGuide.pdf",'dataSubmission.csv','metadataTemplate.rtf')
      wd <- getwd()
      FullPath=path.expand(paste0(wd,"/www/",fileNames))
      
      tmpdir <- tempdir()
      setwd(tmpdir)
      #create excel sheet
      databaseNamingConventionsList <- list(traits = dbGetQuery(db, "SELECT trait_name AS 'trait' FROM traits"),
                                            genotypes = dbGetQuery(db, "SELECT geno_name, source_lat, source_long, CSR_accession AS 'CSR ID', mlg_id AS 'STAGdb MLG ID' FROM genotypes"),
                                            methods = dbGetQuery(db, "SELECT method, notes FROM methods"),
                                            units = dbGetQuery(db, "SELECT name, unit FROM standards"))
      write.xlsx(databaseNamingConventionsList, file="databaseNamingConventions.xlsx")
      fileNames <- c(fileNames, "databaseNamingConventions.xlsx")
      
      #copy files from www folder into zip output
      lapply(FullPath,function(i) file.copy(i,tmpdir))
      
      #create the zip file
      system2("zip", args=(paste(file,fileNames,sep=" ")))
      
      setwd(wd)
      
      showModal(modalDialog(
        title = 'Data Submission Folder Downloaded',
        "A zip file has been downloaded containing all the necessary forms and guides to submit your data to AcDC. Please follow the instructions carefully, and please contact us if you have any questions.",
        easyClose = TRUE,
        footer = NULL
      ))
    }, contentType = "application/zip"
  )
}

# Run the application---- 
shinyApp(ui, server,  onStop(function() {dbDisconnect(db)}))

#' # Regional Population Projections for Israel
#' 
#' ## Copyright
#'  
#' Copyright 2017 National Economic Council, Prime Minister Office of Israel
#' 
#' This program is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#'  
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'  
#' You should have received a copy of the GNU General Public License
#' along with this program.  If not, see <http://www.gnu.org/licenses/>.
#'  
#' ## Dependencies
#' This program was coded in R 3.4, using RStudio 1.0.143.
#' Additional packages: tidyverse, readxl 

library(tidyverse)
library(readxl)
rm(list=ls())

#' ## Load Model Data

load("data/model/vectors.Rdata")
load("data/model/base.Rdata")
load("data/model/fertility.Rdata")
load("data/model/mortality.Rdata")
load("data/model/intmig.Rdata")
load("data/model/extmig.Rdata")
load("data/model/notmarried.Rdata")
load("data/model/housing_plans.Rdata")
load("data/model/CBS2065.Rdata")

#' Load Model Functions

source("scripts/model_fn.R")

#' ## Define globals for model run

#' Leslie matrix:
n_ages <- length(ages) #=17
last_age <- ages[n_ages] #=80+
ffab <- 0.485 # % females at birth (convention used in CBS model).
fertile_ages <- as.character(ages[4:10]) #women give birth between age 15 and 49
fertile_ages_ext <- as.character(ages[3:10]) 
#10-14 age group can give birth because they enter the 15-19 age between two time periods.

#' Housing plans
housing_plans_until <- 2040

#' Housing demand by age
NewHousing_jews20_29 <- 0.32 #Jews 20-29
NewHousing_jews30p <- 0.6 #Jews 30+
NewHousing_arabs20_29 <- 0.21 #Arabs 20-29
NewHousing_arabs30p <- 0.475 #Arabs 20-29

#' Main
proj_length <- length(years[1:6]) # run until 2040
print_progress <- TRUE #should the model print progress during run?

#' Apply fitting to CBS national projections
fit_to_CBS_model <- TRUE

#' ## Prepare base data
base <- merge(proj,base2015)
colnames(base)[1] <- "proj"
base$deaths <- 0
base$year <- years[1]
base$intmig_enter <- 0
base$intmig_exit <- 0
base$yerida <- 0
base$aliya <- 0
base$oldpop <- 0

#' ## Run 
#' Projections adjusted to housing plans
housing_adj_setting <- TRUE 
proj_results <- fn_proj(base, proj_length, fertility, mortality,
                        intmig, aliya, yerida, CBS2065, housing_plans)
pop_results_1 <- proj_results[[1]]
pop_results_1$housing_adjusted <- housing_adj_setting

#' Projections based o historic migration trends
housing_adj_setting <- FALSE 
proj_results <- fn_proj(base, proj_length, fertility, mortality,
                        intmig, aliya, yerida, CBS2065, housing_plans)
pop_results_2 <- proj_results[[1]]
pop_results_2$housing_adjusted <- housing_adj_setting

#' Merge and clean data
pop_results <- rbind(pop_results_1,pop_results_2)
pop_results$nat <- "Jew"
pop_results$nat[pop_results$group %in% groups[4:6]] <- "Arab"
pop_results$nat <- factor(pop_results$nat)
districts$nafa <- factor(districts$nafa, levels=nafas)
pop_results <- left_join(pop_results, districts, by=c("geo"="nafa"))
pop_results$proj <- factor(pop_results$proj, levels=proj)
pop_results$year <- as.numeric(pop_results$year)

#' ## Output
write.csv(pop_results, "output/pop_results.csv", row.names = FALSE)
save(pop_results, file="output/pop_results.Rdata")

#' Prepare for publishing
pop_pub <- pop_results %>% 
  select(proj, housing_adjusted, year, proj, 
         group, district, nafa=geo, gender=sex, age, pop)
pop_pub <- pop_pub %>% 
  mutate(group=recode(group, muslim="arab", christian="arab", 
                      druze="arab", hiloni="jew_not_haredi", dati="jew_not_haredi"))
pop_pub <- pop_pub %>% 
  group_by(proj, housing_adjusted, year, proj, group, 
           district, nafa, gender, age) %>% summarise(pop=sum(pop))

#' Publish data
write.csv(pop_pub, "output/pub_pop_results.csv", row.names = FALSE)



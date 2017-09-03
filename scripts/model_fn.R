#' # Regional Population Projections for Israel - Functions
#' 
#' ### Copyright
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
#' ### Dependencies
#' This program was coded in R 3.4, using RStudio 1.0.143.
#' Additional packages: tidyverse, readxl 
#'
#' ## Leslie Matrix Multiplication Function
#' #### input:
#'- `k0` - population vector (by age). 17 female age groups, 17 male.
#'- `fert` - fertility data, 7 values.
#'- `mort` - mortality data, 36 values: 18 survival rates ($s_x$) for females and males (0-1,1-4,5-9,...,80+).
#' 
#' #### output:
#' `k1` - vector of population by age in next period. 
#' calculated by multiplying $A \cdot k_0$, where $A$ is a structural matrix 
#' with a subdiagonal for aging and deaths, and 2 rows that represent births. 

fn_leslie_mult <- function(k0 ,fert, mort) {
  mort <- arrange(mort, rev(sex),age) #make sure vector is arranged in order (females first)
  fert <- arrange(fert, age)
  sx_f <- mort$rate[mort$sex=="female" & mort$age!=0] #survival rates, not including baby mortality.
  sx_m <- mort$rate[mort$sex=="male" & mort$age!=0]
  s0_f <- mort$rate[mort$age==0 & mort$sex=="female"]
  s0_m <- mort$rate[mort$age==0 & mort$sex=="male"]
  sx_fert <- mort$rate[mort$age %in% fertile_ages_ext & mort$sex=="female"] #survival rates for birthgivers
  
  #creation of matrix (n_ages^2, all zeroes)
  A <- matrix(0,nrow=n_ages*2,n_ages*2)
  colnames(A)<-rep(ages,2)
  rownames(A)<-rep(ages,2)
  
  #subdiagonal - aging and survival rates
  #mort includes age 0-1 that is used for baby mortality. Leslie subdiagonal starts from second object
  A[2:n_ages,1:n_ages-1] <- diag(sx_f[1:(n_ages-1)]) #females
  A[n_ages,n_ages] <- sx_f[n_ages] #for ages 80+ (last value)
  #males
  A[(n_ages+2):(n_ages*2),(n_ages+1):(n_ages*2-1)] <- diag(sx_m[1:(n_ages-1)]) #females
  A[(n_ages*2),(n_ages*2)] <- sx_m[n_ages] #for ages 80+ (last value)
  
  #fertility rates for ages 10-14 to 45-49
  nfx <- c(0,fert$value)
  nfxn <- c(fert$value,0)
  
  # first row formula:
  # (baby mortality)*(asfr per 1 woman, female deaths)*(females at birth)
  # for every age group, half the women are mostly in that age group for the cohort, and half
  # move into the next age group. the formula distributes asfr between two age groups.
  # note: mort[1] = (L1+L0/(5*l0)). multiplication by 5 used for 5 year cohort.
  a1f <- (s0_f*5)*((nfx/1000+nfxn/1000*sx_fert)/2)*(ffab)  
  a1m <- (s0_m*5)*((nfx/1000+nfxn/1000*sx_fert)/2)*(1-ffab)
  A[1,fertile_ages_ext] <- a1f #first row in Leslie for female newborns
  A[n_ages+1,fertile_ages_ext] <- a1m #first row in Leslie for male newborns
  
  #calc k1
  k1 <- as.vector(A %*% k0)
  
  return(k1)
}

#' ## Projection function for natural movements (aging, births and deaths)
#' This function calculates a new cohort in the model w.r.t aging, deaths and births
#' Uses `fn_leslie_mult` for each group in the population.
#'
#' #### Input:
#'- `cohort0` - base pop for all projections. data.frame of proj,group,geo,sex,age,pop (9792 obs)
#'- `t0`,`t1` - base year, next period year. $t_1=t_0+5$
#'- `fertility` - subset of fertility assumptions (asfr) for specific year ($t_1$)
#'- `mortality` - subset of survival assumptions ($s_x$) for specific year ($t_1$)
#'
#' #### Output:
#' `cohort1` - population in $t_1$. data.frame structure same as cohort0, 
#' with additional columns for mortality, fertility, and migration.

fn_proj_natural <- function(cohort0, t0, t1, fertility, mortality) { 
  # Copy data 
  cohort1 <- cohort0
  
  # Subset mortality and fertility data to specific year
  fert_cohort <- filter(fertility, year==t1)
  mort_cohort <- filter(mortality, year==t1)
  
  # arranging and grouping for next steps
  cohort1 <- cohort1 %>% 
    group_by(proj, group,geo) %>% arrange(proj, group,geo,desc(sex),age) #female before male
  mort_cohort <- mort_cohort %>% group_by(proj, group, geo) %>% arrange(proj, group,geo,desc(sex),age)
  fert_cohort <- fert_cohort %>% group_by(proj, group, geo) %>% arrange(proj, group,geo,age)
  
  # Calculate natural movements, using Leslie matrix multiplication function for every group
  cohort1 <- cohort1 %>% mutate(
    pop2=fn_leslie_mult(pop, 
                        fert_cohort[fert_cohort$group==unique(group) & 
                                      fert_cohort$geo==unique(geo) & 
                                      fert_cohort$proj==unique(proj),], 
                        mort_cohort[mort_cohort$group==unique(group) & 
                                      mort_cohort$geo==unique(geo) & 
                                      mort_cohort$proj==unique(proj),]))
  
  # Keep old values, needed for adjustments to housing.
  cohort1 <- cohort1 %>% group_by(proj, group, geo, sex) %>% mutate(oldpop=lag(pop,default = 0))
  
  # Calculate deaths
  cohort1$deaths <- 0
  cohort1$deaths[cohort1$age!=0 & cohort1$age!=last_age] <- 
    cohort1$pop2[cohort1$age!=0 & cohort1$age!=last_age] - 
    cohort1$oldpop[cohort1$age!=0 & cohort1$age!=last_age]
  cohort1$deaths[cohort1$age==last_age] <- 
    cohort1$pop2[cohort1$age==last_age] - 
    (cohort1$pop[cohort1$age==last_age]+cohort1$oldpop[cohort1$age==last_age])
  cohort1 <- cohort1 %>% group_by(proj, group,geo)
  
  # Cleanup
  cohort1$pop <- cohort1$pop2
  cohort1 <- select(cohort1, -pop2)
  return(cohort1)
}

#' ## Projection function for internal migration
#' This function recieves the `cohort1` data frame (population after aging, deaths and births),
#' and calculates the internal migration movements.
#' NOTE: Internal migration changes total population because of differences 
#' in fertility of muslims between different areas.
#' 
#' #### Input:
#'- `cohort1` - Population data before internal migration.
#'- `intmig` - Internal migration assumptions.
#'
#' #### Output:
#' - `cohort1` - Population data after internal migration.

fn_proj_intmig <- function(cohort1, intmig) {
  cohort1_after_intmig <- data.frame()
  for (i in proj){
    cohort1 <- ungroup(cohort1) #clean grouping
    cohort1t <- cohort1 %>% filter(proj==i) %>% select(group, geo, sex, age, pop)
    int_mig_tmp <- merge(intmig, cohort1t, 
                         by.x=c("group","sex","age","exit"), 
                         by.y=c("group","sex","age","geo"), 
                         all.x = TRUE)
    #calc pop*rate, rates are per year. we need 5 year rates.
    int_mig_tmp$mig_abs <- int_mig_tmp$pop*(int_mig_tmp$rate*5) 
    
    #create groups
    exits <- int_mig_tmp %>% group_by(group,sex,age, exit) %>% summarise(pop_exit=sum(mig_abs))
    enters <- int_mig_tmp %>% group_by(group,sex,age, enter) %>% summarise(pop_enter=sum(mig_abs))
    
    #merge data together
    colnames(exits)[4] <- "geo"
    colnames(enters)[4] <- "geo"
    cohort1t<-merge(cohort1t, exits, by=(c("group","sex","age","geo")),all.x=TRUE)
    cohort1t<-merge(cohort1t, enters, by=(c("group","sex","age","geo")), all.x=TRUE)
    rm(exits,enters)
    
    #calc pop after migration
    cohort1t$pop_exit[is.na(cohort1t$pop_exit)] <- 0
    cohort1t$pop_enter[is.na(cohort1t$pop_enter)] <- 0
    cohort1t$pop <- cohort1t$pop + cohort1t$pop_enter - cohort1t$pop_exit
    cohort1t <- cohort1t %>% group_by(group,geo) %>% arrange(group,geo,desc(sex),age)

    # prepare data for return (put into cohort1_after_intmig)
    cohort1t <- ungroup(cohort1t)
    cohort1t$proj <- factor(i, levels=proj)
    cohort1t_to_merge <- cohort1 %>% filter(proj==i) %>% select(-pop, -intmig_enter, -intmig_exit)
    cohort1t_to_merge$proj <- factor(cohort1t_to_merge$proj, levels=proj)
    cohort1t <- left_join(cohort1t, cohort1t_to_merge, by=c("proj", "group", "geo", "sex", "age"))
    cohort1t <- cohort1t %>% select(proj, group, geo, sex, age, pop, deaths, year,
                                    intmig_exit=pop_exit,
                                    intmig_enter=pop_enter,
                                    yerida, aliya, oldpop)
    cohort1t$intmig_exit <- 0 - cohort1t$intmig_exit
    cohort1_after_intmig <- rbind(cohort1_after_intmig, cohort1t)
  }

  return(cohort1_after_intmig)
}

#' ## Projection function for external migration
#' This function adds external (international) migration, divided into immigration and emigration.
#' Based on the CBS assumptions, the total amount of migration is constant. the distribution between
#' geographic regions and population groups is based on the distribution of external migration through
#' the years 2005-2014.
#' 
#' #### Input:
#'- `cohort1` - Population data before internal migration.
#'- `aliya` - External immigration assumptions.
#'- `yerida` - External emmigration assumptions.
#'
#' #### Output:
#'- `cohort1` - Population data after internal migration.

fn_proj_extmig <- function(cohort1, aliya, yerida){
  #Yerida
  cohort1 <- left_join(cohort1,yerida, by=c("group", "geo", "sex", "age"), all=TRUE)
  cohort1$yerida_pop[is.na(cohort1$yerida_pop)]<-0
  cohort1$pop2 <- cohort1$pop + (cohort1$yerida_pop * 5) #emmig is negative, annual
  cohort1$yerida <- cohort1$pop2 - cohort1$pop
  cohort1$pop <- cohort1$pop2
  cohort1 <- select(cohort1, -pop2, -yerida_pop)

  #Aliya
  cohort1 <- merge(cohort1,aliya, by=c("group", "geo", "sex", "age"), all = TRUE)
  cohort1$aliya_pop[is.na(cohort1$aliya_pop)]<-0
  cohort1$pop <- cohort1$pop + (cohort1$aliya_pop * 5) #annual, multiply by 5 for 5-year.
  cohort1$aliya <- cohort1$aliya_pop*5
  cohort1 <- select(cohort1, -aliya_pop)
  
  #Remove migrations of arab population to the West Bank
  cohort1$pop[cohort1$group %in% groups[4:6] & cohort1$geo==nafas[16]] <- 0
  
  return(cohort1)
}

#' ## Main projection function
#' This function calls the functions `fn_proj_natural`, `fn_proj_intmig` and `fn_proj_extmig`
#' to calculate the population change from time $t$ to time $t+5$.
#' 
#'- If `fit_to_CBS_model==TRUE`, then after each period the model is fitted
#' to the national population, given from the CBS 2065 population projections.
#'- If `housing_adj_setting==TRUE`, then housing needs are compared to housing supply, and then 
#' another stage of internal migration is run, in order to disrtibute the population based on the
#' housing plans.
#'- If `print_progress==TRUE`, then progress is printed during run.
#' 
#' #### Input:
#'- `base` - Base population (Population in year 2015)
#'- `proj_length` - Number of periods to run the model.
#'- `fertility`- Fertility assumptions.
#'- `mortality` - Mortality assumptions.
#'- `intmig` - Internal migration assumptions.
#'- `aliya` - External immigration assumptions.
#'- `yerida` - External emmigration assumptions.
#'- `CBS2065` - Population in CBS 2065 projections.
#'- `housing_plans` - Housing supply plans.
#' 
#' #### Output:
#'- `final_results` - Population projection for projection horizon
#'- `housing` - Housing needs (used for analysis)
#'- `housing_adj` - Housing needs after adjustment (used for analysis)

fn_proj <- function(base,proj_length,fertility,mortality,intmig,aliya,yerida, CBS2065, housing_plans){
  final_results <- data.frame(base)
  housing <- data.frame()
  housing_adj <- data.frame() #adjusted housing needs
  
  #run projection
  for (t in years[2:proj_length]) {
    results <- fn_proj_natural(base,t-5,t,fertility, mortality)  # calc natural movements
    results <- fn_proj_intmig(results, intmig)  # calc internal migration
    aliya_cohort <- aliya %>% filter(year==t) %>% select(-year)
    yerida_cohort <- yerida %>% filter(year==t) %>% select(-year)
    if (ext_mig_setting == TRUE) {
      results <- fn_proj_extmig(results, aliya_cohort,yerida_cohort)  # calc external migration    
    }
    results$year <- factor(t, levels=years)
    
    #fit to CBS2065
    if (fit_to_CBS_model==TRUE) results <- fit_model(results, CBS2065)

    #calc new housing needs
    housing_results <- fn_housing_needs(results)
    housing_results$year <- factor(t, levels=years)
    housing <- group_by(housing) #fix grouping
    housing_results <- group_by(housing_results)
    
    housing <- rbind(housing,housing_results)
    
    
    #Housing adjustment
    if (housing_adj_setting==TRUE & t<=housing_plans_until) {
      #adjust population to house building plans  
      results <- fn_adjust_to_plans(results, housing_results,
                                    housing_plans[housing_plans$year==t,],t,intmig)
      
      #calc new housing needs, after adjustments
      housing_results_adj <- fn_housing_needs(results)
      housing_results_adj$year <- factor(t, levels=years)
      housing_adj <- group_by(housing_adj)
      housing_results_adj <- group_by(housing_results_adj)
      housing_adj <- rbind(housing_adj,housing_results_adj)
    }
    
    #prepare for next round
    final_results <- group_by(final_results) #fix grouping
    results <- group_by(results)
    final_results <- rbind(final_results, results)
    base <- results
    
    #print progress
    if (print_progress==TRUE){
      print(paste("Year ", t, " , Housing adjusted: ", housing_adj_setting))
      print(results %>% group_by(proj) %>% summarise(pop=sum(pop)))
    }
    
  }
  
  final_results <- final_results %>% arrange(proj, group,year, geo, sex, age)
  return(list(final_results,housing,housing_adj))
}

#' ## Housing demand calculation
#' This function calculates the empty houses created by immigration and deaths,
#' and the housing needs of new households.
#'
#' #### Input:
#' `cohort` - population data
#'
#' #### Output:
#' `housing` - housing demand

fn_housing_needs <- function(cohort) {

  cohort <- cohort %>% group_by(proj, group, geo,sex,age) %>% arrange(proj, group, geo,sex,age)
  cohort$housing_needs <- 0
  cohort <- left_join(cohort, notmarried, by=c("age", "group", "sex"))
  
  #Jew:
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[1:3] & age==20,
                                        pop*NewHousing_jews20_29))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[1:3] & age==25,
                                        (aliya+yerida+intmig_enter+intmig_exit+deaths)*NewHousing_jews20_29))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[1:3] & age==30,
                                        (0-pop*NewHousing_jews20_29)+pop*NewHousing_jews30p))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[1:3] & age %in% 35:80,
                                        deaths*single+
                                          (intmig_enter+intmig_exit+yerida+aliya)*
                                          NewHousing_jews30p))
  
  #arab
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[4:6] & age==20,
                                        pop*NewHousing_arabs20_29))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[4:6] & age==25,
                                        (aliya+yerida+intmig_enter+intmig_exit+deaths)*
                                          NewHousing_arabs20_29))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[4:6] & age==30,
                                        (0-pop*NewHousing_arabs20_29)+pop*NewHousing_arabs30p))
  cohort <- cohort %>% mutate(housing_needs=
                                replace(housing_needs,
                                        group %in% groups[4:6] & age %in% 35:80,
                                        deaths*single+
                                          (intmig_enter+intmig_exit+yerida+aliya)*
                                          NewHousing_arabs30p))
  
  housing <- cohort %>% group_by(group, geo, proj) %>% 
    select(group, geo, proj, year, housing_needs) %>%
    summarise(housing_needs=sum(housing_needs))
  
  housing$housing_needs <- housing$housing_needs/5  # change to annual values
  return(housing)
  
  
}

#' ## Adjust internal migration to housing plans
#' this is the main function for adjusting migration to plans.
#' It recieves the population and district plans, adjusts the migration to the housing plans
#' and output the adjusted population.  Demand is the housing needs for every population. 
#' It is calculated based on changes in pop in every nafa.
#' supply is the housing plans of the government. I assume all plans will be realized.
#'
#' #### Input:
#'- `cohort` - projection period pop data
#'- `demand` - housing needs
#'- `supply` - housing plans by district
#'- `t` - time
#'- `intmig` - internal migration data
#'
#' #### Output:
#'- `cohort` - cohort population adjusted to plans


fn_adjust_to_plans <- function(cohort,demand, supply, t, intmig) {

  #data is in terms of annual housing units. Transform to 5-year periods:
  demand$housing_needs <- demand$housing_needs*5
  supply$plans <- supply$plans*5
  
  #demand has 288 rows (group*geo*proj). Sum by groups and merge with the supply values
  demand <- demand %>% group_by(geo, proj) %>% summarise(demand=sum(housing_needs)) 
  adjust <- merge(demand, supply, by="geo") %>% select(-year) #merge demand and supply
  
  # Now there is demand and supply for every nafa and proj.
  # The difference between them is the problem
  # If it is positive, then people don't have houses, and need to move.
  # If it is negative, there are empty houses. 
  adjust$diff <- adjust$demand - adjust$plans #Diff>0 is excess demand.
  adjust <- select(adjust, -demand, -plans)
  
  #PART 1
  #Calc how many people need to be removed from excess regions, and where they want to go.
  #####################################################################################
  #these are the places with excess demand. Need to remove people from there:
  excess <- adjust %>% filter(diff>0) #diff>0 there is more demand than supply
  
  # Produce the migration from these places. where people in this place are moving to.
  # excess_mig holds data on values of internal migration (other part in the model).
  # Use this data to calc the rates by which people will move, to solve the excess demand.
  # Merge excess with excess_mig to get the values of migration from the places of excess demand.
  excess_mig <- cohort %>% select(group, geo, sex, age, proj, pop, intmig_exit)
  excess_mig <- as.tbl(excess_mig)
  excess_mig <- left_join(excess, excess_mig, by=c("geo", "proj")) 
  #keep only the places that are in excess demand
  
  # Calc the amount of houses to move out, and the age structure of the out migration.
  # Still don't know where they would have liked to move. 
  
  # Every migration to be broken down into destinations.
  # Using the intmig df, calc "p.rate" - the % of migration by destination for every source of migration
  intmig_dest <- intmig %>% group_by(exit,group,sex,age) %>% 
    mutate(prate=rate/sum(rate)) %>% select(-rate)
  intmig_dest$prate[is.na(intmig_dest$prate)] <-0  # where data is divided by 0.
  
  # Now merging them together:
  excess_mig <- left_join(excess_mig, intmig_dest, by=c("geo"="exit", "group", "sex", "age"))
  
  # Calculate migration flow by destination:
  excess_mig$x_mig <- excess_mig$intmig_exit * excess_mig$prate
  excess_mig <- excess_mig %>% select(-prate, -intmig_exit)
  
  # Now it's known how many people moved from every place to every place (called that x_mig)
  # Using this data to calc stream percentage (calling this y_mig)
  # Grouping is by geo and proj - want the % of (age, sex, group) for every geo and proj
  excess_mig <- excess_mig %>% group_by(geo, proj) %>% mutate(y_mig=x_mig/sum(x_mig))
  excess_mig$y_mig[is.na(excess_mig$y_mig)] <-0
  
  # Add the housing coefs
  excess_mig$housing_coef <- 0
  excess_mig$housing_coef[excess_mig$group %in% groups[1:3] & 
                            excess_mig$age %in% c(20,25)] <- NewHousing_jews20_29
  excess_mig$housing_coef[excess_mig$group %in% groups[1:3] & 
                            excess_mig$age %in% c(seq(30,80,5))] <- NewHousing_jews30p
  excess_mig$housing_coef[excess_mig$group %in% groups[1:3] & 
                            excess_mig$age==30] <- NewHousing_jews30p-NewHousing_jews20_29
  excess_mig$housing_coef[excess_mig$group %in% groups[4:6] & 
                            excess_mig$age %in% c(20,25)] <- NewHousing_arabs20_29
  excess_mig$housing_coef[excess_mig$group %in% groups[4:6] & 
                            excess_mig$age %in% c(seq(30,80,5))] <- NewHousing_arabs30p
  excess_mig$housing_coef[excess_mig$group %in% groups[4:6] & 
                            excess_mig$age==30] <- NewHousing_jews30p-NewHousing_jews20_29
  
  # Calculating total pop to transfer. 
  # tot is calculated by dividing the housing needs by the sum of multiplictions of housing coeficients
  # and the % of migration (by age, sex, group and destination). It gives a population size, that
  # if will be broken down by age, and multiplied by the housing coefs, will fit the difference
  # between demand and supply.
  excess_mig$tot <- 0
  excess_mig <- excess_mig %>% group_by(geo, proj) %>% mutate(tot=diff/c(housing_coef %*% y_mig))
  
  # calculate specific migration
  # This is the amount of people that should be removed, by age, sex and group.
  # it is a result of dividing up tot by all the different y_mig values.
  excess_mig$mig <- -(excess_mig$tot * excess_mig$y_mig) #turn into negative
  excess_mig <- excess_mig %>% select(group, geo, sex, age, proj, enter, mig) #cleaning up
  excess_mig <- ungroup(excess_mig) 
  excess_mig <- arrange(excess_mig)
  #thats it! IT is known from where to take people and where they want to go.
  
  # part 2 - remove people from nafas with excess demand
  ######################################################
  # have to summarize the destinations for each group:
  excess_tot <- excess_mig %>% group_by(group, geo,sex, age, proj) %>% summarize(totmig=sum(mig))
  # merge with cohort, in order to remove people
  cohort <- merge(cohort, excess_tot,by=c("group","geo","sex","age","proj"),all.x=TRUE)
  # for all places that we don't need to remove anyone. there is enough housing
  cohort$totmig[is.na(cohort$totmig)] <- 0 
  
  cohort$pop <- cohort$pop + cohort$totmig #removal of the migration values from pop
  cohort$intmig_exit <- cohort$intmig_exit+cohort$totmig #adding them to the exit_migration flow
  cohort <- select(cohort,-totmig) #cleanup
  adjust$diff[adjust$diff>0] <- 0 #removal the excess demand. 
  
  
  # part 3 - put people where they want to go IF there is place
  ##############################################################
  # Group by destinations, to see the amount of people that want to move to each destination:
  excess_dest <- excess_mig %>% group_by(group,enter,sex,age,proj) %>% 
    summarise(totmig=-sum(mig)) #turn into positive
  excess_dest <- dplyr::rename(excess_dest,geo=enter) #change column name
  
  # calc how much housing will these people need in every nafa
  excess_dest$housing_coef <- 0
  excess_dest$housing_coef[excess_dest$group %in% groups[1:3] & 
                             excess_dest$age %in% c(20,25)] <- NewHousing_jews20_29
  excess_dest$housing_coef[excess_dest$group %in% groups[1:3] & 
                             excess_dest$age %in% c(seq(30,80,5))] <- NewHousing_jews30p
  excess_dest$housing_coef[excess_dest$group %in% groups[1:3] & 
                             excess_dest$age==30] <- NewHousing_jews30p-NewHousing_jews20_29
  excess_dest$housing_coef[excess_dest$group %in% groups[4:6] & 
                             excess_dest$age %in% c(20,25)] <- NewHousing_arabs20_29
  excess_dest$housing_coef[excess_dest$group %in% groups[4:6] & 
                             excess_dest$age %in% c(seq(30,80,5))] <- NewHousing_arabs30p
  excess_dest$housing_coef[excess_dest$group %in% groups[4:6] & 
                             excess_dest$age==30] <- NewHousing_jews30p-NewHousing_jews20_29
  excess_dest$housing_needs <- excess_dest$housing_coef * excess_dest$totmig
  # calc total
  excess_dest_tot <- excess_dest %>% group_by(proj, geo) %>% 
    summarise(housing_needs=sum(housing_needs))
  
  #how much space is there?
  adjust2 <- merge(adjust,excess_dest_tot,by=c("geo","proj"),all=TRUE)
  
  # calc of step 1 of movements - people that want to move to a destination and there is housing available
  # movement subject to minimum. 
  # If there is enough housing, all demand will be supplied.
  # If not, all the available supply will be filled.
  adjust2$step1[-adjust2$diff>adjust2$housing_needs] <- 
    adjust2$housing_needs[-adjust2$diff>adjust2$housing_needs]
  adjust2$step1[-adjust2$diff<adjust2$housing_needs] <- 
    -adjust2$diff[-adjust2$diff<adjust2$housing_needs]
  adjust2$step1[adjust2$diff>0] <- 0 #if there is no place to begin with, don't move anyone.
  
  # calc the % of people that want to move to some place, and there is actually place for them
  # if step1p=1 there is place for all. if step1p=0, there is no place. 
  #if 0<step1p<1, there is partial space
  adjust2$step1p <- adjust2$step1/adjust2$housing_needs  
  adjust3 <- select(adjust2,geo,proj,step1p) 
  
  # merge the % of people that can move in this step with excess_dest
  excess_dest <- merge(excess_dest, adjust3, by=c("geo", "proj"), all=TRUE)
  excess_dest <- select(excess_dest, -housing_needs)
  excess_dest$totmig_step1 <- excess_dest$totmig*excess_dest$step1p #calc streams of people for step 1
  
  # totmig_step1 is the amount of intmig_enter to add to cohort. lets move the people.
  excess_temp <- select(excess_dest,group,geo,sex,age,proj,totmig_step1)
  cohort <- merge(cohort, excess_temp, by=c("group","geo","sex","age","proj"), all.x=TRUE)
  cohort$pop <- cohort$pop+cohort$totmig_step1 #adding people to projection
  cohort$intmig_enter <- cohort$intmig_enter+cohort$totmig_step1 #people were moved!
  rm(excess_temp)
  cohort <- select(cohort,-totmig_step1) 
  
  # update excess_dest: calc how much are left to move, after step 1 was executed
  excess_dest$totmig <- excess_dest$totmig - excess_dest$totmig_step1 #pop left to be moved. 
  excess_dest <- select(excess_dest, -totmig_step1, -step1p) #cleanup
  
  # until this step, it was important where people want to go. 
  # Now, it does not matter, because there is no place for them there.
  # Collapse the destinations in order to prepare for next step of migration
  
  excess_dest <- excess_dest %>% group_by(proj, group,sex,age,housing_coef) %>% 
    summarise(totmig=sum(totmig))
  # Added housing_coef to group_by() to prevent it being deleted
  
  # part 4 - move the rest of people to where there is place.
  # this is the second step in moving the people
  # people will not be moved to where they wanted to go - because there is no place.
  # they will be moved based on the available housing.
  ###########################################################
  
  excess_dest$housing_needs <- excess_dest$housing_coef * excess_dest$totmig 
  #what are the new housing needs?
  #sum the new housing needs
  excess_dest_tot <- excess_dest %>% group_by(proj) %>% summarise(housing_needs=sum(housing_needs))
  #combine with the remaining housing supply
  adjust4 <-adjust2
  adjust4$diff <- adjust4$diff+adjust4$step1
  adjust4 <- select(adjust4,-step1,-step1p,-housing_needs)
  adjust4 <- merge(adjust4,excess_dest_tot,by="proj",all=TRUE)
  
  #for this step, Need two rates:
  # 1. What is the total available housing compared to the total needs.
  adjust4 <- adjust4 %>% group_by(proj) %>% mutate(tot_spc=-sum(diff)/housing_needs)
  # 2. How does the available housing break up between the different areas. 
  adjust4 <- adjust4 %>% group_by(proj) %>% mutate(comp_spc=diff/sum(diff))
  
  #multiplication of two rates will give the composition of destinations in step 2
  adjust4$dest_rate <- adjust4$tot_spc * adjust4$comp_spc
  
  #cleanup
  adjust4 <- adjust4 %>% filter(dest_rate>0) %>% select(proj,geo,dest_rate)
  
  #now there are destinations. I need to merge with excess_dest and calc step2
  excess_dest_step2 <- merge(excess_dest,adjust4,by="proj")
  excess_dest_step2$step2 <- excess_dest_step2$totmig*excess_dest_step2$dest_rate 
  excess_dest_step2 <- excess_dest_step2 %>% select(group,geo,sex,age,proj,step2)
  
  #great. step 2 of migration is ready. lets move people:
  cohort <- merge(cohort, excess_dest_step2, by=c("group", "geo", "sex", "age", "proj"), all.x=TRUE)
  cohort$step2[is.na(cohort$step2)] <- 0 #missing values
  cohort$pop <- cohort$pop+cohort$step2 #adding people!
  cohort$intmig_enter <- cohort$intmig_enter+cohort$step2 
  cohort <- select(cohort, -step2)
  
  #now remove step2 people from excess_dest
  excess_dest_step2 <- excess_dest_step2 %>% group_by(proj, group,age,sex) %>% summarise(step2=sum(step2))
  excess_dest <- merge(excess_dest, excess_dest_step2, by=c("group","sex", "age", "proj"), all.x=TRUE)
  excess_dest$step2[is.na(excess_dest$step2)] <- 0 #just to make sure
  excess_dest$totmig <- excess_dest$totmig - excess_dest$step2
  excess_dest <- select(excess_dest, -step2)
  rm(excess_dest_step2)
  
  #PART 5 - step 3 of migration.
  ###############################
  # Now there is no space for all the rest. 
  # Distribute them based on the size of pop in every nafa (by group)
  excess_dest$housing_needs <- excess_dest$housing_coef * excess_dest$totmig #what are the new housing needs?
  excess_dest_tot2 <- excess_dest %>% group_by(proj) %>% summarise(housing_needs=sum(housing_needs))
  
  # so, need a way to split the rest of the population for migration.
  # Use cohort df to understand the the composition of each nafa by group.
  dest3 <- cohort %>% group_by(proj,group,geo) %>% summarise(pop=sum(pop))
  dest3 <- dest3 %>% group_by(proj, group) %>% mutate(p_pop=pop/sum(pop)) %>% select(-pop)
  
  #merge excess_dest with dest3
  excess_dest <- merge(excess_dest,dest3,by=c("proj", "group"), all.x=TRUE)
  excess_dest$totmig <- excess_dest$totmig*excess_dest$p_pop #split totmig by destination
  excess_dest <- excess_dest %>% select(group, geo, sex, age, proj, totmig)
  rm(dest3)
  
  #merge excess_dest with cohort
  cohort <- merge(cohort, excess_dest, by=c("group", "geo", "sex", "age", "proj"), all.x=TRUE)
  cohort$totmig[is.na(cohort$totmig)] <- 0 #just to make sure
  cohort$pop <- cohort$pop+cohort$totmig
  cohort$intmig_enter <- cohort$intmig_enter+cohort$totmig
  cohort <- select(cohort, -totmig)
  
  #finished. all totmig was distributed.
  
  return(cohort)
  
}

#' ## Fit model aggregates to the CBS National Model
#' After each period in the model, the population is fitted to the total populations of the CBS model.
#' 
#' #### Input:
#'- `results` - cohort data
#'- `CBS2065` - data from CBS
#'
#' #### Output:
#'- `results_fitted` - results after fit to CBS data

fit_model <- function(results, CBS2065){
  results$cbs_group[results$group %in% groups[3]] <- "haredi"
  results$cbs_group[results$group %in% groups[1:2]] <- "jew"
  results$cbs_group[results$group %in% groups[4:6]] <- "arab"
  results$cbs_group <- factor(results$cbs_group)
  
  df1 <- results %>% group_by(year, cbs_group, sex, age, proj) %>% summarise(pop=sum(pop))
  df2 <- left_join(df1, CBS2065, by = c("year", "cbs_group", "sex", "age", "proj")) 
  df2$c <- df2$pop.y/df2$pop.x #create coefs
  df2 <- select(df2, -pop.x, -pop.y)
  results_fitted <- left_join(results, df2,by = c("sex", "age", "proj", "year", "cbs_group"))
  results_fitted$pop <- results_fitted$pop * results_fitted$c
  results_fitted <- select(results_fitted, -c, -cbs_group)
  
  return(results_fitted)
}
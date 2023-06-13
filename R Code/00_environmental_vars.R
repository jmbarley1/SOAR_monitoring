#####SOAR monitoring environmental variables####

require(here)
require(tidyverse)
require(readxl)

#read in master data sheet
data<-read_excel(here('Data','01_SOAR_monitoring_data_master copy.xlsx'),   #Argument "col_types" needed because for some reason, 
                 col_types = c('guess',                                     #R needed to be told what type of data each column is
                               'guess',                                     #And I wanted to skip a couple columns 
                               'guess',
                               'numeric',
                               'numeric',
                               'date',
                               'date',
                               'numeric',
                               'guess',
                               'numeric',
                               'guess',
                               'guess',
                               'numeric',
                               'numeric',
                               'numeric',
                               'numeric',
                               'skip',
                               'numeric',
                               'numeric',
                               'skip',
                               'numeric',
                               'skip',
                               'skip',
                               'skip',
                               'numeric',
                               'skip',
                               'numeric',
                               'skip',
                               'numeric',
                               'skip',
                               'skip',
                               'numeric',
                               'numeric',
                               'skip'))
str(data)

#Changing variables to convenient variable types
data$state<- factor(data$state)
data$collection_agency<- factor(data$collection_agency)
data$site_name<- factor(data$site_name)
data$sampling_method<-factor(data$sampling_method)
data$replicate<- factor(as.character(data$replicate))
data$box_survival<- factor(data$box_survival)
data$SOAR_recruit<- as.numeric(as.character(data$SOAR_recruit))


data<-data %>%                                              #adding needed columns:
  mutate(survival_prop=number_live/total_collected,         #proportional survival of SOAR oysters
         growth= size_after-size_before,                    #Average "growth" of oysters at each site based on the average size of oysters at deployment
         recruit_density= SOAR_recruit*quadrat_size_m2,     #Density of recruitment at each site (those that collected these data)
         live_density_m2= number_live*quadrat_size_m2)      #density of live SOAR oysters


#read in New Hampshire and Massachusetts recruitment data
recruit<-read_excel(here('Data','NH','NH_MA_recruitment.xlsx'))

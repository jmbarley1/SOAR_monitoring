#SOAR monitoring environmental variavbles

require(here)
require(tidyverse)
require(readxl)

data<-read_excel(here('Data','01_SOAR_monitoring_data_master copy.xlsx'), 
                 col_types = c('guess',
                               'guess',
                               'guess',
                               'guess',
                               'guess',
                               'guess',
                               'guess',
                               'skip',
                               'skip',
                               'skip',
                               'skip',
                               'skip',
                               'skip',
                               'guess',
                               'numeric',
                               'numeric',
                               'numeric',
                               'skip',
                               'skip',
                               'numeric',
                               'skip',
                               'numeric',
                               'numeric',
                               'guess',
                               'numeric',
                               'skip',
                               'skip',
                               'guess',
                               'numeric',
                               'guess',
                               'numeric',
                               'numeric',
                               'numeric',
                               'numeric',
                               'guess',
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
data$size_before_error_metric<- factor(data$size_before_error_metric)
data$size_metric<- factor(data$size_metric)
data$substrate<- factor(data$substrate)
data$SOAR_recruit<- as.numeric(as.character(data$SOAR_recruit))


data<-data %>% 
  mutate(survival_prop=number_live/total_collected, 
         growth= size_after-size_before)



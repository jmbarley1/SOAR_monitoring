#EDA on MD monitoring data

require(here)
require(tidyverse)
require(readxl)
data<-read_excel(here('Data','MD','SOAR Monitoring Raw Data.xlsx'), sheet = 'DSOysData')
str(data)
data$Rep<- factor(data$Rep)
data$Substrate<-factor(data$Substrate)
data$OysSize<-factor(data$OysSize)
data$OysterStatus<- factor(data$OysterStatus)

data %>% 
  ggplot(aes(x=Length))+
  geom_histogram()

data$Length

data %>% 
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  ggplot(aes(x=Rep, y=n, fill=OysterStatus))+
  geom_bar(stat = 'identity', position='dodge')

data %>% 
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  group_by(Rep) %>% 
  mutate(prop = n/sum(n)) %>% 
  ggplot(aes(x=Rep, y=prop))+
  geom_point()+
  facet_wrap(~OysterStatus)
  

#same as above, but just with live oysters
data %>% 
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  group_by(Rep) %>% 
  mutate(prop = n/sum(n)) %>% 
  filter(OysterStatus=='Live') %>% 
  ggplot(aes(x=Rep, y=prop))+
  geom_point(size=2)
  
  
  
  
  
  
  
  







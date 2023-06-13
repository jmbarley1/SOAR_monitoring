#EDA on MD monitoring data

#load required libraries
require(here)
require(readxl)
require(plotrix)
require(AICcmodavg)
require(lme4)
require(multcomp)
require(tidyverse)  #have to load tidyverse after multcomp
require(maps)
require(mapdata)
require(raster)

#load environmental variables
source(here('R Code', '00_environmental_vars.R'))

#load color-blind friendly color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colors<-13
cbPalette <- colorRampPalette(c("#56B4E9", "#F0E442"))(colors)

str(data)


#Proportional Survival-----------------------------------------------------------------------------

#Proportional survival by state
data %>%                                    
  filter(!str_detect(site_name, "CO"), 
         survival_prop!="NaN", 
         survival_prop!='NA') %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=survival_prop))+
  geom_histogram()+
  facet_wrap(~state)+
  geom_vline(data=data %>% 
               filter(!str_detect(site_name, "CO"), 
                      survival_prop!="NaN", 
                      survival_prop!='NA') %>% 
               group_by(state) %>% 
               summarize(mean_sur=mean(survival_prop)), aes(xintercept = mean_sur),col='blue',size=1)+
  labs(y='Count', 
       x='Proportion Survival')

#Proportional survival by state and site
data %>%                                    
  filter(!str_detect(site_name, "CO"), 
         survival_prop!="NaN", 
         survival_prop!='NA') %>%
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop, color=state))+
  geom_jitter(alpha=0.6,height=0.09,width=0.3)+
  facet_wrap(~state, scales ='free_x')+
  theme(legend.position = 'none')

#Proportional survival bar graph by state
data %>% 
  filter(!str_detect(site_name, "CO"), 
         survival_prop!="NaN", 
         survival_prop!='NA') %>% 
  group_by(state, collection_agency) %>% 
  summarize(survival_prop=mean(survival_prop*100)) %>% 
  ggplot(aes(x=state, y=survival_prop, fill=collection_agency))+
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_light()+
  theme(legend.position = 'none', 
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 20))+
  labs(x='State', 
       y='Average Percent Survival')+
  scale_fill_manual(values=cbPalette)

#Proportional survival by collection agency and state
data %>%                                    
  #filter(!str_detect(site_name, "CO")) %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=collection_agency, y=survival_prop, color=site_name))+
   geom_jitter(alpha=0.7,position = 
                position_jitterdodge(jitter.width = 1, jitter.height = 0.1, dodge.width = 0.7))+
  facet_wrap(~state, scales ='free_x')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(legend.position = 'none', 
        #axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.x.bottom  =  element_text(size=12))+
  labs(y='Proportaional survival', 
       x='Collection Agency')

#Proportional survival by site name and state
data %>%                                    
  filter(!str_detect(site_name, "CO")) %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop, color=collection_agency))+
  geom_jitter(alpha=0.5,position = 
                position_jitterdodge(jitter.width = 1, jitter.height = 0.1, dodge.width = 0.7), 
              size=2)+
  facet_wrap(~state, scales ='free_x')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x  =  element_blank(), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        strip.text.x = element_text(size = 14))+
  labs(y='Proportaional survival')

#Filtering control sites out of the Maryland dataset
md<-data %>% 
  filter(state=='MD', 
         !str_detect(site_name, "CO"))

str(md)

#Maryland proportional survival by site- scatter plot
md %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=collection_agency, y=survival_prop, color=site_name))+
  geom_jitter(alpha=0.7,position = 
                position_jitter(width = 1, height = 0.1))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  #theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Collection Agency')

#Maryland proportional survival by site- boxplot
md %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=collection_agency, y=survival_prop, color=site_name))+
  geom_boxplot()+
  theme_bw()+
  #theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Collection Agency')

md$site_name<-as.character(md$site_name)

#Maryland proportional survival for each site with mean and standard error bars
#Will spit out a warning that it removed some rows because of missing data. This is just because some rows have 'NA'
md %>% 
  mutate(site=case_when(startsWith(site_name, 'eastern') ~ 'Eastern Bay',                           #changing site names
                        startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                        startsWith(site_name, 'st') ~ 'St. Marys River')) %>%
  group_by(state, collection_agency, site_name, site, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>%                                                  #making sure that we report one average prop. survival per quadrat
  ggplot(aes(x=site, y=survival_prop))+
  geom_jitter(alpha=0.5,position = position_jitter(width = 0.2, height = 0.1), size=2)+             #scatter plot
  geom_pointrange(data= md %>%                                                                      #adding means and error bars
                    mutate(site=case_when(startsWith(site_name,'eastern') ~ 'Eastern Bay', 
                                          startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                                          startsWith(site_name, 'st') ~ 'St. Marys River')) %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency,site) %>% 
                    summarize(survival_prop_mean=mean(survival_prop),
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=site, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error), 
                  color='red')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15))+
  theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Restoration Site')

colors<-3
cbPalette <- colorRampPalette(c("#56B4E9", "#F0E442"))(colors)
range(md$survival_prop)
range(test$survival_prop)

md %>% 
  mutate(site=case_when(startsWith(site_name, 'eastern') ~ 'Eastern Bay',                           #changing site names
                        startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                        startsWith(site_name, 'st') ~ 'St. Marys River'), 
         percent_surv=survival_prop*100) %>%
  group_by(state, collection_agency, site) %>% 
  summarize(percent_surv=mean(percent_surv)) %>%                                                  #making sure that we report one average prop. survival per quadrat
  ggplot(aes(x=site, y=percent_surv, fill=site))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(limits=c(0,100))+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.text = element_text(size=20),
        axis.title  =  element_text(size=26))+
  labs(x='Site',
       y='Oyster Percent Survival')

#MA and NH

data$site_name<-as.character(data$site_name)

#Proportional survival for MA and NH by site with means and standard error bars
#Will spit out a warning that it removed some rows because of missing data. This is just because some rows have 'NA'
#Very similar method to getting this plot to MD
data %>% 
  filter(state=='MA'|
           state=="NH") %>% 
  mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                        startsWith(site_name, 'fair') ~ 'Fairhaven', 
                        startsWith(site_name, 'mart') ~ 'Marthas Vineyard', 
                        startsWith(site_name, 'nan') ~ 'Nannie Island')) %>%
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site, y=survival_prop, color=state))+
  geom_jitter(alpha=0.3,position = position_jitter(width = 0.1, height = 0.1), size=2)+
  geom_pointrange(data= data %>% 
                    filter(state=='MA'|
                             state=='NH') %>% 
                    mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                                          startsWith(site_name, 'fair') ~ 'Fairhaven', 
                                          startsWith(site_name, 'mart') ~ 'Marthas Vineyard', 
                                          startsWith(site_name, 'nan') ~ 'Nannie Island')) %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency, site) %>% 
                    summarize(survival_prop_mean=mean(survival_prop),
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=site, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.title.x = element_text(size=15),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=15), 
        strip.text = element_text(size=15), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=15))+
  labs(y='Proportaional survival', 
       x='Restoration Site')

colors<-4
cbPalette <- colorRampPalette(c("#56B4E9", "#F0E442"))(colors)
data %>% 
  filter(state=='MA'|
           state=="NH") %>% 
  mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                        startsWith(site_name, 'fair') ~ 'Fairhaven', 
                        startsWith(site_name, 'mart') ~ 'Marthas Vineyard', 
                        startsWith(site_name, 'nan') ~ 'Nannie Island'), 
         percent_surv=survival_prop*100) %>%
  drop_na(percent_surv) %>% 
  group_by(state, collection_agency, site) %>% 
  summarize(percent_surv=mean(percent_surv)) %>%                                                  #making sure that we report one average prop. survival per quadrat
  ggplot(aes(x=site, y=percent_surv, fill=site))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(limits=c(0,100))+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.text = element_text(size=20),
        axis.title  =  element_text(size=26))+
  labs(x='Site',
       y='Oyster Percent Survival')
#NY and NJ proportional survival
data$collection_agency<-as.character(data$collection_agency)

#Will spit out a warning that it removed some rows because of missing data. This is just because some rows have 'NA'
#Very similar method to getting this plot to MD
data %>% 
  filter(state=='NY'|
           state=="NJ") %>% 
  mutate(site=factor(case_when(startsWith(collection_agency, 'BOP') ~ 'Governors Island', 
                        startsWith(collection_agency, 'hud') ~ 'Hudson River', 
                        startsWith(collection_agency, 'long') ~ 'Long Island', 
                        startsWith(collection_agency, 'oyster') ~ 'Oyster Bay', 
                        startsWith(collection_agency, 'sedge') ~ 'Sedge Reef', 
                        startsWith(collection_agency, 'baykeeper') ~ 'Earle', 
                        startsWith(collection_agency, 'stock') ~ 'Tuckerton Reef'))) %>%
  drop_na(survival_prop) %>% 
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site, y=survival_prop, color=state))+
  geom_jitter(alpha=0.5,position = position_jitter(width = 0.1, height = 0.1), size=2)+
  geom_pointrange(data= data %>% 
                    filter(state=='NY'|
                             state=='NJ') %>% 
                    drop_na(survival_prop) %>% 
                    mutate(site=fct_relevel(case_when(startsWith(collection_agency, 'BOP') ~ 'Governors Island', 
                                          startsWith(collection_agency, 'hud') ~ 'Hudson River', 
                                          startsWith(collection_agency, 'long') ~ 'Long Island', 
                                          startsWith(collection_agency, 'oyster') ~ 'Oyster Bay', 
                                          startsWith(collection_agency, 'sedge') ~ 'Sedge Reef',
                                          startsWith(collection_agency, 'stock') ~ 'Tuckerton Reef'))) %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency, site) %>% 
                    summarize(survival_prop_mean=mean(survival_prop),
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=site, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.title.x = element_text(size=15),
        axis.text = element_text(size=14),
        axis.title.y = element_text(size=15), 
        strip.text = element_text(size=15), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=15))+
  labs(y='Proportaional survival', 
       x='Restoration Site')

colors<-7
cbPalette <- colorRampPalette(c("#56B4E9", "#F0E442"))(colors)
data %>% 
  filter(state=='NY'|
           state=="NJ") %>% 
  mutate(site=factor(case_when(startsWith(collection_agency, 'BOP') ~ 'Governors Island', 
                               startsWith(collection_agency, 'hud') ~ 'Hudson River', 
                               startsWith(collection_agency, 'long') ~ 'Long Island', 
                               startsWith(collection_agency, 'oyster') ~ 'Oyster Bay', 
                               startsWith(collection_agency, 'sedge') ~ 'Sedge Reef', 
                               startsWith(collection_agency, 'baykeeper') ~ 'Earle', 
                               startsWith(collection_agency, 'stock') ~ 'Tuckerton Reef'), 
                     levels=c('Governors Island','Hudson River', 'Long Island','Oyster Bay', 
                              'Sedge Reef', 'Earle', 'Tuckerton Reef')),
         percent_surv=survival_prop*100)%>%
  drop_na(percent_surv) %>% 
  group_by(state, collection_agency, site) %>% 
  summarize(percent_surv=mean(percent_surv)) %>%                                                  #making sure that we report one average prop. survival per quadrat
  ggplot(aes(x=site, y=percent_surv, fill=site))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(limits=c(0,100))+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.text = element_text(size=20),
        axis.title  =  element_text(size=26))+
  labs(x='Site',
       y='Oyster Percent Survival')

#Just MA proportional survival by site and collection agency
data %>% 
  filter(state=="MA") %>% 
  ggplot(aes(x=site_name, survival_prop))+
  geom_boxplot()+
  geom_jitter(data = data %>% 
                filter(state=='MA') %>% 
                group_by(state, collection_agency, site_name, replicate) %>% 
                summarize(survival_prop=mean(survival_prop)), 
              aes(x=site_name, y=survival_prop, color=site_name), 
               position = position_jitter(width = 0.1, height = 0.1))+
  facet_wrap(~collection_agency, scales = 'free_x')


#NH proportional survival
data %>% 
  filter(state=='NH') %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop))+
  geom_jitter(alpha=0.5,position = position_jitter(width = 0.1, height = 0.1), size=2)+
  geom_pointrange(data = data %>% 
                    filter(state=='NH') %>% 
                    group_by(state, collection_agency, site_name) %>% 
                    summarize(survival_prop_mean=mean(survival_prop),
                              survival_prop_error=sd(survival_prop)),
                  aes(x=site_name, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.text = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_blank())+
  theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Collection Agency')+
  ggtitle('New Hampshire Proportional Survival')


#NY proportional survival by site and collection agency. With mean and standard error bar
data %>% 
  filter(state=='NY') %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop))+
  geom_jitter(alpha=0.7,position = position_jitter(width = 0.2, height = 0.1), size=2)+
  geom_pointrange(data= data %>% 
                    filter(state=='NY', 
                           collection_agency!='BOP') %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency, site_name) %>% 
                    summarize(survival_prop_mean=mean(survival_prop), 
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=site_name, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error), 
                  color='red')+
  facet_wrap(~collection_agency, scales = 'free_x')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  #scale_colour_manual(values=cbPalette)+
  theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=15), 
        strip.text = element_text(size=15))+
  theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Restoration Site')+
  ggtitle('New York Proportional Survival')

#NY proportional survival by collection agency
data %>% 
  filter(state=='NY') %>% 
  group_by(state, collection_agency, site_name) %>%
  drop_na(survival_prop) %>% 
  summarize(survival_prop_mean=mean(survival_prop), 
            survival_prop_error=std.error(survival_prop)) %>%  
  ggplot(aes(x=collection_agency, y=survival_prop_mean, color=collection_agency))+
  geom_pointrange(mapping = aes(ymin = survival_prop_mean-survival_prop_error, 
                                ymax = survival_prop_mean+survival_prop_error), 
                  position = position_jitter(width = 0.2), 
                  size = 0.6)+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=14))+
  labs(y='Proportional Survival', 
       x='Collection Agenccy')


#NJ proportional survival
data %>% 
  filter(state=='NJ') %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=collection_agency, y=survival_prop))+
  geom_jitter(alpha=0.7,position = position_jitter(width = 0.1, height = 0.1), size=2)+
  geom_pointrange(data= data %>% 
                    filter(state=='NJ', 
                           collection_agency!='baykeeper') %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency, site_name) %>% 
                    summarize(survival_prop_mean=mean(survival_prop), 
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=collection_agency, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error), 
                  color='red')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=15), 
        strip.text = element_text(size=15))+
  theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Collection Agency')+
  ggtitle('New Jersey Proportional Survival')

#WA proportional survival
data %>% 
  filter(state=='WA') %>% 
  mutate(site=case_when(startsWith(collection_agency, 'john_adams') ~ 'Baron Point', 
                               startsWith(collection_agency, 'squaxin') ~ 'Squaxin Island')) %>% 
  group_by(state, collection_agency, site_name, site, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site, y=survival_prop))+
  geom_jitter(alpha=0.7,position = position_jitter(width = 0.1, height = 0.1), size=2)+
  geom_pointrange(data= data %>% 
                    filter(state=='WA', 
                           collection_agency!='john_adams') %>% 
                    mutate(site=case_when(startsWith(collection_agency, 'squaxin') ~ 'Squaxin Island')) %>% 
                    drop_na(survival_prop) %>% 
                    group_by(state, collection_agency, site) %>% 
                    summarize(survival_prop_mean=mean(survival_prop), 
                              survival_prop_error=sd(survival_prop)), 
                  aes(x=site, y=survival_prop_mean,
                      ymin = survival_prop_mean-survival_prop_error, 
                      ymax = survival_prop_mean+survival_prop_error), 
                  color='red')+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15), 
        strip.text = element_text(size=15))+
  theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Restoration Site')

colors<-2
cbPalette <- c("#56B4E9", "#009E73")
data %>% 
  filter(state=='WA') %>% 
  mutate(site=case_when(startsWith(collection_agency, 'john_adams') ~ 'Baron Point', 
                        startsWith(collection_agency, 'squaxin') ~ 'Squaxin Island'), 
         percent_surv=survival_prop*100) %>% 
  drop_na(percent_surv) %>% 
  group_by(state, collection_agency, site) %>% 
  summarize(percent_surv=mean(percent_surv)) %>%                                                  #making sure that we report one average prop. survival per quadrat
  ggplot(aes(x=site, y=percent_surv, fill=site))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=cbPalette)+
  scale_y_continuous(limits=c(0,100))+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.text = element_text(size=20),
        axis.title  =  element_text(size=26))+
  labs(x='Site',
       y='Oyster Percent Survival')

#proportional survival for all sites within each state. Color denotes sites
data %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=state, y=survival_prop, color=site_name))+
  geom_jitter(height = .1, width = 0.3)+
  theme(legend.position = 'none')

#looking just at Maryland because it has very little variation

data %>% 
  filter(state=='MD', 
         !str_detect(site_name, "CO")) %>% 
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop, color=site_name, shp=replicate))+
  geom_jitter(alpha=0.5,position=position_jitterdodge(jitter.height = 0.09, jitter.width = 10, dodge.width = 0.05))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = 'none')+
  labs(x='Site Name', 
       y='Proportional Survival')

#looking at replicate- going by state and then collection agency/site
data %>% 
  filter(state=='MD', 
         !str_detect(site_name, "CO")) %>% 
  ggplot(aes(x=replicate, y=survival_prop))+
  geom_jitter()


#Preliminary survival model
#need to first convert survival to 1's and 0's
data<-data %>% 
  mutate(box_survival_mod=case_when(box_survival=="Live" ~ 1, TRUE ~ 0))


survival_mod_list<-list(survival_mod1<- glm(box_survival_mod~state, data=data, family = 'binomial'),
survival_mod2= glm(box_survival_mod~collection_agency, data=data, family = 'binomial'),
survival_mod3= glm(box_survival_mod~site_name, data=data, family = 'binomial'),
survival_mod4= glm(box_survival_mod~state+collection_agency, data=data, family = 'binomial'),
survival_mod5= glm(box_survival_mod~state+site_name, data=data, family = 'binomial'),
survival_mod6= glm(box_survival_mod~collection_agency+site_name, data=data, family = 'binomial'),
survival_mod7= glm(box_survival_mod~state+collection_agency+site_name, data=data, family = 'binomial'))

aictab(survival_mod_list) #mod 3
summary(glm(box_survival_mod~site_name, data=data, family = 'binomial'))
#Tukeys honest difference test- have to use glmer  for this
mod3<-glm(box_survival_mod~site_name, data=data, family = 'binomial')
#summary(glht(mod3, mcp(site_name="Tukey"))) #takes a while. actually, forever

#collinearity-VIF
vif(survival_mod7)# I think that collection agency and site name are highly correlated, therefore throwing an error back
vif(survival_mod5)


#Population estimates-----------------------------------------------------------------------

#MA-saltbox sea farm

#create new MA data set
ma<-data %>% 
  filter(collection_agency=='saltbox_sea_farm') %>% 
  group_by(state, collection_agency, site_name) %>% 
  summarise(live_density_mean=mean(live_density_m2), 
            live_density_error=std.error(live_density_m2, na.rm=TRUE)*quadrat_size_m2)

ma_join<-right_join(data, ma, by="site_name")
ma_join<-ma_join %>% 
  mutate(pop_est=live_density_mean*area_planted_m2) %>% 
  mutate(pop_est_prop=pop_est/number_seeded)

ma_join %>%
  mutate(site=case_when(startsWith(site_name, 'bourne_b') ~ 'Buttermil Bay', 
                        startsWith(site_name, 'bourne_m') ~ 'Monks Cove', 
                        startsWith(site_name, 'fairhaven_r') ~ 'Round Hill',
                        startsWith(site_name, 'fairhaven_s') ~ 'Spindrift',
                        startsWith(site_name, 'fairhaven_t') ~ 'Taylor Culture', 
                        startsWith(site_name, 'fairhaven_w') ~ 'Ward')) %>%#Proportion population est- prop with number seeded
  group_by(state.x, collection_agency.x, site_name, site) %>% 
  summarise(pop_est_prop=mean(pop_est_prop)) %>% 
  ggplot(aes(x=site, y=pop_est_prop))+
  geom_point(width = 0.2, size=2, color='red')+
  geom_point(data=ma_join %>% 
               mutate(site=case_when(startsWith(site_name, 'bourne_b') ~ 'Buttermil Bay', 
                                     startsWith(site_name, 'bourne_m') ~ 'Monks Cove', 
                                     startsWith(site_name, 'fairhaven_r') ~ 'Round Hill', 
                                     startsWith(site_name, 'fairhaven_s') ~ 'Spindrift',
                                     startsWith(site_name, 'fairhaven_t') ~ 'Taylor Culture', 
                                     startsWith(site_name, 'fairhaven_w') ~ 'Ward')) %>% 
               group_by(state.x, collection_agency.x, site) %>%
               drop_na(survival_prop) %>% 
               summarise(pop_est_prop=mean(survival_prop)), 
             mapping = aes(x=site, y=pop_est_prop), size=2, position = 'identity')+
  labs(y='Proportion Population Estimate', 
       x='Site')+
  theme_bw()+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))
  

ma_join %>%         
  group_by(state.y, collection_agency.y, site_name) %>%
  summarise(live_density_mean=mean(live_density_mean), 
            live_density_error=mean(live_density_error)) %>% 
  ggplot(aes(x=site_name, y=live_density_mean, ymin=live_density_mean-live_density_error, 
             ymax=live_density_mean+live_density_error))+
  geom_pointrange()+
  labs(y='Average live oyster density', 
       x='Site')+
  ggtitle('Massachusetts average live oyster density')+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))
  
ma_join %>% 
  group_by(state.x, collection_agency.x, site_name) %>% 
  summarise(number_seeded=mean(number_seeded)) %>% 
  ggplot(aes(x=site_name, y=number_seeded))+
  geom_point(size=2)+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))+
  labs(x='Site', 
       y='Number of oysters seeded')

#MA- Martha's Vineyard
mv<-data %>% 
  filter(collection_agency=='mv_shellfsih_group') %>% 
  group_by(state, collection_agency, site_name) %>% 
  summarise(live_density_mean=mean(live_density_m2), 
            live_density_error=std.error(live_density_m2, na.rm=TRUE)*quadrat_size_m2)

mv_join<-right_join(data, mv, by="site_name")
mv_join<-mv_join %>% 
  mutate(pop_est=(mean(live_density_mean))*area_planted_m2) %>% 
  mutate(pop_est_prop=pop_est/number_seeded)

((mean(mv_join$live_density_mean))*mv_join$area_planted_m2)/16500 ### 2.266


mv_join %>%                                               #Proportion population est- prop with number seeded
  group_by(state.x, collection_agency.x, site_name) %>% 
  summarise(pop_est_prop=mean(pop_est_prop)) %>% 
  ggplot(aes(x=site_name, y=pop_est_prop))+
  geom_jitter(width = 0.2, size=2)+
  labs(y='Proportion Population Estimate', 
       x='Site')+
  ggtitle('Marthas vineyard Population Estimate')+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))

mv_join %>%                                               
  group_by(state.x, collection_agency.x, site_name) %>%
  summarise(live_density_mean=mean(live_density_mean), 
            live_density_error=mean(live_density_error)) %>% 
  ggplot(aes(x=site_name, y=live_density_mean, ymin=live_density_mean-live_density_error, 
             ymax=live_density_mean+live_density_error))+
  geom_pointrange()+
  labs(y='Average live oyster density', 
       x='Site')+
  ggtitle('Marthas Vineyard average live oyster density')+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))

mv_join %>% 
  group_by(state.x, collection_agency.x, site_name) %>% 
  summarise(number_seeded=mean(number_seeded)) %>% 
  ggplot(aes(x=site_name, y=number_seeded))+
  geom_point(size=2)+
  theme(axis.text.x  =  element_text(size = 12), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))+
  labs(x='Site', 
       y='Number of oysters seeded')

#MD
md<-data %>% 
  filter(state=='MD'& !str_detect(site_name, "CO")) %>% 
  group_by(state, collection_agency, site_name) %>% 
  summarise(live_density_mean=mean(live_density_m2), 
            live_density_error=std.error(live_density_m2, na.rm=TRUE)*quadrat_size_m2) 

md_join<-right_join(data, md, by='site_name')
md_join<-md_join %>% 
  mutate(pop_est=live_density_mean*area_planted_m2) %>% 
  mutate(pop_est_prop=pop_est/number_seeded)

md_join %>% 
  mutate(pop_est=live_density_mean*area_planted_m2) %>% 
  group_by(state.x, collection_agency.x, site_name) %>% 
  drop_na(pop_est) %>% 
  summarise(pop_est_prop=mean(pop_est_prop)) %>%
  ggplot(aes(x=site_name, y=pop_est_prop))+
  geom_point(size=2, color='red')+
  geom_point(data=md_join %>% 
                group_by(state.x, collection_agency.x, site_name) %>%
                summarise(pop_est_prop=mean(survival_prop)), mapping = aes(x=site_name, y=pop_est_prop))+
  theme(axis.text.x  =  element_text(size = 12, angle = 45, hjust=1), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))+
  labs(y='Oyster Population Estimate',
       x='Site')+
  ggtitle('Maryland population estimate')
  
md_join %>% 
  group_by(state.x, collection_agency.x, site_name) %>%
  summarise(live_density_mean=mean(live_density_mean), 
            live_density_error=mean(live_density_error)) %>% 
  ggplot(aes(x=site_name, y=live_density_mean, ymin=live_density_mean-live_density_error, 
             ymax=live_density_mean+live_density_error))+
  geom_pointrange()+
  labs(y='Average live oyster density', 
       x='Site')+
  ggtitle('Maryland average live oyster density')+
  theme(axis.text.x  =  element_text(size = 12, angle = 45, hjust=1), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))

md_join %>% 
  group_by(state.x, collection_agency.x, site_name) %>% 
  summarise(number_seeded=mean(number_seeded)) %>% 
  ggplot(aes(x=site_name, y=number_seeded))+
  geom_point(size=2)+
  theme(axis.text.x  =  element_text(size = 12, angle = 45, hjust=1), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))+
  labs(x='Site', 
       y='Number of oysters seeded')

#shinnecock
ny<-data %>% 
  filter(collection_agency=='long_island_shellfish_restoration'|
           collection_agency=='sedge_reef') %>% 
  group_by(state, collection_agency, site_name) %>% 
  drop_na(number_live) %>% 
  summarise(live_density_mean=mean(after_density_live_m2), 
            live_density_error=std.error(live_density_m2, na.rm=TRUE)*quadrat_size_m2) 
ny_join<-right_join(data, ny, by='site_name')
  
ny_join<-ny_join %>% 
  mutate(pop_est=live_density_mean*area_planted_m2) %>% 
  mutate(pop_est_prop=pop_est/number_seeded)

ny_join %>% 
  mutate(pop_est=live_density_mean*area_planted_m2) %>% 
  group_by(state.x, collection_agency.x, site_name) %>% 
  drop_na(pop_est) %>% 
  summarise(pop_est_prop=mean(pop_est_prop)) %>%
  ggplot(aes(x=site_name, y=pop_est_prop))+
  geom_point(size=2, color='red')+
  geom_point(data=ny_join %>% 
               group_by(state.x, collection_agency.x, site_name) %>%
               summarise(pop_est_prop=mean(survival_prop)), mapping = aes(x=site_name, y=pop_est_prop))+
  ylim(0,1)+
  theme(axis.text.x  =  element_text(size = 12, angle = 45, hjust=1), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15), 
        axis.title.x = element_text(size=15))+
  labs(y='Oyster Population Estimate',
       x='Site')+
  ggtitle('Shinnecock Bay population estimate')

#Recruitment---------------------------------------------------------------------------------------
str(data)
which(data$recruit_density>1)
data %>% 
  filter(recruit_density>1) %>% 
  dplyr::select(state, collection_agency, site_name, recruit_density) %>% 
  print(n=300)

data %>% 
  filter(state=='NH',
         recruit_density>0)  #because data were not collected with quadrats

#recruitment density by state
data %>% 
  ggplot(aes(x=state, y=recruit_density))+
  geom_boxplot()

#Just MA and NY
data %>% 
  filter(state %in% c('MA','NY')) %>% 
  ggplot(aes(x=recruit_density, fill=collection_agency))+
  geom_histogram()+
  facet_wrap(~state)

 
data %>% 
  filter(state %in% c('MA','NY')) %>%
  drop_na(recruit_density) %>% 
  ggplot(aes(x=collection_agency, y=recruit_density))+
  geom_boxplot()+
  facet_wrap(~state, scales = 'free_y')

data %>% 
  filter(state %in% c('MA','NY')) %>%
  drop_na(recruit_density) %>% 
  ggplot(aes(x=state, y=recruit_density))+
  geom_bar(stat = 'identity')
  

#MA
data$site_name<-as.character(data$site_name)
data %>% 
  filter(state=='MA') %>% 
  mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                        startsWith(site_name, 'fair') ~ 'Fairhaven', 
                        startsWith(site_name, 'mart') ~ 'Marthas Vineyard')) %>% 
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(recruit_density_mean=mean(recruit_density)) %>% 
  ggplot(aes(x=site, y=recruit_density_mean, color=site))+
  geom_jitter(alpha=0.3,position = position_jitter(width = 0.1, height = 0.1), size=4)+
  geom_pointrange(data=data %>% 
                    filter(state=='MA') %>% 
                    mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                                          startsWith(site_name, 'fair') ~ 'Fairhaven', 
                                          startsWith(site_name, 'mart') ~ 'Marthas Vineyard')) %>% 
                    group_by(state, collection_agency, site) %>% 
                    summarize(recruit_density_mean=mean(recruit_density), 
                              recruit_density_error=sd(recruit_density)), 
                  aes(x=site, y=recruit_density_mean,
                      ymin = recruit_density_mean-recruit_density_error, 
                      ymax = recruit_density_mean+recruit_density_error), 
                  size=1, 
                  linewidth = 1.5)+
  theme_bw()+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=18), 
        legend.position = 'none')+
  labs(x='Restoartion site', 
       y=bquote('Recruitment density'~(m^-2)))
       
colors<-3
cbPalette <- colorRampPalette(c("#F0E442","#56B4E9"))(colors)
data %>% 
  filter(state=='MA') %>% 
  mutate(site=case_when(startsWith(site_name, 'bourne') ~ 'Bourne', 
                        startsWith(site_name, 'fair') ~ 'Fairhaven', 
                        startsWith(site_name, 'mart') ~ 'Marthas Vineyard')) %>% 
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(recruit_density_mean=mean(recruit_density)) %>% 
  ggplot(aes(x=site, y=recruit_density_mean, fill=site))+
  geom_bar(stat = 'identity')+
  theme_bw()+
  theme(axis.text = element_text(size=20), 
        axis.title = element_text(size=26), 
        legend.position = 'none')+
  labs(x='Restoartion site', 
       y=bquote('Recruitment density'~(m^-2)))+
  scale_fill_manual(values = cbPalette)

#NH
data %>% 
  filter(state=='NH') %>% 
  ggplot(aes(x=replicate, y=SOAR_recruit))+
           geom_point()

data %>% 
  filter(state=='NH') %>% 
  ggplot(aes(x=site_name, y=SOAR_recruit))+
  geom_boxplot()+
  geom_jitter(data = data %>% 
                filter(state=='NH') %>% 
                group_by(state, collection_agency, site_name, replicate) %>% 
                summarise(SOAR_recruit=mean(SOAR_recruit)))+
  theme_classic()+
  theme(axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15))+
  xlab('')+
  ylab('Recruitment Count')+
  ggtitle('New Hampshire Recruitment')

recruit$recruitment<-numeric(recruit$recruitment)
str(recruit)

#NH recruitment over time
recruit %>% 
  filter(year>2007) %>% 
  ggplot(aes(x=year, y=recruitment))+
  #geom_point(size=3)+
  theme_bw()+
  geom_line(size=1, color="#0072B2")+
  theme(axis.title = element_text(size = 26), 
        axis.text = element_text(size=20))+
  labs(x='Year', 
       y='Juvenile Oyster Count')

#MD recruitment density
data %>% 
  filter(state=='MD', 
         !str_detect(site_name, "CO")) %>% 
  mutate(site=case_when(startsWith(site_name, 'eastern') ~ 'Eastern Bay', 
                        startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                        startsWith(site_name, 'st') ~ 'St. Marys River')) %>%  
  drop_na(recruit_density) %>% 
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(recruit_density_mean=mean(recruit_density)) %>% 
  ggplot(aes(x=site, y=recruit_density_mean, color=site))+
  geom_jitter(alpha=0.3,position = position_jitter(width = 0.1, height = 0.1), size=4)+
  geom_pointrange(data=data %>% 
                    filter(state=='MD', 
                           !str_detect(site_name, "CO")) %>% 
                    mutate(site=case_when(startsWith(site_name, 'eastern') ~ 'Eastern Bay', 
                                          startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                                          startsWith(site_name, 'st') ~ 'St. Marys River')) %>% 
                    group_by(state, collection_agency, site) %>% 
                    summarize(recruit_density_mean=mean(recruit_density), 
                              recruit_density_error=sd(recruit_density)), 
                  aes(x=site, y=recruit_density_mean,
                      ymin = recruit_density_mean-recruit_density_error, 
                      ymax = recruit_density_mean+recruit_density_error), 
                  size=1, 
                  linewidth = 1.5)+
  scale_y_continuous(limits = c(0, 40), oob = scales::squish)+
  theme_bw()+
  theme(axis.text = element_text(size=16), 
        axis.title = element_text(size=18), 
        legend.position = 'none')+
  labs(x='Restoartion site', 
       y=bquote('Recruitment density'~(m^-2)))

colors<-3
cbPalette <- colorRampPalette(c("#F0E442","#56B4E9"))(colors)
data %>% 
  filter(state=='MD', 
         !str_detect(site_name, "CO")) %>% 
  mutate(site=case_when(startsWith(site_name, 'eastern') ~ 'Eastern Bay', 
                        startsWith(site_name, 'nan') ~ 'Nanticoke River', 
                        startsWith(site_name, 'st') ~ 'St. Marys River')) %>%  
  drop_na(recruit_density) %>% 
  group_by(state, collection_agency, site, site_name, replicate) %>% 
  summarize(recruit_density_mean=mean(recruit_density)) %>% 
  ggplot(aes(x=site, y=recruit_density_mean, fill=site))+
  geom_bar(stat='identity')+
  theme_bw()+
  theme(axis.text = element_text(size=20), 
        axis.title = element_text(size=26), 
        legend.position = 'none')+
  labs(x='Restoartion site', 
       y=bquote('Recruitment density'~(m^-2)))+
  scale_fill_manual(values = cbPalette)


#Size------------------------------------------------------------------------------------------

#Size histograms by state
data %>% 
  filter(!str_detect(site_name, "CO")) %>% 
  ggplot(aes(x=size_after))+
  geom_histogram()+
  geom_vline(data=data %>% 
               filter(!str_detect(site_name, "CO")) %>% 
               drop_na(size_after) %>% 
               group_by(state) %>% 
               summarize(mean_size=mean(size_after)), aes(xintercept = mean_size), col='blue',size=1)+
  facet_wrap(~state)+
  theme_bw()+
  theme(axis.text.x  =  element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        strip.text.x = element_text(size = 14))+
  labs(x='Size after deployment', 
       y='Count')+
  ggtitle('Histogram of oyster size by state')

           
#Growth-------------------------------------------------------------------------
#Growth was calculated by taking the mean of the SOAR oyster deployed at each site and subtracting that from the
#size of each oyster at monitoring. This method proved to be largely unhelpful because it was later figured out 
#that at sites where SOAR oysters were planted on top of already existing reefs, it was hard to distinguish SOAR
#oysters from others during monitoring data collection

#Histogram of growth by state with mean growth denotedd by blue line
data %>% 
  filter(state!='NH') %>% 
  ggplot(aes(x=growth))+
  geom_histogram()+
  geom_vline(data=data %>% 
               filter(!str_detect(site_name, "CO"), state!='NH"') %>% 
               drop_na(growth) %>% 
               group_by(state) %>% 
               summarize(mean_growth=mean(growth)), aes(xintercept = mean_growth),col='blue',size=1)+
  facet_wrap(~state)+
  theme_bw()+
  labs(x='Growth (mm)', 
       y='Count')

#Growth by state, color denotes site
data %>% 
  filter(state!='NH') %>% 
  drop_na(size_before, growth) %>% 
  ggplot(aes(x=collection_agency, y=growth, color=site_name))+
  geom_jitter(position = position_jitterdodge())+
  facet_wrap(~state, scales = 'free_x')+
  theme_bw()+
  theme(legend.position = 'none')+
  geom_hline(yintercept = 0)+
  labs(x='Collection Agency', 
       y='Growth')
  
#looking at growth a different way
#calculating sd from se

data %>% 
  filter(state=="MD") %>% 
  drop_na(size_before) %>% 
  group_by(replicate) %>% 
  summarise(size_before_dist=rnorm(100, mean=size_before, sd=size_before_error)) %>% 
  ggplot(aes(x=size_before_dist,fill='blue'))+
  geom_histogram()+
  geom_histogram(data = data %>% 
                   filter(state=="MD") %>% 
                   drop_na(size_after),
                 aes(x=size_after, fill='red'))+
  facet_wrap(~replicate)
#maybe not the best way to visualize this. lets try strip plot

data %>% 
  filter(state=="MD",) %>% 
  drop_na(size_before) %>% 
  group_by(replicate) %>% 
  summarise(size_before_dist=rnorm(100, mean=size_before, sd=size_before_error), 
            size_after) %>% 
  ggplot(aes(x=replicate, y=size_before_dist), color='blue')+
  geom_jitter(alpha=0.5, position = position_jitter(width = 0.3,), color='blue')+
  geom_jitter(data = data %>% 
                 filter(state=="MD") %>% 
                 drop_na(size_after) %>% 
                 mutate(size_before_dist=size_after), aes(x=replicate, y=size_before_dist, color=replicate), 
                 alpha=0.5, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1))

#Maps---------------------------------------------------------------------------------------------


#practice
#maps package
usa <- map_data("usa")

states<- map_data('state')

ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, fill = region, group = group), color = "white") + 
  coord_fixed(1.3) +
  guides(fill=FALSE)

#east coast

#just NH and MA
new_england <- subset(states, region %in% c("new hampshire", "massachusetts"))
ggplot(data = new_england) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-72,-69.5), ylim = c(41, 43.5))+
  geom_point(data=data %>% 
               filter(state=='NH'|
                        state=='MA'), aes(x=lon, y=lat), color = 'red', size=3)+ 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#just MA
ma<-subset(states, region %in% c('massachusetts'))
ggplot(data = ma) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-72,-69.5), ylim = c(41, 43.5))+
  geom_point(data=data %>% 
               filter(state=='MA'), aes(x=lon, y=lat), color = 'red', size=3)+ 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#just NH
nh<-subset(states, region %in% c('new hampshire'))
ggplot(data = nh) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-72,-69.5), ylim = c(42, 46))+
  geom_point(data=data %>% 
               filter(state=='NH'), aes(x=lon, y=lat), color = 'red', size=3)+ 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#New York and new jersey
new_york<-subset(states, region %in% c('new york', 'new jersey'))
ggplot(data = new_york) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(1.3)+
  geom_point(data=data %>% 
               filter(state=='NY'), aes(x=lon, y=lat), color = 'red', size=3)+
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#Maryland
mid_atlantic<-subset(states, region %in% c('maryland','delaware','virginia'))
ggplot(data = mid_atlantic) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(1.3)+
  geom_point(data=data %>% 
               filter(state=='MD'), aes(x=lon, y=lat), color = 'red', size=3)+
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#Washington
washington<-subset(states, region %in% c('washington'))
ggplot(data = washington) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(1.3)+
  geom_point(data=data %>% 
               filter(state=='WA'), aes(x=lon, y=lat), color = 'red', size=3)+
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank(),
        panel.background = element_blank())

#Early Data Exploration####
#This is early data exploration that I didn't necessarily want to get rid of but is not as polished as the code above

#Maryland###
data<-read_excel(here('Data','MD','SOAR Monitoring Raw Data.xlsx'), sheet = 'DSOysData')
str(data)
data$Rep<- factor(data$Rep)
data$Substrate<-factor(data$Substrate)
data$OysSize<-factor(data$OysSize)
data$OysterStatus<- factor(data$OysterStatus)

data %>% 
  ggplot(aes(x=Length))+
  geom_histogram()+
  theme_classic()

data$Length

data %>% 
  mutate(Rep=case_when(Rep=='001' ~ '1', 
                       Rep=='002' ~ '2',
                       Rep=='003' ~ '3', 
                       Rep=='004' ~ '4', 
                       Rep=='005' ~ '5', 
                       Rep=='006' ~ '6', 
                       Rep=='007' ~ '7', 
                       Rep=='008' ~ '8', 
                       Rep=='009' ~ '9', 
                       Rep=='010' ~ '10')) %>% 
  mutate(Rep = fct_relevel(Rep, '1','2','3','4','5','6','7','8','9','10'),
         OysterStatus = case_when(OysterStatus=='Box' ~ 'Dead',
                                  OysterStatus=='Live' ~ 'Live')) %>%  
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  ggplot(aes(x=Rep, y=n, fill=OysterStatus))+
  geom_bar(stat = 'identity', position='dodge')+
  scale_fill_manual(values = c('#56B4E9','#009E73'))+
  theme_classic()+
  ylab('Count')+
  xlab('Replicate Number')+
  ggtitle('Maryland- all sites')+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18),
        plot.title = element_text(size = 19), 
        legend.position = 'none')
#Survival###
data %>% 
  mutate(Rep=case_when(Rep=='001' ~ '1', 
                       Rep=='002' ~ '2',
                       Rep=='003' ~ '3', 
                       Rep=='004' ~ '4', 
                       Rep=='005' ~ '5', 
                       Rep=='006' ~ '6', 
                       Rep=='007' ~ '7', 
                       Rep=='008' ~ '8', 
                       Rep=='009' ~ '9', 
                       Rep=='010' ~ '10')) %>% 
  mutate(Rep = fct_relevel(Rep, '1','2','3','4','5','6','7','8','9','10')) %>%  
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  group_by(Rep) %>% 
  mutate(prop = n/sum(n)) %>% 
  ggplot(aes(x=Rep, y=prop))+
  geom_point(size=2)+
  facet_wrap(~OysterStatus)+
  theme_bw()+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=18), 
        plot.title = element_text(size=22),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))+
  labs(x='Replicate Number',
       y='Proportion of total oyster sample')
  

#same as above, but just with live oysters
data %>% 
  mutate(Rep=case_when(Rep=='001' ~ '1', 
                       Rep=='002' ~ '2',
                       Rep=='003' ~ '3', 
                       Rep=='004' ~ '4', 
                       Rep=='005' ~ '5', 
                       Rep=='006' ~ '6', 
                       Rep=='007' ~ '7', 
                       Rep=='008' ~ '8', 
                       Rep=='009' ~ '9', 
                       Rep=='010' ~ '10')) %>% 
  mutate(Rep = fct_relevel(Rep, '1','2','3','4','5','6','7','8','9','10')) %>%  
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  group_by(Rep) %>% 
  mutate(prop = n/sum(n)) %>% 
  filter(OysterStatus=='Live') %>% 
  ggplot(aes(x=Rep, y=prop))+
  geom_point(size=2.5)+
  theme_classic()+
  geom_hline(aes(yintercept=mean(prop)), color='red', size=1)+
  geom_hline(aes(yintercept=(mean(prop)+std.error(prop))), lty=2)+
  geom_hline(aes(yintercept=(mean(prop)-std.error(prop))), lty=2)+
  ylab('Proportion of live oysters')+
  xlab('Replicate number')+
  ggtitle('MD survival')+
  labs(caption = 'Proportion of live oysters collected from each replicate bag. Red line denotes the mean across replicates, dashed lines denote SEM.')+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=20))   
  
#same as abovve but with proportion of dead oysters
data %>% 
  mutate(Rep=case_when(Rep=='001' ~ '1', 
                       Rep=='002' ~ '2',
                       Rep=='003' ~ '3', 
                       Rep=='004' ~ '4', 
                       Rep=='005' ~ '5', 
                       Rep=='006' ~ '6', 
                       Rep=='007' ~ '7', 
                       Rep=='008' ~ '8', 
                       Rep=='009' ~ '9', 
                       Rep=='010' ~ '10')) %>% 
  mutate(Rep = fct_relevel(Rep, '1','2','3','4','5','6','7','8','9','10')) %>%  
  group_by(Rep) %>% 
  count(OysterStatus) %>% 
  group_by(Rep) %>% 
  mutate(prop = n/sum(n)) %>% 
  filter(OysterStatus=='Box') %>% 
  ggplot(aes(x=Rep, y=prop))+
  geom_point(size=2.5)+
  theme_classic()+
  geom_hline(aes(yintercept=mean(prop)), color='red', size=1)+
  geom_hline(aes(yintercept=(mean(prop)+std.error(prop))), lty=2)+
  geom_hline(aes(yintercept=(mean(prop)-std.error(prop))), lty=2)+
  ylab('Proportion of dead oysters')+
  xlab('Replicate number')+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=20))
  
  

#Size###
data %>% 
  drop_na(Length) %>% 
  ggplot(aes(x=Length))+
  geom_histogram()+
  geom_vline(aes(xintercept = mean(Length)), color='blue', size=1)+
  ylab('Count')+
  theme_classic()+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=20))
  
#New York###
ny22<-read_excel(here('Data','NY','BOP Governors Is SOAR Oysters.xlsx'), sheet = '2022 data')
str(ny22)  
ny21<-read_excel(here('Data','NY','BOP Governors Is SOAR Oysters.xlsx'), sheet ='pre-deployment')
str(ny21)

#survival
ny22$sample_number<-factor(ny22$sample_number)
#include data for both box and lie oysters per sample

ny22 %>% 
  filter(sample_number=='1') %>%
  select(bag_number, count_live, count_dead) %>% 
  gather(status, count, 2:3) %>% 
  mutate(Status= case_when(status=='count_live' ~ 'Live', 
                           status=='count_dead' ~ 'Dead')) %>% 
  ggplot(aes(x=bag_number, y=count, fill=Status))+
  geom_bar(stat = 'identity', position='dodge')+
  scale_fill_manual(values = c('#56B4E9','#009E73'))+
  labs(x='Replicate Number', 
       y='Count', 
       title = 'New York- Governors Island')+
  theme_classic()+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18),
        plot.title = element_text(size = 19),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))

ny22 %>% 
  filter(sample_number=='1') %>% 
  mutate(prop_live=count_live/25) %>% 
  ggplot(aes(x=bag_number, y=prop_live))+
  geom_bar(stat='identity')+
  labs(x='Replicate Number', 
       y='Proportion of live oysters',
       title = 'Governors Island- Survival')+
  theme_classic()+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size=20))

##just size first
#predeployment

ny21 %>%                                          
  ggplot(aes(x=bag_number, y=shell_height_mm))+
  geom_jitter(height = 1, width=0.3)
 
ny22 %>% 
  ggplot(aes(x=bag_number, y=shell_height_mm))+
  geom_jitter(height=1, width=0.3)

#attempting to calculate mean change in size
ny22_mean<-ny22 %>%                             #After SOAR mean size calculation by bag
  group_by(bag_number) %>% 
  summarise(mean_22=mean(shell_height_mm))

ny21_mean<-ny21%>%                              #before SOAR mean size calculation per bag
  group_by(bag_number) %>% 
  summarise(mean_21=mean(shell_height_mm))

ny_mean<-ny21_mean %>%                          #joining meean size data for both years
  inner_join(ny22_mean, by='bag_number') %>% 
  gather(year, mean, 2:3 ) %>% 
  mutate(year = paste(str_sub(year, -2,-1),bag_number))       #creating a new column that just has the year
  


ny22_sd<-ny22 %>%                               #After SOAR sd of size calculation per bag
  group_by(bag_number) %>% 
  summarise(sd_22=sd(shell_height_mm))

ny21_sd<-ny21 %>%                               #before SOAR sd of size calc per bag
  group_by(bag_number)%>% 
  summarise(sd_21=sd(shell_height_mm))

ny_sd<-ny22_sd %>%                              #joining the sd data for both years
  inner_join(ny21_sd, by='bag_number') %>% 
  gather(year, sd, 2:3) %>% 
  mutate(year = paste(str_sub(year, -2,-1),bag_number))

ny_size<-ny_sd %>%                              #joining the mean and sd data for both years
  full_join(ny_mean, by = c('year')) %>% 
  mutate(bag_number=bag_number.x, 
         year=substr(year, 1, 2)) %>% 
  select(bag_number, year, mean, sd) 
  
ny_size$year<-factor(ny_size$year)
ny_size %>% 
  ggplot(aes(x=bag_number, y=mean, color=year))+
  geom_point(position=position_dodge(width=0.2,preserve = 'total'), size=2.5)+
  geom_linerange(aes(ymax=mean+sd, ymin=mean-sd), position=position_dodge(width = 0.2), size=1.5)+
  scale_color_manual(labels = c('Pre-deployment', 'Post-deployment'), values = c('#56B4E9','#009E73'))+
  theme_classic()+
  ylab('Mean size (mm)')+
  xlab('Replicate number')+
  ggtitle('New York size data- Governors Island')+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18), 
        plot.title = element_text(size=19),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))

#growth
ny_size %>% 
  gather(variable, value, -(bag_number:year)) %>% 
  unite(year_var, year, variable) %>% 
  spread(year_var, value) %>% 
  mutate(growth= `22_mean`-`21_mean`) %>% 
  ggplot(aes(x=bag_number, y=growth))+
  geom_point(size=2.5)+
  labs(x='Replicate Number', 
       y='Mean Growth (mm)', 
       title = 'Governors Island- Growth')+
  theme_classic()+
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18), 
        plot.title = element_text(size=19))



df <- data.frame(month=rep(1:3,2),
                 student=rep(c("Amy", "Bob"), each=3),
                 A=c(9, 7, 6, 8, 6, 9),
                 B=c(6, 7, 8, 5, 6, 7))

df %>% 
  gather(variable, value, -(month:student)) %>%
  unite(temp, student, variable) %>%
  spread(temp, value)

#color-blind friendly pallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ny_size %>% 
  ggplot(aes(x=bag_number, y=mean_size_21, color='red'))+
  geom_jitter(width=0.2, size=2)+
  geom_linerange(aes(ymax=mean_size_21+sd_size_21, ymin=mean_size_21-sd_size_21))+
  geom_jitter(data=ny_size, aes(x=bag_number, y=mean_size_22, color='blue'), size=2, width=0.2)+
  geom_linerange(aes(ymax=mean_size_22+sd_size_22, ymin=mean_size_22-sd_size_22))+
  scale_color_manual(labels = c('After deployment', 'Pre-deployment'), values = c('#E69F00','#009E73'))
  



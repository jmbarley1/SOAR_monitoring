#EDA on MD monitoring data

require(here)
require(readxl)
require(plotrix)
require(AICcmodavg)
require(lme4)
require(multcomp)
require(tidyverse)  #have to load todyverse after multcomp

source(here('R Code', '00_environmental_vars.R'))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
str(data)
levels(data$site_name)

#survival####

data %>%                                    #by population
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


data %>%                                    #by state and site
  filter(!str_detect(site_name, "CO"), 
         survival_prop!="NaN", 
         survival_prop!='NA') %>%
  group_by(state, collection_agency, site_name, replicate) %>% 
  summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=site_name, y=survival_prop, color=state))+
  geom_jitter(alpha=0.6,height=0.09,width=0.3)+
  facet_wrap(~state, scales ='free_x')+
  theme(legend.position = 'none')

data %>%                                    #by collection_agency state
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
        axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.x.bottom  =  element_text(size=12))+
  labs(y='Proportaional survival', 
       x='Collection Agency')

md<-data %>% 
  filter(state=='MD')
str(md)
md %>% 
  filter(!str_detect(site_name, "CO")) %>% 
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
data %>% 
  filter(state=='NH') %>% 
  #group_by(state, collection_agency, site_name, replicate) %>% 
  #summarize(survival_prop=mean(survival_prop)) %>% 
  ggplot(aes(x=collection_agency, y=survival_prop))+
  geom_jitter(alpha=0.7,position = 
                position_jitter(width = 1, height = 0.1))+
  scale_y_continuous(limits = c(0, 1), oob = scales::squish)+
  theme_bw()+
  #theme(legend.position = 'none')+
  labs(y='Proportaional survival', 
       x='Collection Agency')
data %>% 
  


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
#summary(glht(mod3, mcp(site_name="Tukey"))) #takes a while. actually forever

#collinearity-VIF
vif(survival_mod7)# I think that collectiong agengy and site name are highly correlated, therefore throwing and error back
vif(survival_mod5)

#recruitment####
str(data)
which(data$SOAR_recruit>0) #wtf
is.numeric(data$SOAR_recruit)
data %>% 
  ggplot(aes(x=state, y=SOAR_recruit))+
  geom_point()

data %>% 
  ggplot(aes(x=SOAR_recruit))+
  geom_histogram()+
  facet_wrap(~state)




asdfasd
#growth####
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
  filter(state=="MD") %>% 
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

asfsd
#early EDA####
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
  



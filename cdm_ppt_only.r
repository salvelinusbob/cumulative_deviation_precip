
require(tidyverse)
require(lubridate)
require(geoknife)
require(stringr)
require(zoo)


site_lat_dd <- 45.627547 
site_lon_dd <- -89.437777 

############################################
#This section returns the monthly PRISM PPT data for the site_no lat/long 
#using the geoknife package and then combines it with the groundwater level data
############################################


stencil<-simplegeom(c(site_lon_dd,site_lat_dd))
ppt_job<-geoknife(stencil, (list(times = as.POSIXct(c('1895-01-01','2018-01-01')),url = 'http://cida.usgs.gov/thredds/dodsC/prism_v2', variables = 'ppt')), wait = TRUE)
ppt_data = result(ppt_job)[,c(1,2)] 
colnames(ppt_data)<-c("obs_mo", "ppt_mo_obs")
ppt_data$obs_mo<-as.Date(as.POSIXct(ppt_data$obs_mo))



############################################
#select rolling mean length... 60 was optimal in WI
############################################

rolling_mean_length <- 60

ppt_data$obs_mo <- ymd(ppt_data$obs_mo)
ppt_data <- ppt_data %>% 
  mutate(ppt_mo_mean = rollmean(ppt_mo_obs, k=rolling_mean_length, align = "right", fill = FALSE))

ppt_data <- ppt_data %>%
  filter(obs_mo> (paste0(1894+(rolling_mean_length/12),"-12-01")))

ppt_data <- ppt_data %>%
  mutate(ppt_cdm = cumsum(ppt_mo_obs-ppt_mo_mean),
         ppt_cdm_z = ((ppt_cdm-mean(ppt_cdm, na.rm = TRUE))/sd(ppt_cdm, na.rm = TRUE)))


ppt_data %>% 
  ggplot(aes(x=obs_mo))+
  geom_line(aes(y=ppt_cdm_z, color = "CDM"),  size = .5, color ="steel blue")+
  theme(legend.position = c(.15,.9), 
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",
                                         linetype = "solid",
                                         color="light gray"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="light gray"),
        axis.line = element_line(color = "gray"))+
  labs(y = "CDM60 Z-score", 
       x="Observation Month")+ 
  scale_x_date(date_breaks = "10 years", 
               date_labels = "%Y",
               limits = as.Date(c('1906-01-01','2018-01-01')))


#####################################################

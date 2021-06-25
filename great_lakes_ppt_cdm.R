library(reshape2)
library(ggplot2)
library(data.table)
require(tidyverse)
require(lubridate)
require(dataRetrieval)
require(geoknife)
require(stringr)
require(tibbletime)
require(zoo)


Lake_Mi<-fread('https://www.glerl.noaa.gov/data/dashboard/data/levels/1918_PRES/miHuron1918.csv', skip = 2, header=TRUE)
Lake_Mi<-melt(Lake_Mi, id.vars = "year")
Lake_Mi$date<-paste(Lake_Mi$year, Lake_Mi$variable, "01", collapse = NULL)
Lake_Mi$date<-as.Date(Lake_Mi$date, "%Y %B %d")
ggplot(Lake_Mi, aes(Lake_Mi$date, Lake_Mi$value)) +geom_line()

mihursup_huc08 <- read_csv("~/DNR/GIS/DNR Data/Watersheds/mihursup_huc08.csv")
huc08 <- as.vector(mihursup_huc08$HUC8)


#horribly generic and hamfisted area, but a used a bounding box to get prism area within CONUS that was close ot watershed boundary of lake MI

upper_gl_bb <- Polygon(cbind(c(-93,-90,-83,-86),c(47,48,43,41)))
upper_gl_bb <- Polygons(list(upper_gl_bb), "bb1")
upper_gl_bb <- SpatialPolygons(list(upper_gl_bb),
                               proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(upper_gl_bb)
stencil <- simplegeom(upper_gl_bb)
ppt_job<-geoknife(stencil, (list(times = as.POSIXct(c('1895-01-01','2020-01-01')),url = 'http://cida.usgs.gov/thredds/dodsC/prism_v2', variables = 'ppt')), wait = TRUE)
ppt_data = result(ppt_job)[,c(1,2)] 
colnames(ppt_data)<-c("obs_mo", "ppt_mo_obs")
ppt_data$obs_mo<-as.Date(as.POSIXct(ppt_data$obs_mo))


ppt_lvl_obs <- full_join(ppt_data, Lake_Mi, by = c("obs_mo" = "date"))

ppt_lvl_obs <- ppt_lvl_obs  %>% 
  mutate(ppt_mo_mean = rollmean(ppt_mo_obs, k=60, align = "right", fill = FALSE))%>%
  filter(obs_mo> ("1899-12-01"))

ppt_lvl_obs <- ppt_lvl_obs %>%
  mutate(ppt_cdm = cumsum(ppt_mo_obs-ppt_mo_mean),
         ppt_cdm_z = ((ppt_cdm-mean(ppt_cdm, na.rm = TRUE))/sd(ppt_cdm, na.rm = TRUE)),
         lvl_z = ((value -mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)))

ppt_lvl_obs %>%   
  ggplot(aes(x=lvl_z,y=ppt_cdm_z, color = obs_mo))+
  geom_point()+
  theme(legend.position = c(.85,.25), 
        legend.background = element_rect(fill="white",
                                         linetype = "solid",
                                         color="light gray"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="light gray"),
        axis.line = element_line(color = "gray"))+
  labs(y = "Precip 60 mo. CDM (z)", 
       x="Monthly Mean Lake Level (z)", 
       color="Year")

ppt_lvl_obs %>% 
  ggplot(aes(x=obs_mo))+
  geom_line(aes(y=ppt_cdm_z, color = "CDM"),  size = .5)+
  geom_point(aes(y=lvl_z, color = "Water Level"),  size = 1.)+
  scale_color_manual(values= c("steel blue", "dark blue"))+
  theme(legend.position = c(.15,.9), 
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",
                                         linetype = "solid",
                                         color="light gray"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="light gray"),
        axis.line = element_line(color = "gray"))+
  labs(y = "CDM60 and Mean Lake Mi Level Z-score", 
       x="Observation Month")+
  scale_x_date(date_breaks = "10 years", 
               date_labels = "%Y",
               limits = as.Date(c('1906-01-01','2021-01-01')))+
  ggtitle("Cumulative Deviation Monthly Precip vs Lake Michigan/Huron Water Level Z-score")

ppt_lvl_obs <- ppt_lvl_obs %>% 
  mutate(z_diff = ppt_cdm_z - lvl_z)

ppt_lvl_obs %>% 
  ggplot(aes(x=obs_mo))+
  geom_point(aes(y=z_diff, color = "ppt_cdm_z minus lvl_z"),  size = 1.)+
  scale_color_manual(values= c("steel blue", "dark blue"))+
  theme(legend.position = c(.15,.9), 
        legend.title = element_blank(),
        legend.background = element_rect(fill="white",
                                         linetype = "solid",
                                         color="light gray"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="light gray"),
        axis.line = element_line(color = "gray"))+
  labs(y = "CDM60 minus Mean Lake Mi Level Z-score", 
       x="Observation Month")+
  scale_x_date(date_breaks = "10 years", 
               date_labels = "%Y",
               limits = as.Date(c('1906-01-01','2021-01-01')))

ppt_lvl_obs <- ppt_lvl_obs %>% 
  mutate(z_diff = ppt_cdm_z - lvl_z)

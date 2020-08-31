
require(tidyverse)
require(lubridate)
require(dataRetrieval)
require(geoknife)
require(stringr)
require(tibbletime)
require(zoo)


##############################################
#Enter the USGS/NWIS site number to be elavuated here:
rm(list = ls())
site_no <- "05356000"
##############################################
#This section queries NWIS using the dataRetireval package and returns 
#the monthly mean derived from all water level data available
#############################################
sf_data_avail<-whatNWISdata(siteNumber = site_no, parameterCd = c("00060")) 

sf_site_desc<-readNWISsite(c(sf_data_avail$site_no))

# if(any(sf_data_avail$data_type_cd =="gw")){sf_data<-readNWISgwl(sf_data_avail$site_no, convertType = FALSE)} 
# sf_data <- sf_data %>% select(site_no, obs_date = lev_dt, obs = lev_va)
# if(exists("sf_data")){sf_data$obs_date<-as.Date(sf_data$obs_date, "%Y-%m-%d")}
# if(exists("sf_data")) {sf_data$obs <- as.numeric(sf_data$obs)}

if(any(sf_data_avail$data_type_cd =="dv" & sf_data_avail$parm_cd=="00060")){
  dv_sf_data<-readNWISdv(sf_data_avail$site_no, parameterCd =  "00060" , statCd = "00003")}
head(dv_sf_data)

if(exists("dv_sf_data")){ 
  dv_sf_data <- dv_sf_data %>%
    select(site_no, obs_date = Date, obs = X_00060_00003)%>% 
    mutate(obs_mo = ymd(paste0(str_sub(obs_date, 1, 7), "-01")))}

sf_monthly<- dv_sf_data %>% 
  group_by(obs_mo) %>% 
  summarize(site_no = first(site_no),
            sf_mo_mean = mean(obs))%>% 
  mutate(sf_mo_mean = rollmean(sf_mo_mean, k = 12, align = "right", fill = FALSE)) %>% 
  dplyr::select(site_no, obs_mo, sf_mo_mean) 



############################################
#This section returns the monthly PRISM PPT data for the site_no lat/long 
#using the geoknife package and then combines it with the groundwater level data
############################################
stencil<-simplegeom(c(sf_site_desc$dec_long_va, sf_site_desc$dec_lat_va))
ppt_job<-geoknife(stencil, (list(times = as.POSIXct(c('1895-01-01','2020-01-01')),url = 'http://cida.usgs.gov/thredds/dodsC/prism_v2', variables = 'ppt')), wait = TRUE)
ppt_data = result(ppt_job)[,c(1,2)] 
colnames(ppt_data)<-c("obs_mo", "ppt_mo_obs")
ppt_data$obs_mo<-as.Date(as.POSIXct(ppt_data$obs_mo))


ppt_sf_monthly <- left_join(ppt_data, sf_monthly, by =c("obs_mo"))
ppt_sf_monthly$obs_mo <- ymd(ppt_sf_monthly$obs_mo)

#########################################
#This section determines the rolling mean length with the highest correlation
#between precip and groundwater level. The default evaulates 12 month intervals from 1 to 40 years.
#########################################
roll_steps<-seq(12, 480, 12)
cumsumna <- function(x){if(is.null(x)){
  return(x)}
  nas <- is.na(x)
  s <- cumsum(ifelse(nas, 0, x))
  s[nas] <- NA
  return(s)
}

col_names<-map_chr(roll_steps, ~paste0("rollmean_", .x))

rollers<-map(roll_steps, ~rollify(mean, window = .x)) %>% 
  set_names(nm=col_names)
ppt_data_roll <-ppt_data%>%
  filter(obs_mo>"1904-12-01") %>%
  mutate_at("ppt_mo_obs", funs(!!!rollers))
ppt_data_roll<-ppt_data_roll %>% 
  mutate_at(vars(contains("rollmean")),funs(dev=(ppt_mo_obs-.)))
ppt_data_roll<-ppt_data_roll %>%
  arrange(obs_mo) %>%
  mutate_at(vars(contains("dev")),funs(cdm= cumsumna(.)))

ppt_roll <- ppt_data_roll %>% 
  dplyr::select(obs_mo, ends_with("_cdm"))

sf_ppt_roll <- left_join(ppt_sf_monthly, ppt_roll, by = "obs_mo")

roll_corrs <-  sf_ppt_roll %>% 
  summarize_at(vars(ends_with("_cdm")), funs(cor(., sf_mo_mean, use= "pairwise.complete.obs")))

corr_max_mean <-  roll_corrs %>% 
  gather(max_mean_len, corr_max) %>% 
  filter(corr_max == max(corr_max)) %>% 
  mutate(max_mean_len = substr(max_mean_len, 10,12),
         max_mean_len = as.integer(gsub("_", "", max_mean_len)))


########################################
#This section calculates z scores for the optimal precip mean length plots them 
#against groundwater level z-scores. The optimal mean length and correlation between ppt 
#and groundwater level are also returned.
##########################################
ppt_sf_monthly <- ppt_sf_monthly %>% 
  mutate(ppt_mo_mean = rollmean(ppt_mo_obs, k=corr_max_mean$max_mean_len, align = "right", fill = FALSE))

ppt_sf_monthly <- ppt_sf_monthly %>%
  filter(obs_mo> (paste0(1894+(corr_max_mean$max_mean_len/12),"-12-01")))

ppt_sf_monthly <- ppt_sf_monthly %>%
  mutate(ppt_cdm = cumsum(ppt_mo_obs-ppt_mo_mean),
         ppt_cdm_z = ((ppt_cdm-mean(ppt_cdm, na.rm = TRUE))/sd(ppt_cdm, na.rm = TRUE)),
         sf_z = ((sf_mo_mean -mean(sf_mo_mean, na.rm = TRUE))/sd(sf_mo_mean, na.rm = TRUE)),
         sf_mo_meters = sf_mo_mean*.3048)


ppt_sf_monthly %>%   
  ggplot(aes(x=sf_mo_mean,y=ppt_cdm, color = obs_mo))+
  geom_point()+
  theme(legend.position = c(.85,.25), 
        legend.background = element_rect(fill="white",
                                         linetype = "solid",
                                         color="light gray"),
        legend.key = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="light gray"),
        axis.line = element_line(color = "gray"))+
  labs(y = "CDM (mm)", 
       x="Monthly Mean 12 Month Rollling Mean Discharge (cfs)", 
       color="Year", 
       title = paste0("Site:", site_no, "\nMax Mean Length (months) = ", 
                      corr_max_mean$max_mean_len, 
                      "\nPPT/Discharge Pearson R = ",
                      round(corr_max_mean$corr_max, 2)))


ppt_sf_monthly %>% 
  ggplot(aes(x=obs_mo))+
  geom_line(aes(y=ppt_cdm_z, color = "CDM"),  size = .5)+
  geom_point(aes(y=sf_z, color = "Water Level"),  size = 1.)+
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
  labs(y = "CDM60 and 12-Mo Discharge Z-score", 
       x="Observation Month", 
       title = paste0("Site:", site_no, "\nMax Mean Length (months) = ", 
                      corr_max_mean$max_mean_len, 
                      "\nPPT/Discharge Pearson R = ",
                      round(corr_max_mean$corr_max, 2)))+
  scale_x_date(date_breaks = "10 years", 
               date_labels = "%Y",
               limits = as.Date(c('1906-01-01','2018-01-01')))


#####################################################


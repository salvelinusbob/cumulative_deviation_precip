
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
site_no <- 444709089265301
##############################################
#This section queries NWIS using the dataRetireval package and returns 
#the monthly mean derived from all water level data available
#############################################
gw_data_avail<-whatNWISdata(siteNumber = site_no, parameterCd = c("72019")) 

gw_site_desc<-readNWISsite(c(gw_data_avail$site_no))

if(any(gw_data_avail$data_type_cd =="gw")){gw_data<-readNWISgwl(gw_data_avail$site_no, convertType = FALSE)} 
gw_data <- gw_data %>% select(site_no, obs_date = lev_dt, obs = lev_va)
if(exists("gw_data")){gw_data$obs_date<-as.Date(gw_data$obs_date, "%Y-%m-%d")}
if(exists("gw_data")) {gw_data$obs <- as.numeric(gw_data$obs)}

if(any(gw_data_avail$data_type_cd =="dv" & gw_data_avail$parm_cd=="72019")){dv_gw_data<-readNWISdv(gw_data_avail$site_no, parameterCd =  "72019" , statCd = "00001")}

if(exists("dv_gw_data")){ 
  dv_gw_data <- dv_gw_data %>% select(site_no, obs_date = Date, obs = X_72019_00001)}

if(exists("dv_gw_data")){
  gwl_data <- bind_rows(gw_data, dv_gw_data) %>% 
  mutate(obs_mo = ymd(paste0(str_sub(obs_date, 1, 7), "-01")))
}else {
  gwl_data <- gw_data %>% 
    mutate(obs_mo = ymd(paste0(str_sub(obs_date, 1, 7), "-01")))
  }

gwl_monthly<- gwl_data %>% 
  group_by(obs_mo) %>% 
  summarize(gwl_mo_mean = mean(obs)*(-1))

gwl_monthly <- gwl_monthly %>% 
  dplyr::select(site_no, obs_mo, gwl_mo_mean) 



############################################
#This section returns the monthly PRISM PPT data for the site_no lat/long 
#using the geoknife package and then combines it with the groundwater level data
############################################
stencil<-simplegeom(c(gw_site_desc$dec_long_va, gw_site_desc$dec_lat_va))
ppt_job<-geoknife(stencil, (list(times = as.POSIXct(c('1895-01-01','2018-01-01')),url = 'http://cida.usgs.gov/thredds/dodsC/prism_v2', variables = 'ppt')), wait = TRUE)
ppt_data = result(ppt_job)[,c(1,2)] 
colnames(ppt_data)<-c("obs_mo", "ppt_mo_obs")
ppt_data$obs_mo<-as.Date(as.POSIXct(ppt_data$obs_mo))


ppt_gwl_monthly <- left_join(ppt_data, gwl_monthly, by =c("obs_mo"))
ppt_gwl_monthly$obs_mo <- ymd(ppt_gwl_monthly$obs_mo)

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

gwl_ppt_roll <- left_join(ppt_gwl_monthly, ppt_roll, by = "obs_mo")

roll_corrs <-  gwl_ppt_roll %>% 
  summarize_at(vars(ends_with("_cdm")), funs(cor(., gwl_mo_mean, use= "pairwise.complete.obs")))

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
ppt_gwl_monthly <- ppt_gwl_monthly %>% 
  mutate(ppt_mo_mean = rollmean(ppt_mo_obs, k=corr_max_mean$max_mean_len, align = "right", fill = FALSE))

ppt_gwl_monthly <- ppt_gwl_monthly %>%
  filter(obs_mo> (paste0(1894+(corr_max_mean$max_mean_len/12),"-12-01")))

ppt_gwl_monthly <- ppt_gwl_monthly %>%
  mutate(ppt_cdm = cumsum(ppt_mo_obs-ppt_mo_mean),
         ppt_cdm_z = ((ppt_cdm-mean(ppt_cdm, na.rm = TRUE))/sd(ppt_cdm, na.rm = TRUE)),
         gwl_z = ((gwl_mo_mean -mean(gwl_mo_mean, na.rm = TRUE))/sd(gwl_mo_mean, na.rm = TRUE)),
         gwl_mo_meters = gwl_mo_mean*.3048)


ppt_gwl_monthly %>%   
  ggplot(aes(x=gwl_mo_mean,y=ppt_cdm, color = obs_mo))+
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
       x="Monthly Mean Groundwater Level (m)", 
       color="Year", 
       title = paste0("Mean Length (months) = ", 
                      corr_max_mean$max_mean_len, 
                      "\nPPT/GWL Pearson R = ",
                      round(corr_max_mean$corr_max, 2)))


ppt_gwl_monthly %>% 
  ggplot(aes(x=obs_mo))+
  geom_line(aes(y=ppt_cdm_z, color = "CDM"),  size = .5)+
  geom_point(aes(y=gwl_z, color = "Water Level"),  size = 1.)+
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
  labs(y = "CDM60 and Mean Groundwater Level Z-score", 
       x="Observation Month", 
       title = paste0("Mean Length (months) = ", 
                      corr_max_mean$max_mean_len, 
                      "\nPPT/GWL Pearson R = ",
                      round(corr_max_mean$corr_max, 2)))+
  scale_x_date(date_breaks = "10 years", 
               date_labels = "%Y",
               limits = as.Date(c('1906-01-01','2018-01-01')))


#####################################################


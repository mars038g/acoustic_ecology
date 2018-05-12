rm(list = ls())
setwd('c:/users/sean hardison/documents/snap_analyses')

#Libraries
library(dplyr);library(tidyr)
library(lme4);library(nlme)

#Get data
snap_count <- read.csv('snap_count.csv')
snap_long <- gather(snap_count, Site, Val, C1:T3, factor_key=TRUE) #wide to long
snap_long$Recording.Date <- as.Date(snap_long$Recording.Date, format = "%d-%b-%y") #convert date from factor to date class

#add treatment factor
snap_long <- snap_long %>% mutate(Site, Treatment = plyr::mapvalues(Site, from = c("C1","C2","C3","T1","T2","T3"),
                                                        to = c("C","C","C","T","T","T"))) %>%
  group_by(Recording.Date) %>%
  arrange(Recording.Date) %>%
  filter(!is.na(Val))

dates <- data.frame(Recording.Date = seq(as.Date("2016-May-1",format = "%Y-%b-%d"), #dates df to merge on 
             as.Date("2016-Jul-10",format = "%Y-%b-%d"), by = "1 day"))

#merge on dates df and get weeks
snap_long_m <- dates %>% left_join(snap_long, by = "Recording.Date")
snap_long_wks <- snap_long_m %>% mutate(Week = format(Recording.Date, format = "%W")) %>%
  filter(!is.na(Val))

#--------------Add water temperature as a fixed effect----------------#
qual <- read.csv('WaterQual_Summer2016.csv') #Data sampled from CMS dock

#take average by day
temp <- qual %>% group_by(Julian.Day) %>% dplyr::summarise(mean_temp = mean(Temp..deg.C.)) %>%
  filter(Julian.Day <= 192)
plot(mean_temp~Julian.Day, data = temp)

#map mean daily temp to dates df
dates <- cbind(dates, temp)

#merge temperature data with snap count data
snap_long_m <- dates %>% left_join(snap_long, by = "Recording.Date")
snap_long_wks <- snap_long_m %>% mutate(Week = format(Recording.Date, format = "%W")) %>%
  filter(!is.na(Val))

#fit models incorporating temperatures as a fixed effect and time as random. One with AR(1) and one without
mod2_AR <- lme(Val ~ Treatment*mean_temp, random = list(~+1|Site,~+1|Week), 
              data = snap_long_wks, cor = corAR1())
mod2_noAR <- lme(Val ~ Treatment*mean_temp, random = list(~+1|Site,~+1|Week), 
              data = snap_long_wks)



anova(mod2_AR, mod2_noAR) #The model with AR(1) error structure is a significantly better fit

anova(mod2_AR)
summary(mod2_AR)

#Above models DO NOT include nested effects, which is needed to get correct degrees of freedom. 

#--------------Include temperature and random effect with nested structure----------------#
snap_nre <- snap_df %>% 
  mutate(Site, block = plyr::mapvalues(Site, from = c("C1","C2","C3","T1","T2","T3"),
                                       to = c("1","2","3","1","2","3"))) %>%
  mutate(Site, subject = plyr::mapvalues(Site, from = c("C1","C2","C3","T1","T2","T3"),
                                         to = c("1","3","5","2","4","6")))

nested_mod <- lme(Val ~ Treatment*mean_temp, random = ~+1|block/subject,
    data = snap_nre, correlation = corAR1())

anova(nested_mod)

#--------------Fit GLS models (incorrect df)----------------#

#Linear with AR(1) error structure
linear_AR<- gls(snaps ~ temp, 
   data = gls_dat,
   correlation = nlme::corAR1(form = ~1|temp))

newdata <- data.frame(temp = newtemp,
                      temp2 = newtemp^2)

lm_pred <- AICcmodavg::predictSE(linear_AR, 
                                 newdata = newdata,
                                 se.fit = TRUE)
lm_pred <- lm_pred$fit
fit_temp <- seq(min(gls_dat$temp),max(gls_dat$temp),length.out = length(lm_pred))
lines(fit_temp,lm_pred, type = "b")

#Quadratic with AR(1) error structure
gls_dat$temp2 <- gls_dat$temp^2
quadratic_AR<- gls(snaps ~ temp + temp2, 
                data = gls_dat,
                correlation = nlme::corAR1(form = ~1|temp))

newdata <- data.frame(temp = newtemp,
                      temp2 = newtemp^2)

quad_pred <- AICcmodavg::predictSE(quadratic_AR, 
                                 newdata = newdata,
                                 se.fit = TRUE)
quad_pred <- quad_pred$fit
lines(fit_temp,quad_pred, type = "b")

#Compare models 
anova(update(linear_AR, method = "ML"),
      update(quadratic_AR, method = "ML"))

df_aicc <- data.frame(aicc = c(AICc(quadratic_AR),AICc(linear_AR)))

#find significance of model fit - first derive constant model
constant_norm <- nlme::gls(snaps ~ 1, 
            data = gls_dat)
anova(update(constant_norm, method = "ML"),
      update(quadratic_AR, method = "ML"))
anova(update(constant_norm, method = "ML"),
      update(linear_AR, method = "ML"))




rm(list = ls())
setwd('c:/users/sean hardison/documents/snap_analyses')
library(dplyr);library(ggplot2);library(lmerTest)
library(nlme);library(lme4);library(AICcmodavg)

snap_df <- read.csv("Weekly_snap_and_temp.csv")

#visualize snap data with nested effects structure - Correct df

snap_nre <- read.csv("Snap_nested_effects_struc.csv")
nested_mod <- lme(Val ~ mean_temp, random = ~+1|block/subject,
                  data = snap_nre, correlation = corAR1())

temp_pred <- seq(min(snap_nre$mean_temp),max(snap_nre$mean_temp),length.out = length(snap_nre$mean_temp))

#Generate confidence intervals from LME
new.dat <- data.frame(mean_temp = temp_pred)
new.dat$pred <- predict(nested_mod, level = 0, list(mean_temp = temp_pred)) #predict LME

#Compute standard error for predictions
Designmat <- model.matrix(eval(eval(nested_mod$call$fixed)[-2]), new.dat[-ncol(new.dat)])
predvar <- diag(Designmat %*% nested_mod$varFix %*% t(Designmat))
new.dat$SE <- sqrt(predvar) 
new.dat$SE2 <- sqrt(predvar+nested_mod$sigma^2)

#Plot data with 95% CI and prediction intervals

ggplot(new.dat, aes(x = mean_temp, y = pred)) +
  geom_line(aes(x = mean_temp, y = pred)) +
#   geom_ribbon(aes(ymin=pred-2*SE2,ymax=pred+2*SE2),alpha=0.2,fill="red") +
   geom_ribbon(aes(ymin=pred-2*SE,ymax=pred+2*SE),alpha=0.2,fill="blue") +
  geom_point(data = snap_nre, aes(x = mean_temp, y = Val))  +
  ylab("Snap Count") +
  xlab("Mean Temperature (°C)") +
  ylim(-30,350)+
  theme(axis.text=element_text(size=15),
                      text = element_text(size = 15))

#If we were to fit a GLS  
snap_df$temp2 <- snap_df$mean_temp^2
quadratic_AR<- gls(Val ~ mean_temp + temp2, 
                   data = snap_df,
                   correlation = nlme::corAR1(form = ~1|mean_temp))

newdata <- data.frame(mean_temp = snap_df$mean_temp,
                      temp2 =snap_df$mean_temp^2)

quad_pred <- AICcmodavg::predictSE(quadratic_AR, 
                                   newdata = newdata,
                                   se.fit = TRUE)
quad_trend <- quad_pred$fit
quad_lwr <- quad_trend - quad_pred$se.fit*1.96
quad_upr <- quad_pred$se.fit*1.96 + quad_trend


temp_pred <- seq(min(snap_df$mean_temp),max(snap_df$mean_temp),length.out = length(quad_pred))

#snap count and mean temperature
p <- ggplot(snap_df, aes(x=mean_temp, y=Val)) + geom_point() +
  geom_line(aes(y = quad_trend)) + 
  geom_ribbon(aes(ymin=quad_lwr, ymax=quad_upr), alpha=0.2) +
  ylab("Snap Count") +
  xlab("Mean Temperature (°C)") + 
  theme_gray()
p+theme(axis.text=element_text(size=15),
        text = element_text(size = 15))

#Plot snap count by week and treatment group
snap_df <- snap_df %>% mutate(Treatment, Treatment = plyr::mapvalues(Treatment, from = c("C","T"), to = c("NS","S")))
bx <- ggplot(snap_df, aes(as.factor(Week), Val))
bx + geom_boxplot(aes(color = Treatment), lwd = 0.95) +
  theme(axis.text=element_text(size=12),
        text = element_text(size = 12)) + 
  labs(x = "Week", y = "Snap Count")

count <- snap_df %>% group_by(Week,Treatment) %>% dplyr::summarise(count = n())

  

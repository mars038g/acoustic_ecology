## finding RMS SPL of all acoustic samples
library(zoo)
library(seewave)
library(tuneR)
library(Hmisc)

cal = function(wave, f, sens = NULL, lvl = NULL){
  input = inputw(wave = wave, f = f)
  wave = input$w
  f = input$f
  if(!is.null(lvl)){
    if (lvl == 50){
      wave = wave*3.05e-6
      if (!is.null(sens)){
        if (sens == "H1") 
          wave = wave/(10^(-163.8/20))
        else if (sens == "H2") 
          wave = wave/(10^(-164/20)) 
      }
    }
    else if (lvl == 40){
      wave = wave*4.76e-6
      if (!is.null(sens)){
        if (sens == "H1") 
          wave = wave/(10^(-163.8/20))
        else if (sens == "H2") 
          wave = wave/(10^(-164/20))
      }
    }
  }
  rm(input)
  wave = Wave(wave, samp.rate = f, pcm = TRUE, bit = 16)
}

#Wave function - Quick way to remove input from the right channel and return a matrix into S4 format
w <- function(s){
  Wave(s, right = numeric(0), samp.rate = 48000, bit = 16, pcm = TRUE)
}

#example
s1 <- readWave("HC_sound_lvl40.wav")
s1 <- as.matrix(s1@left)
s1 <- w(s1)
s1 <- cal(s1, sens = "H2", lvl = 40)

#wd
setwd("f:/acoustic data")
date_list <- list.files()
date_list <- date_list[2:22]
date_list <- date_list[-c(6,7)]

#H2 <- C, H1 <- T
# for each sampling date
for (i in 1:length(date_list)){
  folder = date_list[i]
  setwd(paste0("f:/acoustic data/",folder,"/Uncut/amp filter edit/"))
  samples <- list.files()
  
  #for each sample in each sampling date
  for (j in 1:length(samples)){
    sample <- readWave(samples[j])
    sample <- as.matrix(sample@left)
    sample <- w(sample)
    
    #if sample recorded at lvl 40 do this:
    if (grepl("40",folder)){
      if (grepl("C",samples[j])){
        sample <- cal(sample, sens = "H2", lvl = 40)
        assign(paste(samples[j],folder), rms(sample@left))
        print(paste(samples[j],"lvl 40"))
        rm(sample)
      } else if (grepl("T",samples[j])){
        sample <- cal(sample, sens = "H1", lvl = 40)
        assign(paste(samples[j],folder), rms(sample@left))
        print(paste(samples[j],"lvl 40"))
        rm(sample)
      }
    
    #otherwise do this:
    } else {
      if (grepl("C",samples[j])){
        sample <- cal(sample, sens = "H2", lvl = 50)
        assign(paste(samples[j],folder), rms(sample@left))
        print(paste(samples[j],"lvl 50"))
        rm(sample)
      } else if (grepl("T",samples[j])){
        sample <- cal(sample, sens = "H1", lvl = 50)
        assign(paste(samples[j],folder), rms(sample@left))
        print(paste(samples[j],"lvl 50"))
        rm(sample)
      }
    }
  }
}

#get values from environment
env <- environment()
env <- as.data.frame(ls(env))
env <- env[grepl(".wav",env$`ls(env)`),]

#add values into data.frame
dat <- data.frame(sample = as.character(env),
                  rms = NA)
for (i in 1:length(env)){
  dat$rms[i] = get(as.character(env[i]))
}

#break up data.frame into components by sampling site
library(dplyr)
library(stringr)
c1 <- dat %>% filter(grepl("CR1",sample)) %>% mutate(Var = "rms c1") %>% dplyr::rename(date = sample)
c2 <- dat %>% filter(grepl("CR2",sample)) %>% mutate(Var = "rms c2") %>% dplyr::rename(date = sample)
c3 <- dat %>% filter(grepl("CR3",sample)) %>% mutate(Var = "rms c3") %>% dplyr::rename(date = sample)

t1 <- dat %>% filter(grepl("TR1",sample)) %>% mutate(Var = "rms t1") %>% dplyr::rename(date = sample)
t2 <- dat %>% filter(grepl("TR2",sample)) %>% mutate(Var = "rms t2") %>% dplyr::rename(date = sample)
t3 <- dat %>% filter(grepl("TR3",sample)) %>% mutate(Var = "rms t3") %>% dplyr::rename(date = sample)

#combine into long format
ndf <- rbind(c1, c2, c3,
             t1, t2, t3)

#extract dates
ndf$date <- str_sub(ndf$date, -7,-1)
ndf$date <- gsub("([A-Za-z]+)"," ",ndf$date)
ndf$date <- gsub("_","-",ndf$date)
ndf$date <- gsub("-40","",ndf$date)
ndf$date <- gsub(" ","", ndf$date)

#convert to dB
ndf_db <- ndf %>% mutate(spl_db = 20*log10(rms), date = as.Date(date, "%m-%d")) %>% group_by(Var, date) 

plot(NULL, xlim = c(as.Date("2018-05-17"),as.Date("2018-07-10")),ylim = c(90,125),
     yaxt = "n", xlab = "Time", ylab = "")


for (i in 1:length(unique(ndf_db$Var))){
  
  data <- ndf_db[ndf_db$Var == unique(ndf_db$Var)[i],]$spl_db
  time <- ndf_db[ndf_db$Var == unique(ndf_db$Var)[i],]$date
  ts <- zoo(data, time)
  points(ts, type = "o", col = i, lwd = 3)
  
}

#get mean and sd by sampling date
ndf_date <- ndf_db %>% group_by(date) %>%
  dplyr::summarise(mean = mean(spl_db),
                   sd = sd(spl_db),
                   n = n())
#break into control and treatment

c_db <- ndf_db %>% filter(grepl("c",Var)) %>%
  group_by(date) %>%
  dplyr::summarise(mean = mean(spl_db),
                   sd = sd(spl_db),
                   n = n()) %>%
  filter(!n < 2)
t_db <- ndf_db %>% filter(grepl("t",Var))%>%
  group_by(date) %>%
  dplyr::summarise(mean = mean(spl_db),
                   sd = sd(spl_db),
                   n = n()) %>%
  filter(!n <= 2)

c_ts <- zoo(c_db$mean, c_db$date)
t_ts <- zoo(t_db$mean, t_db$date)


#plot mean RMS SPL at control and treatment sites
plot(NULL, xlab = "Time" ,xlim = c(as.Date("2018-05-17"),as.Date("2018-07-10")),xaxt = "n",
     ylab = expression(paste("RMS SPL (dB re 1 ",mu,"Pa)")), ylim = c(100, 120), las = 1)
abline(h = mean(ndf_db$spl_db), lty = 2, col = "grey80", lwd = 2)
points(t_ts, type = "o", lwd = 3, col = "darkorange")

errbar(t_db$date, t_db$mean, t_db$mean+t_db$sd,t_db$mean-t_db$sd, add = TRUE, lwd = 2,errbar.col = "darkorange")
points(t_ts, type = "p", pch = 16, lwd = 3, col = "darkorange")

points(c_ts, type = "o", lwd = 3, col = "black", xlab = "Time")
errbar(c_db$date, c_db$mean, c_db$mean+c_db$sd,c_db$mean-c_db$sd, add = TRUE, lwd = 2)

legend("topleft",c("Shrimp Added","No Shrimp Added"),
       col = c("darkorange","black"), cex = 1.35, lwd = 3, bty = "n")
axis(1, at = seq(as.Date("2018-05-17"),as.Date("2018-07-10"),"2 weeks"),
     labels = seq(as.Date("2016-05-17"),as.Date("2016-07-10"),"2 weeks"))

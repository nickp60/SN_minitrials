library(tidyverse)
library(lavaan)
sn_raw <- read.csv("~/GitHub/SN_minitrials/data/clean/trial_data.csv")
sn <- sn_raw %>%
  select(-Methane) %>%
  rowwise() %>% 
  # technically methane per day
  mutate(methane= max(Day.3, Day.6, Day.9) / HRT) %>%
  select(-Day.3, -Day.6, -Day.9)

ggplot(sn, aes(x=Temp, y=Ammonia, color=pH)) + geom_point() 

gsn <- sn %>% 
  gather(key="day", value="methane", Day.3, Day.6, Day.9) %>%
  transform(day = as.numeric(gsub("(.*)\\.(.*)", "\\2", day))) %>%
  mutate(match = day == HRT)
  select(-day) %>%
  distinct()

# get final day value for biogas
sn_raw <- 






tsn <- sn %>% 
  transform(
    E..coli = log(E..coli + 1, 10),
    Coliforms = log(Coliforms + 1, 10),
    Enterococci = log(Enterococci + 1, 10),
    Ammonia = Ammonia/max(Ammonia),
    methane = methane/max(methane),
    Temp = Temp/max(Temp),
    Recipe2 = Recipe2/max(Recipe2),
    sCOD = sCOD/max(sCOD)
  )

summary(tsn)
ggplot(tsn, aes(x=Temp, y=E..coli, color=pH)) + geom_point()

model_1 <- '
  # latent variables
  bugs =~ Coliforms + E..coli + Enterococci
  emissions =~  methane + Ammonia + sCOD
  health =~   sCOD + Recipe2 + methane
  # regressions
  Coliforms ~ OLR + Temp
  E..coli ~ OLR + Temp
  Enterococci ~ OLR + Temp  
  # residual covariances

'

fitA<-sem(model_1, data=tsn)
varTable(fitA)
lavInspect(fitA, "cov.lv")
 
summary(fitA)



HS.model <- ' visual  =~ x1 + x2 + x3 
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
fit
summary(fit, fit.measures=TRUE)

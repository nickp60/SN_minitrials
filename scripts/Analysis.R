library(tidyverse)
library(lavaan)
library(semPlot)
library(GGally)
setwd("~/GitHub/SN_minitrials/scripts/")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

main_data <- read.csv("../data/clean/trial_data.csv")
day0 <- read.csv("../data/clean/day0.csv", skip=3, skipNul = T, stringsAsFactors = F)
day0 <- day0 %>% select(-fs_conc, -fs_vol, -inoc_conc, -inoc_vol) %>% spread(key = measurement, value = Value)
day0<- rbind(
  day0,
  day0 %>% transform(Temp=replace(Temp,Temp==37, 19)),
  day0 %>% transform(Temp=replace(Temp,Temp==37, 55))
) 
names(main_data)
names(day0)
sn_raw <- rbind(main_data, day0)

sn <- sn_raw %>%
  rowwise() %>% 
  transform(
    Recipe=as.factor(Recipe),
    Ecoli =  log10(Ecoli + 1),
    Coliforms = log10(Coliforms + 1),
    Enterococci = log10(Enterococci + 1)) %>%
  transform(
    Ecoli =  ifelse(Ecoli < 1, 1, Ecoli), 
    Coliforms = ifelse(Coliforms < 1, 1, Coliforms), 
    Enterococci =  ifelse(Enterococci < 1, 1, Enterococci))

ggplot(sn, aes(x=VS, y=TS)) + geom_point() + facet_grid(~OL)  + geom_smooth()
ggplot(sn, aes(x=VS, y=TS)) + geom_point() + facet_grid(~Temp) + geom_smooth()
ggplot(sn, aes(x=VS, y=TS, color=Ammonia)) + geom_point() + facet_grid(~Temp) + geom_smooth()
ggplot(sn, aes(x=VS, y=TS, color=Ammonia)) + geom_point() + facet_grid(~Temp + Recipe) + geom_smooth()
ggplot(sn, aes(x=Time)) + geom_histogram() + facet_grid(~OL + Time)

#ggpairs(sn)
ggplot(sn, aes(x=Time, y=Ecoli, color=interaction(OL , Temp))) + 
  geom_smooth(method="lm", formula = "y~x") + 
  facet_wrap(~Temp+ OL, ncol = 4) + 
  geom_jitter()


tsn <- sn %>%
  transform(
    E..coli =   E..coli/max(E..coli),
    Coliforms = Coliforms/max(Coliforms),
    Enterococci = Enterococci/max(Enterococci),
    Ammonia = Ammonia/max(Ammonia),
    methane = methane/max(methane),
    Temp = Temp/max(Temp),
    pH = pH/max(pH),
    OLR = OLR/max(OLR),
    Recipe = as.character(Recipe),
    sCOD = sCOD/max(sCOD),
    TS = TS/max(TS),
    VS = VS/max(VS),
    HRT = HRT/max(HRT)
  )

summary(tsn)
ggplot(tsn, aes(x=Temp, y=E..coli, color=Ammonia)) +
  geom_jitter(width = .01, height = .01)

model_1 <- '
  # latent variables
  bugs =~ Coliforms + E..coli + Enterococci
  health =~ Ammonia + sCOD
  nutrients =~ VS + OLR 
  # regressions
  #TS ~ VS + Recipe + HRT + OLR
  VS ~ Recipe + Temp + HRT + OLR
  pH ~ Temp + Recipe + VS + Ammonia 
  bugs ~ Temp + pH + Ammonia + HRT + Recipe + OLR
  methane ~ Recipe + sCOD + VS + HRT
  Ammonia ~ Recipe + sCOD + VS + HRT 
  # residual covariances
  #Temp ~~ Temp
  #OLR ~~ OLR
  #HRT ~~ HRT
'

fitA<-sem(model_1, data=tsn)
semPaths(fitA, curvePivot = TRUE, title = T)
varTable(fitA)
lavInspect(fitA, "cov.lv")
lavInspect(fitA, "optim.gradient")
lavInspect(fitA)
summary(fitA)


#install.packages("bnlearn")
library(bnlearn)
data(coronary)
bn_df <- data.frame(coronary)
res <- hc(bn_df)
plot(res)

ggplot(sn, aes(y=Enterococci, x=Time)) + 
  facet_wrap(~Temp + OL, ncol=3) +
  geom_jitter(width = .1, height = .1) + geom_smooth(method="lm",formula = "y~x")
ggplot(sn, aes(y=Time, x=OL)) + geom_jitter(width = .1, height = .1) + geom_smooth(method="lm",formula = "y~x")
ggplot(sn, aes(y=Time, x=TS)) + geom_jitter(width = .1, height = .1) + geom_smooth(method="lm",formula = "y~x")
ggplot(sn, aes(y=Time, x=TS)) + geom_jitter(width = .1, height = .1) + geom_smooth(method="lm",formula = "y~x")
ggpairs(sn %>% filter(Time==3), aes(color=Recipe, alpha=.2))
#ggscatmat(sn[, c("Temp", "OLR", "pH", "alpha", "Recipe")],  color="Recipe", alpha = .2)
ggscatmat(sn %>% select(-Temp, -OLR, -alpha),  color="Recipe", alpha = .3)
sn_bn <- sn %>% 
  #select(-TS) %>%
  transform(
    Temp = as.numeric(Temp),
    HRT = as.numeric(HRT),  
    Ammonia = as.numeric(Ammonia), 
    Recipe = as.factor(Recipe))


sn_res <- hc(sn_bn)#, whitelist = data.frame(from=c(), to=c()))
plot(sn_res)
sn_scaled_res <- hc(tsn %>% transform(Recipe=as.factor(Recipe)))
plot(sn_scaled_res)
fittedbn <- bn.fit(sn_scaled_res, data = tsn %>% transform(Recipe=as.factor(Recipe)))
print(fittedbn$Temp)







#install.packages("rstan", dependencies = T)

fit <- stan(
  file = "../models/model_min.stan",
  verbose = T,
  chains=4,
  data = list(
    nbugs=3,
    nrows=nrow(sn),
    Temp=(sn$Temp -19)/18 + 1,
    Ammonia=sn$Ammonia,
    pH=sn$pH,
    Ecoli=sn$Ecoli,
    Coliforms=sn$Coliforms,
    Time=sn$Time,
    OL=sn$OL,
    Enterococci=sn$Enterococci,
    sCOD=sn$sCOD,
    VS=sn$VS,
    TS=sn$TS,
    Recipe=as.numeric(sn$Recipe),
    Methane=sn$Methane
  )
)

plot(fit, pars=c("recipe_effect"))
plot(fit, pars=c("load_effect"))
plot(fit, pars=c("temp_effect"))
plot(fit, pars=c("lp__"))
plot(fit, pars=c("sigma"))
plot(fit, chains=T)

stan_model(fit)
mle = optimizing( stan_model(file = "./models/model_min.stan"), data=c("Methane", "VS"))



# ::::::::::::::::::::: SN minibatch analysis of FIB :::::::::::::::::::::::::::::

# this script is for creating dot plots of delta removal 
# where d removal refers to the difference between required removal
# and achieved removal. 
# required removal = log cfu of FIB in feedstock -

# data used can be found nickp60/SN_minitrials/data/clean/All Raw Data Mini Trial (NW).xlsx
# where FIB are on sheet 2

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
source('setFactorOrder.R') # function for defining factor order; function can be found at end of this script

## === reading in data & calculating required, achieved and delta removal  =======

data = read.table('SN-miniBatchSheet6.txt', header=T, sep='\t')
# or can read in directly with readxl
# library(readxl)
#data = read_excel("All Raw Data Mini Trial (NW).xlsx", sheet=1)

head(data) # here is what the raw data look like: 
# Recipe Temp OLR HRT   pH Coliforms E.coli Enterococci
# 1      2   19 2.0   3 6.98     98000  72700     1990000
# 2      2   55 2.0   3 7.14         0      0        1890
# 3      2   37 0.5   9 7.44      2000   2000      770000
# 4      2   19 0.5   6 7.45    921000 260000     2420000
# 5      2   37 0.5   3 7.40    953000 953000      242000
# 6      2   55 0.5   6 7.45         0      0           0

# read in FS FIB data (found in sheet 3 of SN's "All Raw Data Mini Trial.xlsx")
FSdata = read.table('SN_FIB_FS.txt', header=T, sep='\t')

datay = subset(data, select=-c(pH)) # not doing anything with pH variable so remove
test = merge(datay, FSdata, by="Recipe") # add FS (feedstock) data in by merging

# now create new columns for each FIB with log of cfu/ml number.
test$Colif_log = with(test, log10(Coliforms))
test$Ecoli_log = with(test, log10(E.coli))
test$Entero_log = with(test, log10(Enterococci))

# do the same for feed stock: new columns for each FIB with log of cfu/ml number.
test$Colif_FSlog = with(test, log10(Coliforms_FS))
test$Ecoli_FSlog = with(test, log10(E.coli_FS))
test$Entero_FSlog = with(test, log10(Enterococci_FS))

is.na(test)<-sapply(test, is.infinite) # because some numbers are zero, when logged they
                                        # will become Inf which will bring wrong numbers &
                                        # errors later. I will re-assign these as NA
test[is.na(test)]<-0 # and now convert NAs to 0 

# create new column for each FIB to calculate the 'required removal' which is the 
# number of FIB from each feedstock that need to die in order to reach EU limit of 3log 

test$Colif_required = with(test, Colif_FSlog - 3) # log number of colifs in fs - 3 = no. that need to be removed. 
test$Ecoli_required = with(test, Ecoli_FSlog - 3)
test$Entero_required = with(test, Entero_FSlog - 3)

# create new column for each FIB to calculate the 'achieved removal' which is the 
# number of FIB from each feedstock that died. ie log cfu in feed - log cfu in sample
test$Colif_achieved = with(test, Colif_FSlog - Colif_log)
test$Ecoli_achieved = with(test, Ecoli_FSlog - Ecoli_log)
test$Entero_achieved = with(test, Entero_FSlog - Entero_log)

# create new column for each FIB to calculate the 'delta removal' which is the 
# difference between the required removal and the achieved removal
# in other words. where achieved matched required exactly, delta will be zero.
# a number less than zero means more removal than was required
# a number more than zero means required removal not met. 
test$Colif_dremoval = with(test, Colif_required - Colif_achieved)
test$Ecoli_dremoval = with(test, Ecoli_required - Ecoli_achieved)
test$Entero_dremoval = with(test, Entero_required - Entero_achieved)


## choose only columns that need plotting 
# columns with metadata (1:4) and columns with delta removal for 3 FIBs (23:25)

less = test[,c(1:4,23:25)] # call new small df less
long = melt(less, id.vars = c("Recipe", "OLR", "HRT", "Temp")) # now convert this to long
                                                                # format for plotting
head(long) # this is what it looks like now:
# with default col naames of variable and value

# Recipe OLR HRT Temp       variable     value
# 1      2 2.0   3   19 Colif_dremoval  1.991226
# 2      2 2.0   3   55 Colif_dremoval -3.000000
# 3      2 0.5   9   37 Colif_dremoval  0.301030
# 4      2 0.5   6   19 Colif_dremoval  2.964260
# 5      2 0.5   3   37 Colif_dremoval  2.979093
# 6      2 0.5   6   55 Colif_dremoval -3.000000


library(stringr)
# make a column called 'FIB' which defines which bacteria we mean. 
# we'l use str_extract to extract this from our 'variable' column
long$FIB = str_extract(long$variable, "Colif|Ecoli|Entero")
# make another column called HRT_OLR which is both these variables concatenated (with
# _ inbetween two) to make life easier when plotting subsets of data later on
long$HRT_OLR = paste(long$HRT,long$OLR,sep="_")


## ==== preparations for plotting ============

#  will be facetting by recipe, and want better labels in facet box. 
# I will define the labels i want here. 

FacetLabs= c(
  '2'="Feed 1\n3:1 (FW:S)", 
  '3'="Feed 2\n1:2 (FW:S)", 
  '4'="Feed 3\n1:3 (FW:S)")

## now will do plots by making separate ones for each temperature

# use grep to just get dataset for 19'C -> it will be called long19
long19 = long[grep("19",long$Temp),]

# use set factor order to define the order I want  the FIBs to plotted in
long19[["FIB"]] <- setFactorOrder(long19[["FIB"]], c("Colif", "Ecoli", "Entero"))

# using dplyr to get average and standard dev for each of replicates
# the output will be put in a df called plotData

plotData = long19 %>%
  group_by(Recipe, HRT_OLR,FIB) %>% # group data by these 3 variables
  summarise(MEAN = mean(value), # then calculate mean
            SD = sd(value))     # and sd of the reps within the groupings


## ==== plot 'template' for 19'C  ============

# create a plot 'template' called myplot that will then be used to plot subsets
# of data where data is subset by HRT_OLR

library(viridis) # will use greyscale compatible palette virids.
library(scales)

myplot= ggplot(plotData,aes(x=FIB, y=MEAN, ymin=MEAN-SD, ymax=MEAN+SD, fill=FIB))+ # using values calc in lines 133 - 136
  geom_pointrange(size=0.75)+ # i want points with error bars aka geom_pointrange
  geom_hline(yintercept = 0, linetype=2) + # will draw a line at delta = 0 for clarity
  geom_point(colour="black", shape=21, size = 7) + # define outline colour and size of points wehre shape =21 is a fillable shape
  facet_wrap(~Recipe, labeller = as_labeller(FacetLabs)) + # facet wrap by recipe
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 5), limits=c(-4,4.5)) +
  # first plot without limits() so you know spread of data & can define correct limits.
  scale_fill_viridis(discrete=TRUE,name="", 
                     labels=c("Coliforms", "E. coli", "Enterococci"))+
  scale_x_discrete(limits = rev(levels(plotData$FIB)), # do reverse as i will flip axes
                   labels=c("Coliforms", expression(italic("E. coli")), "Enterococci"))+
  coord_flip() + theme_bw()  + # flip axes & define theme
  theme(axis.text.y=element_text(size=12, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title.x=element_text(size=12, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

myplot

# call plot for each HRT_OLR subset in 19'C dataset. eg for this one its 3_2
# which means 3 days hrt and 2 g olr. label as appropriate
p1 = myplot %+% subset(plotData, HRT_OLR %in% c("3_2")) + 
  labs(x="", title="3 days at 2g OLR", y="")
p1 # look at plot to check

## do same with all other subsets of data, giving new name to plot each time.
p2 = myplot %+% subset(plotData, HRT_OLR %in% c("6_0.5")) + 
  labs(x="", title="6 days at 0.5g OLR", y="")
p2

p3 = myplot %+% subset(plotData, HRT_OLR %in% c("6_3.5")) + 
  labs(x="", title="6 days at 3.5g OLR", y=expression(Delta*" removal (log cfu/ml)"))
p3

p4 = myplot %+% subset(plotData, HRT_OLR %in% c("9_2")) + 
  labs(x="", title="9 days at 2g OLR", y=expression(Delta*" removal (log cfu/ml)"))
p4

## ==== Arrange and annotate plots  for 19'C  ============

# using ggpubr

deg19plots = ggarrange(p1, NULL, p2, p3, p4,NULL, # order i want plots to appear where null is empty space
                       ncol = 2, nrow = 3,  align = "hv", # 2 columns and 3 rows
                       widths = c(1, 1), heights = c(1, 1), # size of plots
                       labels = c("A", "","B", "C", "D", ""),
                       legend = "none") # legend a bit unnecessary as y labs define col

deg19plots # check it
## now add common title
deg19title = annotate_figure(deg19plots,
                             top = text_grob(expression(paste(" ", Delta, " removal at 19", degree, "C",  " (required - achieved FIB removal)")),
                                             color = "black", face = "bold", size = 17))
# now save as pdf (use preview to double check best size for plot)
ggsave("SN_minibatch_19C_dRemoval.pdf",deg19title, width=12, height=8.5, units="in")

## ==== Repeat the same for 37'C data ============


long37 = long[grep("37",long$Temp),] # select data
# define order
long37[["FIB"]] <- setFactorOrder(long37[["FIB"]], c("Colif", "Ecoli", "Entero"))

# preparing data for plot by grouping and averaging etc
plotData37 = long37 %>%
  group_by(Recipe, HRT_OLR,FIB) %>%
  summarise(MEAN = mean(value),
            SD = sd(value))

# define plot and make any appropriate changes of limits etc
myplot= ggplot(plotData37,aes(x=FIB, y=MEAN, ymin=MEAN-SD, ymax=MEAN+SD, fill=FIB))+
  geom_pointrange(size=0.75)+
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(colour="black", shape=21, size = 7) +
  facet_wrap(~Recipe, labeller = as_labeller(FacetLabs)) +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 5), limits=c(-4.0,4.5)) +
  scale_fill_viridis(discrete=TRUE,name="", 
                     labels=c("Coliforms", "E. coli", "Enterococci"))+
  scale_x_discrete(limits = rev(levels(plotData$FIB)), 
                   labels=c("Coliforms", expression(italic("E. coli")), "Enterococci"))+
  coord_flip() + theme_bw()  +
  theme(axis.text.y=element_text(size=12, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title.x=element_text(size=12, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5))
myplot

# call plot for each HRT_OLR subset in 37C dataset (note these are different
# to subsets of data in 19C plot due to box benkin (sp) design of expt)
q1 = myplot %+% subset(plotData37, HRT_OLR %in% c("3_0.5")) + 
  labs(x="", title="3 days at 0.5g OLR", y="")
q1

q2 = myplot %+% subset(plotData37, HRT_OLR %in% c("3_3.5")) + 
  labs(x="", title="3 days at 2.5g OLR", y="")
q2

q3 = myplot %+% subset(plotData37, HRT_OLR %in% c("6_2")) + 
  labs(x="", title="6 days at 2g OLR", y="")
q3

q4 = myplot %+% subset(plotData37, HRT_OLR %in% c("9_0.5")) + 
  labs(x="", title="9 days at 0.5g OLR", y=expression(Delta*" removal (log cfu/ml))"))
q4

q5 = myplot %+% subset(plotData37, HRT_OLR %in% c("9_3.5")) + 
  labs(x="", title="9 days at 3.5g OLR", y=expression(Delta*" removal (log cfu/ml))"))
q5

# arrange plots together as one
deg37plots = ggarrange(q1, q2, q3, NULL, q4, q5,
                       ncol = 2, nrow = 3,  align = "hv", 
                       widths = c(1, 1), heights = c(1, 1),
                       labels = c("A", "B","C", "", "D", "E"),
                       legend="none")

# add common title
deg37title = annotate_figure(deg37plots,
                             top = text_grob(expression(paste(" ", Delta, " removal at 37", degree, "C",  " (required - achieved FIB removal)")),
                                             color = "black", face = "bold", size = 17))
# write out .
ggsave("SN_minibatch_37C_dRemoval.pdf",deg37title, width=12, height=8.5, units="in")

### ==== and lastly for 55 degrees  ===

long55 = long[grep("55",long$Temp),]

# averaging and sd this 55C subset
plotData55 = long55 %>%
  group_by(Recipe, HRT_OLR,FIB) %>%
  summarise(MEAN = mean(value),
            SD = sd(value))

# define plot hcange ylims as appropriate. I kept as c(-4,4.5) for ALL plots
myplot= ggplot(plotData55,aes(x=FIB, y=MEAN, ymin=MEAN-SD, ymax=MEAN+SD, fill=FIB))+
  geom_pointrange(size=0.75)+
  geom_hline(yintercept = 0, linetype=2) +
  geom_point(colour="black", shape=21, size = 7) +
  facet_wrap(~Recipe, labeller = as_labeller(FacetLabs)) +
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 5), limits=c(-4,4.5)) +
  scale_fill_viridis(discrete=TRUE,name="", 
                     labels=c("Coliforms", "E. coli", "Enterococci"))+
  scale_x_discrete(limits = rev(levels(plotData$FIB)), 
                   labels=c("Coliforms", expression(italic("E. coli")), "Enterococci"))+
  coord_flip() + theme_bw()  +
  theme(axis.text.y=element_text(size=12, colour="black"),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title.x=element_text(size=12, colour="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

myplot

# do for subsets of HRT_OLR. This time they are same subsets as for 19C data. 
r1 = myplot %+% subset(plotData55, HRT_OLR %in% c("3_2")) + 
  labs(x="", title="3 days at 2g OLR", y="")
r1

r2 = myplot %+% subset(plotData55, HRT_OLR %in% c("6_0.5")) + 
  labs(x="", title="6 days at 0.5g OLR", y="")
r2

r3 = myplot %+% subset(plotData55, HRT_OLR %in% c("6_3.5")) + 
  labs(x="", title="6 days at 3.5g OLR", y=expression(Delta*" removal (log cfu/ml)"))
r3

r4 = myplot %+% subset(plotData55, HRT_OLR %in% c("9_2")) + 
  labs(x="", title="9 days at 2g OLR", y=expression(Delta*" removal (log cfu/ml)"))
r4

# arrange
deg55plots = ggarrange(r1, NULL, r2, r3, r4,NULL,
                       ncol = 2, nrow = 3,  align = "hv", 
                       widths = c(1, 1), heights = c(1, 1),
                       labels = c("A", "","B", "C", "D", ""),
                       legend="none")

deg55plots 
## adding title 

deg55title = annotate_figure(deg55plots,
                             top = text_grob(expression(paste(" ", Delta, " removal at 55", degree, "C",  " (required - achieved FIB removal)")),
                                             color = "black", face = "bold", size = 17))

# write out

ggsave("SN_minibatch_55C_dRemoval.pdf",deg55title, width=12, height=8.5, units="in")


## ====== end ======= 


setFactorOrder <- function(x, order=sort(levels(x))) { 
  # Returns a factor ordered by `order`.  
  # If order is missing, defaults to `levels(x)` if available, else to `sort(unique(x))`
  # Useful for ggplot and elsewhere were ordering is based on the order of the levels
  
  if (!is.factor(x)) {
    warning("`x` is not a factor. Will coerce.")
    levs <- sort(unique(x))
    if (missing(order))
      order <- levs
  } else {
    levs <- levels(x)
  }
  
  # any values in order, not in levels(x)
  NotInx <- setdiff(order, levs)
  
  if (length(NotInx)) {
    warning ("Some values not in x:\n", paste(NotInx, collapse=", "))
  }
  
  # levels(x) not explicitly named in order
  Remaining <-  setdiff(levs, order)
  
  order <- c(setdiff(order, NotInx), Remaining)
  
  factor(x, level=order)
}


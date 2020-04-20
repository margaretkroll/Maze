#This assumes data collected in IbexFarm, https://github.com/addrummond/ibex, using the Maze Ibex implementation
#found at https://github.com/vboyce. However, the analysis and plotting are more widely applicable.

#The bayesian analysis uses brms, a wrapper for STAN: https://github.com/paul-buerkner/brms

setwd("yourwd")
library(ggplot2)
library(ggthemes)
library(ordinal)
library(lme4)
library(magrittr)
library(lmerTest)
library(tidyr)
library(plotrix)
library(reshape)
library(stringr)
library(dplyr)
library(optimx)
library(dfoptim)
library(brms)
library(rstan) 
rstan_options(auto_write=TRUE)
library(rstanarm)
options(mc.cores=parallel::detectCores())



#read in data
allData <- read.csv('IbexResults.csv', header = 0, comment.char = "#")

#using IP address as Participant #, change this if you're running in a lab
colnames(allData) <- c("Time","Participant","Controller", "Number", "Element", "Condition", "Item", "Seq", "Word","Foil","Button", "Correct", "RT", "FullSent")

#remove any participants you need to 
allData <- subset(allData, allData$Participant != "Participant#")

allData %>% 
  group_by(Participant) %>% 
  summarize(n=n()) -> participant_Sum

write.csv(allData,'allData.csv')


#Participant check

#pull all exp. items and fillers (using example conditions NL, ND, EL, and ED)
doi_all <- subset(allData, allData$Condition == "NL" | allData$Condition == "ND" | allData$Condition == "EL" | allData$Condition == "ED" | allData$Condition == "F")

doi_all %>% 
  group_by(Participant) %>% 
  summarize(n=n()) -> participant_Sum_doi

#pull only the maze data (remove questions)
maze_all <- filter(doi_all, Controller == "Maze")


#Exlude participants by total sentences missed
#subset into only trials for which the participants made it through the whole sentence
#set the variable to identify a period in a string
period <- "[.]"

#use str_detect to label the cells with a period in the word column (this will identify only the final words in each trial)
#my items here were not the same length, if all the sentence are the same length just subset by the final region
maze_all$final <- str_detect(maze_all$Word, period)

#add a new column 'TrialSucc' such that 1 = the participant successfully completed the final word of the sentence in the trial
maze_all$TrialSucc <- if_else(maze_all$final == TRUE & maze_all$Correct == "yes", 1, 0)

#pull only final trials
maze_finalTrials <- subset(maze_all, maze_all$final == TRUE)

maze_finalTrials %>% 
  group_by(Participant, Correct) %>% 
  summarize(n=n()) -> participantsAllMaze

write.csv(participantsAllMaze,'participantAllMaze.csv')

#use spread to reshape data to get average correct by participant
participantsAllMaze <- spread(participantsAllMaze, Correct, n)
participantsAllMaze$yes[is.na(participantsAllMaze$yes)] <- 0
participantsAllMaze$no[is.na(participantsAllMaze$no)] <- 0
participantsAllMaze$total <- participantsAllMaze$no + participantsAllMaze$yes
participantsAllMaze$percCorr <- participantsAllMaze$yes / participantsAllMaze$total

#can exclude extra participants here or in the line above


#Clean Maze RT data

#pull exp. items maze only
doi_plaus <- subset(allData, allData$Condition == "NL" | allData$Condition == "ND" | allData$Condition == "EL" | allData$Condition == "ED")

#pull only Maze trials
doi_plaus <- subset(doi_plaus, doi_plaus$Controller == "Maze")

#pull only the question data (remove maze)
question_all <- filter(doi_plaus, Controller == "Question")
question_all %>% 
  group_by(Participant, Foil) %>% 
  summarize(n=n()) -> participantsQuestions

#subset the participants by Correct (remove trials with incorrect response)
usableDOI <- subset(doi_plaus, doi_plaus$Correct == "yes")

#make the RTs numeric
usableDOI$RT <- as.numeric(as.character(usableDOI$RT))
is.numeric(usableDOI$RT)

#remove sequence 0, since it's not informative
usableDOI <- subset(usableDOI, usableDOI$Seq != 0)

#then trim off the top 1% of the RTs
#this trim cut-off is a judgment call
upperTrim2 <- quantile(x = usableDOI$RT, probs = .997)
upperTrim2
usableDOI <- subset(usableDOI, usableDOI$RT < 4323.987)

#trim fast end: anything below 200ms because these are physically bounded by how long it takes to push a button
usableDOI <- subset(usableDOI, usableDOI$RT > 199)


#add in log transformed RTs if you want them
#log transform instead of z-score because z-score doesn't actually make your data normal;
#it just transforms around the mean, which is skewed for RTs bc RTs are not normally distributed
#or use an exgaussian in your analysis (below)
usableDOI[,"RT_log"] <- NA
usableDOI$RT_log <- log2(usableDOI$RT)


#label critical regions
#this is only for data whose critical regions don't line up with the Maze trial Seq output
usableDOI$Region <- if_else(usableDOI$Seq == 8 & (usableDOI$Condition == "ND" | usableDOI$Condition == "NL"), 1, usableDOI$Region)
usableDOI$Region <- if_else(usableDOI$Seq == 9 & (usableDOI$Condition == "ND" | usableDOI$Condition == "NL"), 2, usableDOI$Region)
usableDOI$Region <- if_else(usableDOI$Seq == 10 & (usableDOI$Condition == "ND" | usableDOI$Condition == "NL"), 3, usableDOI$Region)

usableDOI$Region <- if_else(usableDOI$Seq == 7 & (usableDOI$Condition == "ED" | usableDOI$Condition == "EL"), 1, usableDOI$Region)
usableDOI$Region <- if_else(usableDOI$Seq == 8 & (usableDOI$Condition == "ED" | usableDOI$Condition == "EL"), 2, usableDOI$Region)
usableDOI$Region <- if_else(usableDOI$Seq == 9 & (usableDOI$Condition == "ED" | usableDOI$Condition == "EL"), 3, usableDOI$Region)

write.csv(usableDOI,'usableDOI.csv')

#useful to check and make sure everything is properly lined up
usableDOI %>% 
  group_by(Region, Word) %>% 
  summarize(n=n()) -> regionCheck

#create experimental column
#only if the experiment was run in two parts
typeof(usableDOI$Item)
usableDOI$Item <- as.numeric(as.character(usableDOI$Item))
is.numeric(usableDOI$Item)

usableDOI$Experiment <- if_else(usableDOI$Item < 25, 1, 2)


#make factor levels
#remember baseline defaults to first alphabetical factor in each level
#this uses contrast coding
#this also labels factors based on the condition labels, which I always name after the levels in the experimental design
usableDOI$Condition %>%
substr(1,1) %>% factor -> usableDOI$Ellipsis
contrasts(usableDOI$Ellipsis) <- contr.sum(2) * (1/2)
#told it there were two levels, default would be -1 and 1, so multiplying by 1/2 leads to .5 and -.5
#.5 is a good default bc of scaling, don't make it too far from zero
levels(usableDOI$Ellipsis) <- c("Ellipsis", "NoEllipsis")

usableDOI$Condition %>%
substr(2,2) %>% factor -> usableDOI$Locality
contrasts(usableDOI$Locality)  <- contr.sum(2) * (-1/2)
levels(usableDOI$Locality) <- c("Distant", "Local")


#For adding in word frequencies and lengths for the analysis

#read in experiment words that need frequency score
expWordsMaze <- read.csv('criticalWordsMaze.csv', header = 1)

#in case I got word frequencies from the Subtlex corpus
#read in large SUBTLEX file with frequency scores for words in corpus
subtlex <- read.csv('SUBTLEX.csv', header = 1)

expWordsMaze$critical_word <- as.character(expWordsMaze$critical_word)
is.character(expWordsMaze$critical_word)
subtlex$Word <- as.character(subtlex$Word)
is.character(subtlex$Word)

#use dplyr join to merge the list of experimental words and the corpus of all the Subtlex words
#will return only the experimental words, with their frequency
wordFreqMaze <- left_join(expWordsMaze, subtlex, by = c("critical_word" = "Word"))

#find word character length
wordFreqMaze[,"WordLength"] <- NA
wordFreqMaze$WordLength <- str_length(wordFreqMaze$critical_word)

write.csv(wordFreqMaze,'wordFreqMaze_critical.csv')

#add word frequencies and length to main dataframe with dplyr
usableDOI$Word <- (as.character(usableDOI$Word))
finalPlaus <- left_join(usableDOI, wordFreqMaze, by = c("Word" = "critical_word"))


#By region RT analysis

#find Region 2 means
region2 <- subset(finalPlaus, finalPlaus$Region == 2)

#compute RTs means over items
region2 %>% 
  group_by(Locality, Ellipsis, Item) %>% 
  summarize(m=mean(RT), med=median(RT), sd = sd(RT), sd2 = sd(RT)*2, n=n(), se=std.error(RT)) -> summary.by.items_R2

summary.by.items_R2 %>% 
  group_by(Locality, Ellipsis) %>% 
  summarize(RT=mean(m), med=median(med), sd = sd(m), sd2 = sd(m)*2, n=n(), se=sd/sqrt(n)) -> summary.R2

summary.R2$Region <- "NP"



#always check your contrasts again, bc dplyr can mess them up
contrasts(region2$Ellipsis)
contrasts(region2$Locality)

#brms analysis

#always run a fixed effects only model first
#then compare the results with the full  model to ensure no craziness has occurred
brm1 <- brm(RT_log ~ Locality*Ellipsis, 
            data = region2, 
            family = gaussian(link=identity),
            control = list(adapt_delta = 0.8),
            warmup = 1000, iter = 5000, chains = 4) #these are good defaults, can change if warnings occur
summary(brm1, waic = TRUE)

#note: if you increase adapt_delta to .99 and still get divergence warnings, the warnings will send you to STAN
#help pages telling you to re-paramaterize your model. brms already runs this parameter as a default, so your 
#problem is elsewhere. try using a different family to fit your model.

plot_brmR2 <- plot(brmR2)
#will give you all your lovely plots

#for log RTs use family=gaussian; for raw RTs use family=exgaussian
#the exgaussian family already accounts for the gaussian+exponential tail of reaction times
#shifted_log normal also works for raw RTs

#fun full model with random slopes and intercepts
brmR2  <- brm(RT ~ Locality*Ellipsis + WordLength + Zipf.value + #fixed effects
              (0 + Locality*Ellipsis|Participant) + #random slopes for Participant
              (0 + Locality*Ellipsis|Item) +   #random slopes for Items
              (1|Participant) + #random Participant intercept
              (1|Item),   #random Item intercept
             data = region2, 
             family = exgaussian(link=identity),
             control = list(adapt_delta = 0.8),
             warmup = 1000, iter = 5000, chains = 4)
summary(brmR2, waic = TRUE)

plot(brmR2)
#run a posterior check to see how well your model fits the observed data
#shinystan is also available for this, but I've never been able to get it to work well
pp_check(brmR2_mixedLog, nsamples = 300)


#combine all your item regions with dplyr join into a single dataframe for plotting
region_RTs <- full_join(summary.R1, summary.R2)
region_RTs <- full_join(region_RTs, summary.R3)


png("MyMazeRT_plot.png", width=31.75, height=17.46, units="cm", res=300, bg="transparent")
#name your plot, adjust the height, resolution, and make it transparent

#thanks to Steven Foley for sharing his version of this code with me
# Reorder factors so they don't appear alphabetically
region_RTs$Region %<>% factor(levels=c('MC-2', 'MC-1', 'ClauseEnd', 'Quant','NP', 'be', 'Adverb', 'Spill'))

region_RTs$Locality %<>% factor(levels=c('Distant','Local'))

# Make an all-caps legend for the Ellipsis factor
Ellipsis_names <- c('Ellipsis' = "ELLIPSIS", 'NoEllipsis' = "NO ELLIPSIS")

region_plot_NPE3_Maze <- region_RTs %>% 
  ggplot(aes(x=Region, y=RT,
             col=Ellipsis, group=Locality, alpha=Locality))
region_plot_NPE3_Maze + 
  # Facets will be labelled with what we defined
  facet_grid(Ellipsis ~ ., labeller = as_labeller(Ellipsis_names)) +
  # Fiddle with the point/line sizes to make them an appropriate size
  # This gives points at mean RTs and confidence intervals around them
  geom_pointrange(size=1, aes(ymin=RT-se,
                                  ymax=RT+se)) +
  geom_line(lwd=1) + theme_minimal() + 
  scale_shape_circlefill() +
  # More informative/readable region lables
  #scale_x_discrete(labels=c('HdN' = "HdN", 'WhP' = "[ whP", 'Adj' = "Adj", 
                            #'CoArg' = "CoArg", 'XP1' = "XP1", 'XP2' = "XP2", 
                            #'V0' = "Verb ]", 'V1'= "Spill1", 'V2' = "Spill2", 
                            #'V3' = "Spill3")) +
  ylab("Mean RTs (ms)") +
  theme(text=element_text(size=18, family = "Times New Roman"), 
        element_line(size = 2),
        panel.grid = element_line(size=0.5),
        axis.title.x = element_blank()) +
  # Define pretty colors
  scale_color_manual(values=c("#A50D69","#7A0D87","#500DA5")) +
  # Alpha is transparency
  scale_alpha_manual(values=c(1, 0.5),
                     name="Locality",
                     breaks=c("Distant","Local"),
                     labels=c("Distant","Local")) +
  guides(color=FALSE) +
  labs(title ="Cataphoric NPE reading times in Maze task")

dev.off() #will save the file in your R directory



#Failure analysis (where do people fail on the Maze trials)

#this will pull all your trials (using example conditions NL, ND, EL, and ED)
failAll <- subset(allData, allData$Condition == "NL" | allData$Condition == "ND" | allData$Condition == "EL" | allData$Condition == "ED" | | allData$Condition == "F")

#subset all items into only the maze trials
failAll_maze <- subset(failAll, failAll$Controller == "Maze")

#summarize 
failAll_maze %>% 
  group_by(Seq, Correct) %>% 
  summarize(n=n()) -> failAll_mazeSummary

#use spread to reshape data to get average correct by participant on all experimental items
failAll_mazeSummary <- spread(failAll_mazeSummary, Correct, n)
failAll_mazeSummary$total <- failAll_mazeSummary$no + failAll_mazeSummary$yes
failAll_mazeSummary$percCorr <- failAll_mazeSummary$yes / failAll_mazeSummary$total

#where do people fail in the critical items?
failPlaus <- subset(allData, allData$Condition == "NL" | allData$Condition == "ND" | allData$Condition == "EL" | allData$Condition == "ED")

#subset all items into only maze trials and summarize correct
failPlaus_maze <- subset(failPlaus, failPlaus$Controller == "Maze")
failPlaus_maze %>% 
  group_by(Seq, Correct) %>% 
  summarize(n=n()) -> failPlaus_mazeSummary

#look at the total trials with failures
failPlaus_maze_No <- subset(failPlaus_maze, failPlaus_maze$Correct=="no")
failPlaus_maze_No$Seq <- as.numeric(as.character(failPlaus_maze_No$Seq))
is.numeric(failPlaus_maze_No$Seq)

#trim the items because not all of mine were the same length
failPlaus_maze_No_trim <- subset(failPlaus_maze_No, failPlaus_maze_No$Seq < 12)

#summarize for plotting
failPlaus_maze_No_trim %>%
  group_by(Seq, Correct) %>%
  summarize(n=n()) -> failPlausMaze_plot

#Plot the total Maze trials failures per region
png("AllMazeFails.png", width=20.75, height=17.46, units="cm", res=200, bg="transparent")

failPlausMaze_plot %>% ggplot(aes(x=Seq)) -> FailAll

FailAll + geom_bar(stat="identity", aes(y=n), fill="#A50D69")  + theme_minimal(base_size=10, base_family="Times New Roman") + scale_fill_manual(values=c("#A50D69","#7A0D87", "#500DA5")) + ylab("Total Failed Trials") + xlab("Item Region") + 
  theme(legend.text = element_text(size=12, color="black")) + theme(legend.key.height =unit(2,"line")) +
  theme(axis.text.x=element_text(size=12, color ="black")) + theme(axis.text.y=element_text(size=12, color="black")) +
  theme(legend.text = element_text(size=12, color = "black")) + theme(legend.title = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=15, color = "black")) + theme(axis.title.y = element_text(size=15, color = "black"))  +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey60"))  +
  labs(title="Participant total trial failures by region in maze task") + theme(title=element_text(size=15)) + ylim(0, 510) + scale_x_discrete(limits=0:11)

dev.off()

#, x=reorder(Seq, -n)
#this will reorder your x-axis

#Look at only the trials where the failure occurred (non-cumulative fails)
#we can see the trials that people fail on because the button they fail on is correct = 'no' and RT = numeric

failPlaus_maze_fail <- subset(failPlaus_maze_No_trim, failPlaus_maze_No_trim$Correct == "no" & failPlaus_maze_No_trim$RT != "None")
failPlaus_maze_fail %>% 
  group_by(Seq) %>% 
  summarize(n=n()) -> failPlaus_maze_failSummary

failPlaus_maze_failSummary$Seq <- as.numeric(as.character(failPlaus_maze_failSummary$Seq))
is.numeric(failPlaus_maze_failSummary$Seq)

png("NPE3_MazeFail.png", width=20.75, height=17.46, units="cm", res=200, bg="transparent")

failPlaus_maze_failSummary %>% ggplot(aes(x=Seq)) -> Seq

Seq + geom_bar(stat="identity", aes(y=n), fill="#A50D69") + theme_minimal(base_size=10, base_family="Times New Roman") + scale_fill_manual(values=c("#A50D69","#7A0D87", "#500DA5")) + ylab("Total Participant Count") + xlab("Region of Failure") + 
  theme(legend.text = element_text(size=12, color="black")) + theme(legend.key.height =unit(2,"line")) +
  theme(axis.text.x=element_text(size=12, color ="black")) + theme(axis.text.y=element_text(size=12, color="black")) +
  theme(legend.text = element_text(size=12, color = "black")) + theme(legend.title = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=15, color = "black")) + theme(axis.title.y = element_text(size=15, color = "black"))  +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey60"))  +
  labs(title="Participant trial failure by region in maze task") + theme(title=element_text(size=15)) + ylim(0, 63) + scale_x_discrete(limits=0:18)

dev.off()

#Make cumulative distribution plots

maze_plaus$Seq <- as.numeric(as.character(maze_plaus$Seq))
is.numeric(maze_plaus$Seq)
cumdistPlot <- subset(maze_plaus, maze_plaus$Seq < 12)

#repeated code just to get the labels
cumdistPlot$Condition %>%
substr(1,1) %>% factor -> cumdistPlot$Ellipsis
contrasts(cumdistPlot$Ellipsis) <- contr.sum(2) * (1/2)
levels(cumdistPlot$Ellipsis) <- c("Ellipsis", "NoEllipsis")

cumdistPlot$Condition %>%
substr(2,2) %>% factor -> cumdistPlot$Locality
contrasts(cumdistPlot$Locality)  <- contr.sum(2) * (-1/2)
levels(cumdistPlot$Locality) <- c("Distant", "Local")

#this df has all of the failed trials
cumDistPlotFails <- subset(cumdistPlot, cumdistPlot$Correct == "no")

#this df has only the regions on which the failure occurred
cumDistPlotFailedReg <- subset(cumDistPlotFails, cumDistPlotFails$RT != "None")

#this plots the regions on which the failure occurred
png("NPE3_maze_cumDistFailsReg3.png", width=14, height=18, units="cm", res=300, bg="transparent")

cumDistPlotFailedReg %>% ggplot(aes(Correct, col=Ellipsis, linetype=Locality)) -> distFR

distFR + stat_ecdf(aes(x=Seq), geom="point") + scale_x_discrete(limits=0:11) + theme_minimal(base_size=8,base_family="Times New Roman")  + scale_color_manual(values=c("#09406E","#A50D69")) + ylab("Percent Failure") + xlab("Region") + 
  theme(legend.text = element_text(size=11, color="black")) + theme(legend.key.height =unit(1.5,"line")) +
  theme(axis.text.x=element_text(size=11, color ="black")) + theme(axis.text.y= element_text(size=11, color="black")) +
  theme(legend.text = element_text(size=9, color = "black")) + theme(legend.title = element_text(size=9, color = "black")) + 
  theme(axis.title.x = element_text(size=11,color = "black")) + theme(axis.title.y = element_text(size=11,color = "black")) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey88")) + guides(shape = guide_legend(override.aes = list(size=2))) +
  guides(color = guide_legend(override.aes = list(size=.25))) +  labs(title="My maze task", caption="N = #") + theme(title = element_text(size=13)) +
  stat_ecdf(aes(x=Seq), geom="line") + scale_y_continuous(breaks=pretty_breaks(n=10))

dev.off()
#geom="step" gives you straight lines instead of points, take off stat_ecdf(aes(x=Seq), geom="line" for this option
#limits and breaks should be adjusted based on the particular data

#same point with a different visualization, stacked bar plot with condition fill, by region
cumDistPlotFailedReg %>%
  group_by(Condition, Seq, Correct) %>%
  summarize(n=n()) -> MazeFailsBar

png("NPE3_MazeFailBar.png", width=20.75, height=17.46, units="cm", res=200, bg="transparent")

MazeFailsBar %>% ggplot(aes(x=Seq)) -> failBar

failBar + geom_bar(stat="identity", aes(y=n, fill = Condition)) + theme_minimal(base_size=10, base_family="Times New Roman") + scale_fill_manual(values=c("#A50D69","#7A0D87", "#500DA5","#525ecc")) + ylab("Total Failed Trials") + xlab("Item Region") + 
  theme(legend.text = element_text(size=12, color="black")) + theme(legend.key.height =unit(2,"line")) +
  theme(axis.text.x=element_text(size=12, color ="black")) + theme(axis.text.y=element_text(size=12, color="black")) +
  theme(legend.text = element_text(size=12, color = "black")) + theme(legend.title = element_text(size=12, color = "black")) + 
  theme(axis.title.x = element_text(size=15, color = "black")) + theme(axis.title.y = element_text(size=15, color = "black"))  +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color="grey90"))  +
  labs(title="My maze task") + theme(title=element_text(size=15)) + ylim(0, 90) + scale_x_discrete(limits=0:11) 

dev.off()





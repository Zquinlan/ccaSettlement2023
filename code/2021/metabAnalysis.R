# READING -- Libraries ----------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)
library(readxl)
library(randomForest)
library(lme4)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)

#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename
mutate <- dplyr::mutate

zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

angular_transform <- function(x) {
  asin(sqrt(x))
}

genTheme <-  function(x){
  theme(plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
        panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))
}

flagBackgroundNetworks <- function(data, networkColumn, backgroundColumn, minConsensus = 0.5) {
  require("tidyverse")
  
  flagNetworks <- data%>%
    mutate(count = 1)%>%
    group_by({{networkColumn}}, {{backgroundColumn}})%>%
    summarize(count = sum(count))%>%
    ungroup()%>%
    group_by({{networkColumn}})%>%
    mutate(subnetworkPercentReal = count/sum(count),
           backgroundNetworks = case_when(subnetworkPercentReal > 0.5 ~ 'real',
                                          TRUE ~ 'background'))%>%
    filter({{backgroundColumn}} == 'real')%>%
    select({{networkColumn}}, subnetworkPercentReal, backgroundNetworks)
  
  join <- enquo(networkColumn)
  
  flagExport <- data%>%
    left_join(flagNetworks, by = quo_name(join))%>%
    mutate(backgroundNetworks = case_when(is.na(backgroundNetworks) ~ 'background',
                                          TRUE ~ backgroundNetworks))
  
}

# Reading in raw data -----------------------------------------------------
rawMetab <- read_csv('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/metabolomicsQuant.csv')
# rawMetab <- read_csv('~/Documents/SDSU_Scripps/ccaExudatePaper/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-955b2e3b-download_cytoscape_data/quantification_table/quantification_table-00000.csv')

sampleNames <- read_xlsx('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/Orbitrapsequence.xlsx')%>%
  rename(sequenceCode = 1,
         CuracaoCode = 2)%>%
  select(-3)

curacaoNames <- read_csv('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/curacaoNames.csv')%>%
  select(c(1:2, Experiment))%>%
  mutate(CuracaoCode = gsub('Cu', 'Zach_Cu22', CuracaoCode))%>%
  rename(sample = SampleCode)

nodeInfo <- read_tsv('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/nodeInfo.tsv')

allNetworks <- nodeInfo%>%
  select(`cluster index`, componentindex)%>%
  rename(featureNumber = 1,
         network = 2)

conciseRaw <- read_csv('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/conciseConsensus.csv')%>%
  select(-1)%>%
  rename(featureNumber = scan)

libraryMatches <- read_tsv('~/Documents/GitHub/fluffer/data/2021/raw/metabolomics/libraryMatches.tsv')

sirius <- read_tsv('~/Documents/GitHub/fluffer/data/raw/2021/metabolomics/formula_identifications.tsv')%>%
  separate(id, c('num1', 'num2', 'featureNumber'), sep = '_')%>%
  select(molecularFormula, ZodiacScore, featureNumber)

# Cleaning raw data -------------------------------------------------------
cleanMetab <- rawMetab%>%
  select(-c(`row m/z`:`neutral M mass`))%>%
  gather(sequenceCode, xic, 2:ncol(.))%>%
  mutate(sequenceCode = gsub('.mzXML Peak area', '', sequenceCode))%>%
  rename(featureNumber = 1)%>%
  left_join(sampleNames, by = 'sequenceCode')%>%
  left_join(curacaoNames, by = 'CuracaoCode')%>%
  filter(Type %like% '%Blank' | Experiment == 'Fluffer')%>%
  mutate(sample = case_when(is.na(sample) ~ CuracaoCode,
                                  TRUE ~ sample))

networks <- nodeInfo%>%
  select(`cluster index`, componentindex)%>%
  rename(featureNumber = 1,
         network = 2)


# Removing blanks ---------------------------------------------------------
nrow(cleanMetab%>% 
       select(featureNumber)%>%
       unique())

blankRemoval <- cleanMetab%>%
  left_join(networks, by = 'featureNumber')%>%
  group_by(featureNumber)%>%
  mutate(sampleXic = case_when(!Type %like% '%Blank' & !sample %like% 'Blank%' ~ xic,
                               TRUE ~ NA_real_),
         meanSample = mean(sampleXic, na.rm = TRUE),
         #Finding blank and transient features (transient = mean - 1sd)
         blankXic = case_when(Type %like% '%Blank%' ~ xic,
                              TRUE ~ NA_real_),
         maxBlank = max(blankXic, na.rm = TRUE),
         transientNum = sum(sampleXic > 5E4, na.rm = TRUE), 
         # Real columns
         background = case_when(meanSample*0.5 > maxBlank ~ 'real',
                                TRUE ~ 'background'),
         transient = case_when(transientNum >= 3 ~ 'real',
                               TRUE ~ 'transient'))%>%
  ungroup()%>%
  flagBackgroundNetworks(network, background)

nrow(blankRemoval%>% 
       filter(background == 'real')%>%
       select(featureNumber)%>%
       unique())

nrow(blankRemoval%>% 
       filter(background == 'real',
              backgroundNetworks == 'real')%>%
       select(featureNumber)%>%
       unique())

nrow(blankRemoval%>% 
       filter(background == 'real',
              transient == 'real',
              backgroundNetworks == 'real')%>%
       select(featureNumber)%>%
       unique())

postBlanks <- blankRemoval%>%
  filter(background == 'real',
         transient == 'real',
         backgroundNetworks == 'real')%>%
  mutate(log10 = log10(xic +1))


# Cleaning -- min filter ---------------------------------------------
preFilter <- postBlanks%>%
  filter(!Type %like% '%Blank', !sample %like% 'Blank%',
         !sample %like% '%FC%')%>%
  select(featureNumber, network, sample, xic, log10)

minVal <- mean((preFilter%>% filter(log10 > 0))$log10, na.rm = TRUE) 

prevVal <- mean((preFilter%>% filter(log10 > 0))$log10, na.rm = TRUE) - sd((preFilter%>% filter(log10 > 0))$log10, na.rm = TRUE)

pdf('plots/minFilterHistrogram.pdf', height = 12, width = 12)
preFilter%>%
  filter(log10 > 0)%>%
  ggplot(aes(log10)) +
  geom_histogram(bins = 100)+
  # scale_x_log10() +
  geom_vline(xintercept = minVal, color = 'black', alpha = 0.4) +
  genTheme()
dev.off()


postFilter <- preFilter%>%
  group_by(network)%>%
  filter(max(log10) > minVal)%>%
  # mutate(prevalent = sum(mean(log10) > prevVal, na.rm = TRUE))%>%
  # filter(prevalent > 2)%>%
  # select(-prevalent)%>%
  ungroup()

nrow(postFilter%>%
       select(network)%>%
       unique())

nrow(postFilter%>%
       select(featureNumber)%>%
       unique())

# featureList <- postFilter%>%
#   select(featureNumber)%>%
#   unique()
# 
# 
# write_tsv(featureList, 'analysis/featureNumbers.txt')

statsData <- postFilter%>%
  separate(sample, c('experiment', 'water', 'replicate'), sep = '_')

# Pre-Stats -- Random Forest ----------------------------------------------
set.seed(123591)

#NetworBased
randomForestData <- statsData%>%
  ungroup()%>%
  mutate(netFeat = case_when(network == -1 ~ -1*as.numeric(featureNumber),
                             TRUE ~ network))%>%
  group_by(water, replicate, netFeat)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(xic, featureNumber, network))%>%
  spread(netFeat, log10)%>%
  select(-replicate)%>%
  mutate(water = as.factor(water))

#featureBased
# randomForestData <- statsData%>%
#   ungroup()%>%
#   select(-c(xic, experiment))%>%
#   spread(featureNumber, log10)%>%
#   select(-replicate)%>%
#   mutate(water = as.factor(water))


names(randomForestData) <- make.names(names(randomForestData))

rfModel <- randomForest(water ~ ., randomForestData, 
                                importance = TRUE, proximity = TRUE, nPerm = 10000,
                                ntree = 29188)

rfMda <- rfModel$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")


pdf('~/Documents/GitHub/fluffer/data/plots/MDAPlot.pdf', width = 10, height = 10)
rfMda%>%
  # filter(MeanDecreaseAccuracy>0)%>%
  ggplot(aes(x= reorder(feature, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  geom_hline(aes(yintercept = mean(MeanDecreaseAccuracy)+sd(MeanDecreaseAccuracy)))
dev.off()

# ggplot(rfMda, aes(x= reorder(feature, -MeanDecreaseGini), y = MeanDecreaseGini)) +
#   geom_point(stat = "identity") +
#   geom_hline(aes(yintercept = mean(MeanDecreaseGini)+sd(MeanDecreaseGini)))

rfFeatures <- (rfMda%>%
                 mutate(feature = gsub("X", "", feature))%>%
                 filter(MeanDecreaseAccuracy >= mean(MeanDecreaseAccuracy)+sd(MeanDecreaseAccuracy)))$feature%>%
  as.vector()

postRfNetworks <- allNetworks%>%
  mutate(netFeat = case_when(network == -1 ~ -1*as.numeric(featureNumber),
                             TRUE ~ network))%>%
  filter(netFeat %in% rfFeatures)

postRfData <- statsData%>%
  inner_join(postRfNetworks, by = c('featureNumber', 'network'))
  

# Stats -- ttest ----------------------------------------------------------
# lmerTest <- postRfData%>%
#   group_by(network)%>%
#   nest()%>%
#   mutate(mixedModel = map(data, ~ mutate(.x, water = as.factor(water),
#                                        featureNumber = as.factor(featureNumber))%>%
#                           lmer(log10 ~ water + (1|featureNumber) , data =., 
#                                control =  lmerControl(check.nlev.gtr.1 = "ignore",
#                                                       check.conv.singular = 'ignore',
#                                                       check.nobs.vs.nRE = 'ignore'))%>%
#                           car::Anova()%>%
#                           .[['Pr(>Chisq)']]))%>%
#   select(-data)%>%
#   ungroup()

anova <- postRfData%>%
  ungroup()%>%
  mutate(netFeat = case_when(network == -1 ~ -1*as.numeric(featureNumber),
                             TRUE ~ network))%>%
  group_by(water, replicate, netFeat)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(netFeat)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(log10 ~ water, data = .x)%>%
                       tidy()%>%
                       select(p.value)))%>%
  unnest(data)%>%
  mutate(fdr = p.adjust(p.value, method = 'BH'))%>%
  filter(fdr < 0.05)


sigFeatures <- anova$netFeat%>%
  unique()

# PLOTTING -- Significant Features -------------------------------------------------------------
ccaFeatures <- postRfData%>%
  select(-xic)%>%
  mutate(netFeat = case_when(network == -1 ~ -1*as.numeric(featureNumber),
                             TRUE ~ network))%>%
  spread(water, log10)%>%
  group_by(netFeat, replicate)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(netFeat)%>%
  filter(mean(PPL) > mean(FSW, na.rm = TRUE))
  

ccaFeatureList <- ccaFeatures%>%
  select(netFeat)%>%
  unique()


sigFeatureData <- postRfData%>%
  mutate(netFeat = case_when(network == -1 ~ -1*as.numeric(featureNumber),
                             TRUE ~ network))%>%
  filter(netFeat %in% sigFeatures)%>%
  inner_join(ccaFeatureList, by = c('netFeat'))%>%
  left_join(conciseRaw, by = c('featureNumber', 'network'))%>%
  mutate(conciseConsensus = case_when(is.na(conciseConsensus) ~ 'Unknown',
                                      TRUE ~ conciseConsensus))%>%
  unite(label, c('conciseConsensus', 'network'), sep = '  ', remove = FALSE)%>%
  unite(sample, c('water', 'replicate'), sep = '_', remove = FALSE)

sigData <- sigFeatureData%>%
  group_by(label, netFeat, network, sample, water, replicate, conciseConsensus)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()

labelLevels <- sigData$label%>% 
  as.factor()

conciseLevels <-sigData$conciseConsensus%>%
  as.factor()

pdf('plots/pplMetabolites.pdf', width = 20, height = 30)
sigData%>%
  group_by(label)%>%
  mutate(zscore = zscore(log10))%>%
  ggplot(aes(sample, label, fill = zscore)) +
  geom_tile() +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  scale_y_discrete(limits = rev(levels(labelLevels))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))
dev.off()

pdf('plots/pplMetaboliteBars.pdf', width = 20, height = 45)
sigData%>%
  ggplot(aes(conciseConsensus, log10, fill = water, color = water)) +
  geom_violin(alpha = 0.2, position = 'nudge') +
  geom_point(size = 3) +
  scale_color_manual(values = c('#78B7C5', "#E1AF00")) + 
  scale_fill_manual(values = c('#78B7C5', "#E1AF00")) +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(conciseLevels))) +
  genTheme()+ 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 8))
dev.off()




majorFeatures <- (sigData%>%
  filter(log10 >  mean(log10) + sd(log10)))$featureNumber

selectSigData <- sigData%>%
  filter(featureNumber %in% majorFeatures)

selectLabelLevels <- selectSigData$label%>% 
  as.factor()

selectConciseLevels <-selectSigData$conciseConsensus%>%
  as.factor()

selectSigData%>%
  group_by(label)%>%
  mutate(zscore = zscore(log10))%>%
  ggplot(aes(sample, label, fill = log10)) +
  geom_tile() +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  scale_y_discrete(limits = rev(levels(selectLabelLevels))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15))


# Chemical superclasses ---------------------------------------------------
superclasses <- sigData%>%
  separate(conciseConsensus, c('superclass', 'class', 'subclass'), sep = ';')
  # group_by(sample, water, superclass)%>%
  # summarize_if(is.numeric, sum)%>%
  # ungroup()

superclassLevels <- superclasses$superclass%>% 
  as.factor()

pdf('plots/pplSuperclassViolins.pdf', width = 15, height = 15)
superclasses%>%
  ggplot(aes(superclass, log10, fill = water, color = water)) +
  geom_violin(alpha = 0.2, position = 'nudge') +
  geom_point(size = 3) +
  scale_color_manual(values = c('#78B7C5', "#E1AF00")) + 
  scale_fill_manual(values = c('#78B7C5', "#E1AF00")) +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(superclassLevels))) +
  genTheme()+ 
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 8))
dev.off()

superclassPercents <- superclasses%>%
  group_by(superclass, sample)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(sample)%>%
  mutate(percent = log10/sum(log10))%>%
  ungroup()%>%
  filter(sample %like% 'PPL%')

superclassPercents%>%
  ggplot(aes(percent, superclass)) +
  geom_point(size = 4) +
  gen_theme()

superclassAverages <- superclassPercents%>%
  group_by(superclass)%>%
  mutate(sd = sd(percent))%>%
  summarize_if(is.numeric, mean)


# cytoscape ---------------------------------------------------------------

cytoData <- sigFeatureData%>%
  separate(conciseConsensus, c('superclass', 'class', 'subclass'), sep = ';')%>%
  select(featureNumber, network, water, replicate, xic, superclass, class, subclass)%>%
  group_by(featureNumber, network, water, superclass, class, subclass)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  left_join(libraryMatches%>%
              select('#Scan#', 'Compound_Name'), by = c('featureNumber' = '#Scan#'))%>%
  spread(water, xic)%>%
  mutate(size = log10(PPL))

write_csv(cytoData, '~/Documents/GitHub/fluffer/data/analysis/cytoscapeSig.csv')


# Significant points ------------------------------------------------------
ignoreCCA <- ccaFeatureList$netFeat

insignificantPoints <- postRfData%>%
  filter(!netFeat %in% ignoreCCA)%>%
  select(-xic)%>%
  group_by(featureNumber, network, water)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(water, log10)%>%
  mutate(superclass = 'Insignificant')

pdf('~/Documents/GitHub/fluffer/data/2021/plots/PPLdominance.pdf', width = 15, height = 10)
cytoData%>%
  filter(network != 351)%>%
  mutate(FSW10 = log10(FSW +1),
         PPL10 = log10(PPL + 1))%>%
  ggplot() +
  geom_point(data = insignificantPoints, aes(FSW, PPL), color = 'grey', size = 4, alpha = 0.6, shape = 15) +
  geom_point(aes(FSW10, PPL10, color = superclass), size = 4) +
  genTheme() +
  xlim(0,8) +
  ylim(0,8) +
  labs(x = bquote(Water ~control ~production ~log[10](XIC)), y = bquote(Crustose ~coralline ~algae ~production ~log[10](XIC)), color = 'Superclass Concensus')
dev.off()

pplDominant <- c('3226', '413', '3891', '409', '6479')

pdf('~/Documents/GitHub/fluffer/data/2021/plots/PPLNetworkProduction.pdf', width = 15, height = 10)
sigFeatureData%>%
  filter(network %in% pplDominant)%>%
  separate(conciseConsensus, c('superclass', 'class', 'subclass'), sep = ';')%>%
  group_by(network, water, replicate, superclass, class, subclass)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(network = as.character(network),
         log10 = log10(xic))%>%
  ggplot(aes(network, log10, color = water)) +
  geom_boxplot(size = 2) +
  genTheme() +
  scale_color_manual(values = c('#78B7C5', "#E1AF00")) 
dev.off()


# EXPORT -- SIgnificant export data ---------------------------------------
exportSigData <- sigFeatureData%>%
  group_by(featureNumber, network, water, netFeat, conciseConsensus, conciseConsensusLevel, conciseConsensusScore, matchSource)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(log10, numberOfNodes, netFeat))%>%
  group_by(network, water)%>%
  mutate(networkXIC = sum(xic))%>%
  select(-xic)%>%
  spread(water, networkXIC)%>%
  ungroup()%>%
  left_join(sirius%>%
              mutate(featureNumber = as.numeric(featureNumber))%>%
              filter(ZodiacScore > 0.5), by = 'featureNumber')%>%
  left_join(libraryMatches%>%
              select(Compound_Name, Instrument, Smiles, INCHI, `#Scan#`)%>%
              rename(libraryMatchID = Compound_Name,
                     featureNumber = `#Scan#`,
                     libraryInstrument = Instrument, 
                     librarySmiles = Smiles, 
                     libraryINCHI = INCHI),
            by = 'featureNumber')%>%
  left_join(rawMetab%>%
              select(1:3)%>%
              rename(featureNumber = 1,
                     `m/z` = `row m/z`,
                     rententionTime = 3), by = 'featureNumber')%>%
  # select(-featureNumber)%>%
  rename(CCANetworkXIC = PPL,
         WaterNetworkXIC = FSW,
         subnetwork = network)%>%
  filter(subnetwork != 351,
         matchSource != 'Library' | !is.na(libraryMatchID))%>%
  unique()%>%
  select(CCANetworkXIC, WaterNetworkXIC, featureNumber, subnetwork, `m/z`, rententionTime, everything())


write_csv(exportSigData, '~/Documents/GitHub/fluffer/data/2021/analysis/ccaExometaboliteComposition.csv')










##

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

gen_theme <-  function(x){
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


# READING -- raw dataframes -----------------------------------------------
## surface area measurements
surfaceAreasRaw <- read_csv('~/Documents/GitHub/fluffer/data/2021/raw/surfaceArea.csv')%>%
  rename('W' = 3)

# settlementCount1 <- read_csv('~/Documents/GitHub/fluffer/data/raw/SettlementCount1.csv')%>%
#   mutate(prcntSettlement = `Attached/Settled`/Total,
#          Concentration = as.factor(Concentration))

setwd('~/Documents/GitHub/fluffer/data/2021/raw/settlement')
settlement <- dir(path = "~/Documents/GitHub/fluffer/data/2021/raw/settlement", pattern = ".csv")%>%
  map(read_csv)%>%
  reduce(bind_rows)%>%
  mutate(prcntSettlement = `Attached/Settled`/Total,
         Concentration = as.character(Concentration),
         Concentration = case_when(Treatment %like% '%control' ~ 'control',
                                          TRUE ~ as.character(Concentration)))
setwd('~/Documents/GitHub/fluffer/data/2021/')

# ANALYSIS -- surface area averages ---------------------------------------
surfaceAreas <- surfaceAreasRaw%>%
  mutate(L = L/10,
         W = W/10,
         H = H/10,
         surfaceArea = case_when(Shape == 'Cylinder' ~ (2*pi*L*(W/2)),
                                 Shape == 'Rectangle; -bottom and end caps' ~ (L*W + L*H*2),
                                 Shape == 'Rectangle' ~ (2*(W*L + H*L + H*W)),
                                 Shape == 'Triangle' ~ (((L/2)*H)*2 + 2*sqrt((L/2)^2 + H^2) + L*W)))

averageSurfaces <- surfaceAreas%>%
  group_by(Replicate)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

mean(averageSurfaces$surfaceArea)
sd(averageSurfaces$surfaceArea)

c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")


# ANALYSIS -- Average settlement ------------------------------------------
averageControl <- settlement%>%
  filter(Concentration == 'control',
         Date == '11022021')%>%
  mutate(std = sd(prcntSettlement))%>%
  summarize_if(is.numeric, mean)

paste(averageControl$prcntSettlement, '±', averageControl$std)

averageFull <- settlement%>%
  filter(Date == '11022021')%>%
  group_by(Treatment, Concentration)%>%
  mutate(n = 1,
         n = sum(n),
         std = sd(prcntSettlement)/sqrt(n))%>%
  summarize_if(is.numeric, mean)



# qqplots -----------------------------------------------------------------
settlementStats <- settlement%>%
  filter(Tray != '1' & Well != 'A3')%>%
  mutate(asin = asin(sqrt(prcntSettlement)),
         percent = prcntSettlement*100)
  
pdf('plots/settlementDistributionQqplots.pdf', width = 10, height = 10)
car::qqPlot(settlementStats$prcntSettlement,
            xlab = 'Normal Quantiles', ylab = 'Percent Settlement Quantiles')

car::qqPlot(settlementStats$asin,
            xlab = 'Normal Quantiles', ylab = 'Angulart Transformed Percent Settlement Quantiles')
dev.off()

# STATS -- ANOVA -------------------------------------------------------------------
resinVector <- as.vector(c("PPL", "HLB", "C18"))

#having water control be at all concentrations
waterControl10 <- settlementStats%>%
  filter(Treatment == 'Water control')%>%
  mutate(Concentration = case_when(Concentration == 'control' ~ '10',
                                   TRUE ~ Concentration))

waterControl1000 <- settlementStats%>%
  filter(Treatment == 'Water control')%>%
  mutate(Concentration = case_when(Concentration == 'control' ~ '1000',
                                   TRUE ~ Concentration))

# Dunnetts test with binding the control treatments 
WaterControlDunn <- settlementStats%>%
  bind_rows(waterControl10, waterControl1000)%>%
  filter(Concentration != 'control' | Concentration == 'control' & Treatment == 'Water control',
         !Treatment %like% 'Fire%',
         Date == '11022021')%>%
  mutate(Concentration = case_when(Concentration == 'control' ~ '100',
                                   TRUE ~ Concentration),
         Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, 'Water control', 'C18', 'HLB', 'PPL', 'CCA exudate'))%>%
  group_by(Concentration)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Treatment, data = .x)%>%
                      glht(linfct = mcp(Treatment = "Dunnett"))%>%
                      summary()%>%
                      tidy()))%>%
  unnest(data)%>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

# Fractionation dunnetts
fractionationDunn <- settlementStats%>%
  filter(!Treatment %in% c('Water control', 'CCA exudate'),
         !Treatment %like% 'Fire%',
         Date == '11022021')%>%
  mutate(Treatment = gsub(' control', '', Treatment),
         Concentration = as.factor(Concentration),
         Concentration = fct_relevel(Concentration, 'control', '10', '100', '1000'))%>%
  group_by(Treatment)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ Concentration, data = .x)%>%
                      glht(linfct = mcp(Concentration = "Dunnett"))%>%
                      summary()%>%
                      tidy()%>%
                      # select(-c(4:7))%>%
                      mutate(lhs = gsub(" - control", " - resin control", lhs))))%>%
  unnest(data)%>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

  
  

# VIZUALIZATIONS -- settlement Fractions--------------------------------------------
pdf('~/Documents/GitHub/fluffer/data/2021/plots/fractionationSettlement.pdf', width = 15, height = 10)
settlementStats%>%
  mutate(Treatment = case_when(Treatment == 'CCA exudate' ~ 'Non-fractionated exometabolites',
                               TRUE ~ as.character(Treatment)),
         Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, 'Water control', 'C18 control', 'HLB control', 'PPL control', 'Non-fractionated exometabolites'),
         Concentration = case_when(Concentration == '10' ~ '10%',
                                   Concentration == '100' ~ '100%',
                                   Concentration == '1000' ~ '1000%',
                                   TRUE ~ as.character(Concentration)))%>% 
  filter(!Treatment %like% 'Fire%',
         Date == '10312021')%>%
  ggplot(aes(Treatment, percent, color = Concentration)) +
  geom_boxplot() +
  # geom_point(stat = 'identity', position = position_dodge2(width = 1)) +
  scale_color_manual(values = c("#669900", "#FF850A", "#A30029", "#3B9AB2")) +
  gen_theme() +
  theme(legend.position = 'top') +
  labs(y = 'Percent Settlement', x = 'Exometabolite treatment')

settlementStats%>%
  mutate(Treatment = case_when(Treatment == 'CCA exudate' ~ 'Non-fractionated exometabolites',
                               TRUE ~ as.character(Treatment)),
         Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, 'Water control', 'C18 control', 'HLB control', 'PPL control', 'Non-fractionated exometabolites'),
         Concentration = case_when(Concentration == '10' ~ '10%',
                                   Concentration == '100' ~ '100%',
                                   Concentration == '1000' ~ '1000%',
                                   TRUE ~ as.character(Concentration)))%>% 
  filter(!Treatment %like% 'Fire%',
         Date == '11012021')%>%
  ggplot(aes(Treatment, percent, color = Concentration)) +
  geom_boxplot() +
  # geom_point(stat = 'identity', position = position_dodge2(width = 1)) +
  scale_color_manual(values = c("#669900", "#FF850A", "#A30029", "#3B9AB2")) +
  gen_theme() +
  theme(legend.position = 'top') +
  labs(y = 'Percent Settlement', x = 'Exudate treatment')

settlementStats%>%
  mutate(Treatment = case_when(Treatment == 'CCA exudate' ~ 'Non-fractionated exometabolites',
                               TRUE ~ as.character(Treatment)),
         Treatment = as.factor(Treatment),
         Treatment = fct_relevel(Treatment, 'Water control', 'C18 control', 'HLB control', 'PPL control', 'Non-fractionated exometabolites'),
         Concentration = case_when(Concentration == '10' ~ '10%',
                                   Concentration == '100' ~ '100%',
                                   Concentration == '1000' ~ '1000%',
                                   TRUE ~ as.character(Concentration)))%>% 
  filter(!Treatment %like% 'Fire%',
         Date == '11022021')%>%
  ggplot(aes(Treatment, percent, color = Concentration)) +
  geom_boxplot() +
  # geom_point(stat = 'identity', position = position_dodge2(width = 1)) +
  scale_color_manual(values = c("#669900", "#FF850A", "#A30029", "#3B9AB2")) +
  gen_theme() +
  theme(legend.position = 'top') +
  labs(y = 'Percent Settlement', x = 'Exudate treatment')
dev.off()


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


# LOADING -- Dataframes ---------------------------------------------------
raw_scoring_df <- read_csv("~/Documents/GitHub/fluffer/data/2019//raw/CRACK_Scoring.csv")%>%
  filter(!Sample_code == "PP_TI2_AP",
         !Sample_code == "PS_TI4_MF",
         !Sample_code == "PP_TI2_MF",
         !Sample_code == "HB_TI4_MF")

raw_scoring_df[is.na(raw_scoring_df)] <- 0


# CLEANING -- Scoring dataframe to only settled vs not settled-------------------------------------------
tidy_scoring <- raw_scoring_df%>%
  select(-c(2:7))%>%
  mutate(percent = (Sum_settled/Alive)*100,
         percent.asin = asin(sqrt((Sum_settled/Alive))))%>%
  separate(Sample_code, c("CCA_species", "metabolite_pool", "Coral_species"), sep = "_")%>%
  separate(metabolite_pool, c('metabolite_pool', 'replicate'), sep = 2)%>%
  mutate(CCA_species = case_when(CCA_species == "HB" ~ "Hydrolithon borgesenii",
                               CCA_species == "PS" ~ "Paragoniolithon solubile",
                               CCA_species == "PP" ~ "Porolithon pachydernum",
                               CCA_species == "WA" ~ "Water control",
                               TRUE ~ as.character(CCA_species)))%>%
  mutate(metabolite_pool = case_when(metabolite_pool == "EX" ~ "Exometabolites",
                                     metabolite_pool == "TI" ~ "Tissue metabolites",
                                     TRUE ~ as.character(metabolite_pool)))%>%
  mutate(Coral_species = case_when(Coral_species == "AP" ~ "Acropora palmata",
                                   Coral_species == "DL" ~ "Diploria labyrinthiformis",
                                   Coral_species == "MF" ~ "Montastrea faveolata",
                                   TRUE ~ as.character(Coral_species)))

tidy_scoring$CCA_species <- as.factor(relevel(factor(tidy_scoring$CCA_species), "Water control"))

averageSettled <- tidy_scoring%>%
  group_by(CCA_species, metabolite_pool, Coral_species)%>%
  mutate(n = 1,
         n = sum(n),
         sterr = sd(percent)/sqrt(n))%>%
  summarize_if(is.numeric, mean)

# qqPlots -----------------------------------------------------------------

pdf('plots/settlementDistribution.pdf', width = 10, height = 10)
car::qqPlot(tidy_scoring$percent,
            xlab = 'Normal Quantiles', ylab = 'Percent Settlement Quantiles')
car::qqPlot(tidy_scoring$percent.asin,
            xlab = 'Normal Quantiles', ylab = 'Angulart Transformed Percent Settlement Quantiles')
dev.off()


# PRE-STATS -- Set seed ---------------------------------------------------
set.seed(29580)

# STATS -- TWO-Way ANOVA Scoring --------------------------------------------------
twoway_scoring_aov <- tidy_scoring%>%
  group_by(Coral_species)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(percent.asin~CCA_species*metabolite_pool, data = .x)%>%
                       tidy()))%>%
  unnest(data)


# STATS -- Dunnetts scoring -----------------------------------------------
cca_species_vector <- as.vector(c("Hydrolithon borgesenii", "Paragoniolithon solubile", "Porolithon pachydernum"))

dunnett_model_scoring <- tidy_scoring%>%
  filter(CCA_species != 'Porolithon pachydernum')%>%
  group_by(Coral_species, metabolite_pool)%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(percent.asin ~ CCA_species, .x)%>%
                    glht(linfct = mcp(CCA_species = "Dunnett"))),
      dunnett_summary = map(dunnett, ~ summary(.x)%>%
                              tidy()))%>%
  dplyr::select(-c(data, dunnett))%>%
  unnest(dunnett_summary)%>%
  dplyr::select(-c(4:7))%>%
  mutate(lhs = gsub(" - Water control", "", lhs),
         FDR = p.adjust(p.value, method = 'BH'))

# GRAPHING -- Scoring Bar Charts ------------------------------------------
no_pachy <- tidy_scoring%>%
  filter(!CCA_species == "Porolithon pachydernum")%>%
  mutate(CCA_species = case_when(CCA_species %like% 'Hydrolithon%' ~ 'H. boergesenii',
                                 CCA_species %like% 'Paragon%' ~ 'P. solubile',
                                 TRUE ~ as.character(CCA_species)),
         CCA_species = as.factor(CCA_species),
         CCA_species = fct_relevel(CCA_species, 'Water control', 'H. boergesenii', 'P. solubile'))

pdf("~/Documents/GitHub/fluffer/data/2019/plots/perent_settlement.pdf",width = 15, height = 13)
no_pachy%>%
  ggplot(aes(x = CCA_species, y = percent))+
  geom_boxplot(aes(fill= metabolite_pool),
               stat = "boxplot", position = "dodge2") +
  geom_point(position=position_dodge(width=0.75), aes(group=metabolite_pool))+
  scale_fill_manual(values = c('#EBCC2A', '#006658')) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    # panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',colour = "black"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.text.x = element_text(angle = 75, hjust = 1,face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    strip.text = element_text(face = "italic")
  ) +
  facet_wrap(~ reorder(Coral_species, percent), nrow = 1) +
  xlab("CCA species") +
  ylab("Percent Settled")

no_pachy%>%
  filter(metabolite_pool == 'Tissue metabolites')%>%
  ggplot(aes(x = CCA_species, y = percent))+
  geom_boxplot(aes(fill = 'Grey'),
               stat = "boxplot", position = "dodge2") +
  geom_point(position=position_dodge(width=0.75), aes(group=metabolite_pool)) +
  scale_fill_manual(values = 'Grey') +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    # panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',colour = "black"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.text.x = element_text(angle = 75, hjust = 1, size = 30, face = "italic"),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = 'None'
  ) +
  facet_wrap(~ reorder(Coral_species, percent), nrow = 1) +
  xlab("CCA species") +
  ylab("Percent Settled")

no_pachy%>%
  mutate(order = case_when(Coral_species == 'Montastrea faveolata' ~ 2,
                           Coral_species == 'Acropora palmata' ~ 1,
                           TRUE ~ 3))%>%
  filter(metabolite_pool == 'Exometabolites')%>%
  ggplot(aes(x = CCA_species, y = percent))+
  geom_boxplot(aes(fill = 'Grey'),
               stat = "boxplot", position = "dodge2") +
  geom_point(position=position_dodge(width=0.75), aes(group=metabolite_pool))+
  scale_fill_manual(values = 'Grey') +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    # panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',colour = "black"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.text.x = element_text(angle = 75, hjust = 1, size = 30, face = "italic"),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 25),
    strip.text = element_text(face = "italic"),
    legend.position = 'None'
  ) +
  facet_wrap(~ reorder(Coral_species, order), nrow = 1) +
  xlab("CCA species") +
  ylab("Percent Settled")

dev.off()


# WRITING -- stats_dataframes ---------------------------------------------
write_csv(dunnett_model_scoring, "~/Documents/GitHub/fluffer/data/2019/analysis/settlement_dunnet_pvalues.csv")

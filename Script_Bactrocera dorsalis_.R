#*******************#
#### NGB - FRUIT ####
#*******************#

#Packages----
library(tidyverse)
library(phyloseq)
library(ggrepel)
library(vegan)
library(microDecon)
library(paletteer)

#Load data----
#**Nouvelle taxo----
silva_data = "databioinfo/tax_slv_ssu_138.1.txt.gz"
illumina_data = "dataFRUIT/NGB_16S_Illumina_214_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table3"

#Digestion de la base de références SILVA----
#**Functions----
load_tax_slv_data <- function(tax_slv_data) {
  column_names <- c("taxonomic_path", "taxid", "taxonomic_rank", "remark", "release")
  columns_to_keep <- c("taxonomic_path", "taxonomic_rank")
  read_tsv(tax_slv_data,
           col_names = column_names,
           col_types = "cifcf",
           col_select = all_of(columns_to_keep),
           show_col_types = FALSE)
}

. %>%
  mutate(taxonomic_path = str_remove(taxonomic_path, ";$"),
         taxonomic_path = str_remove(taxonomic_path, ".*;")) %>%
  arrange(taxonomic_path) %>%
  count(taxonomic_path, taxonomic_rank) -> link_taxa_and_ranks

. %>%
  count(taxonomic_path) %>%
  filter(n > 1) -> find_duplicated_ranks


. %>%
  filter(! taxonomic_path %in% c("Incertae Sedis", "uncultured") &
           ! (taxonomic_path == "SAR" & taxonomic_rank == "phylum") &
           ! (taxonomic_path == "Stramenopiles" & taxonomic_rank == "subphylum") &
           ! (taxonomic_path == "Labyrinthulomycetes" & taxonomic_rank == "class")) %>%
  select(-n) -> ban_some_taxa_and_ranks

load_assignments <- function(assignment_data){
  columns_to_keep <- c("OTU", "taxonomy")
  read_tsv(assignment_data,
           col_select = all_of(columns_to_keep),
           show_col_types = FALSE)
}

#**Main----
# find potential issues
load_tax_slv_data(silva_data) %>%
  link_taxa_and_ranks %>% 
  find_duplicated_ranks

##   taxonomic_path          n
##   <chr>               <int>
## 1 Incertae Sedis         10
## 2 Labyrinthulomycetes     2
## 3 SAR                     2
## 4 Stramenopiles           2
## 5 uncultured              8

## Besides 'uncultured' and 'Incertae Sedis',
## all conflicts come from these four entries:
## Eukaryota;Amorphea;Amoebozoa;SAR;	46930	phylum		138
## Eukaryota;Amorphea;Amoebozoa;SAR;Stramenopiles;	46931	subphylum		138
## Eukaryota;Amorphea;Amoebozoa;SAR;Stramenopiles;Labyrinthulomycetes;	46932	class		138
## Eukaryota;Amorphea;Amoebozoa;SAR;Stramenopiles;Labyrinthulomycetes;Sorodiplophrys;	46933	genus		138

## handle these specific cases
load_tax_slv_data(silva_data) %>%
  link_taxa_and_ranks %>%
  ban_some_taxa_and_ranks -> taxa2ranks

taxa2ranks %>%
  head(n = 20)

## how many taxonomic levels in Silva?
taxa2ranks %>%
  select(taxonomic_rank) %>%
  distinct()  # 20 levels


## annotate our taxonomic paths
good_taxonomic_levels <- c("domain", "phylum", "class",
                           "order", "family", "genus")
#*******************----
#Metadata----
"dataFRUIT/coordFRUIT_corr_eco.csv" %>% 
  read_delim(show_col_types = FALSE, delim = ";") -> df_lieux

. %>% 
  mutate(temp = row_number(), 
       temp = str_pad(temp, width=2, side="left", pad = "0")) %>% 
  unite(Label, temp, Label, sep = "_") -> number_Label

"dataFRUIT/metadata_illumina.csv" %>% 
  read_csv2(show_col_types = FALSE) %>% 
  filter(experiment == "RUN4") %>% 
  rename(locality = locality...9,
         population = locality...10) %>% 
  arrange(env, locality) %>%
  number_Label %>% 
  left_join(., df_lieux, by = c("env","locality")) %>% 
  mutate(row_label = Label,
         grp_plu = as.character(grp_plu)) %>%
  column_to_rownames(. ,var = "row_label") -> my_mtd

#Estimation de la distribution de nos variables à travers nos échantillons 
ggplot(df_lieux, aes(x=join_Tmoy, y =join_pluie_moy, col = env)) +
  geom_point() +
  geom_text_repel(aes(label = env, fontface = 'italic'), size =
                    3, max.overlaps = 100) +
  geom_segment(aes(x = 20, y = 0, xend = 20, yend = 4000), color = "red", linetype = "dashed", linewidth = 1) +
  geom_segment(aes(x = 22, y = 0, xend = 22, yend = 4000), color = "red", linetype = "dashed", linewidth = 1) +
  geom_segment(aes(x = 18, y = 1000, xend = 24, yend = 1000), color = "red", linetype = "dashed", linewidth = 1) + 
  geom_segment(aes(x = 18, y = 2000, xend = 24, yend = 2000), color = "red", linetype = "dashed", linewidth = 1)

#Basic statistics----
my_mtd %>%
  select(sample_ID) %>%
  pull(sample_ID) -> all_RUN4_samples

my_mtd %>%
  filter(type=="ENV") %>% 
  select(sample_ID) %>%
  pull(sample_ID) -> only_environmental_samples 

#**Nombre de séquences----
columns_to_keep <- c("OTU","references",only_environmental_samples)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% View()
  select(all_of(columns_to_keep)) %>%
  mutate(.before = 1, total = rowSums(.[,3:ncol(.)],na.rm=T)) %>% 
  filter(total > 0) %>% 
  select(-total) -> illumina_data_RUN4 

illumina_data_RUN4 %>% dim()
#Dans les échantillons environnementals, il y a 559 OTU

#**Nombre de séquences non assignées---- 
# illumina_data_RUN4 %>%
#   filter(is.na(references)) %>% view()
#Dans le RUN4, il y a 1 seule OTU associée à aucune référence. Elle a très peu de reads.

#*******************----
#Building a phyloseq object based on OTUs----
#**OTU----
columns_to_keep <- c("OTU", "references", all_RUN4_samples)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  filter(!is.na(references)) %>% 
  select(all_of(columns_to_keep)) %>% 
  replace(is.na(.),0) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>% 
  select(-total,-references) %>% 
  column_to_rownames(var = "OTU") -> my_otu

colnames(my_otu) <- my_mtd$Label

#**Tax à l'ancêtre commun----
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(OTU,taxonomy) %>% 
  filter(!is.na(taxonomy)) %>% 
  filter(OTU %in% rownames(my_otu)) %>% 
  separate_rows(taxonomy, sep = "\\|") %>% 
  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy" = "taxonomic_path")) %>% 
  filter(taxonomic_rank %in% good_taxonomic_levels) %>% #enlève les mauvais niveaux taxo dont les NA
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>% 
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy) %>%  
  column_to_rownames(var = "OTU") %>% 
  cbind(.,OTU = rownames(my_otu)) -> my_tax

#*******************----
#MicroDecon----
# install.packages("remotes")
# remotes::install_github("donaldtmcknight/microDecon")
library("microDecon")

my_mtd %>% 
  filter(type=="NEG") %>% 
  pull(Label) -> neg_samples

my_mtd %>% 
  filter(type=="ENV") %>% 
  pull(Label) -> env_samples

neg <- my_otu[,neg_samples]
env <- my_otu[, env_samples]

data.frame(OTU = rownames(my_otu)) %>% 
  cbind(., neg) %>% 
  cbind.data.frame(., env) -> df_decon

decontaminated <- decon(data = df_decon,numb.blanks=length(colnames(neg)),numb.ind=length(colnames(env)), taxa=F) 

#**OTUs à retirer----
decontaminated$OTUs.removed %>% 
  pull(OTU)-> otus_to_remove
length(otus_to_remove) #il y en a 95

#**OTUs seulement présentes dans les puits négatifs----
my_otu %>% 
  mutate(.before = 1, total_env = rowSums(.[,env_samples])) %>%
  filter(total_env==0) %>% 
  rownames(.)-> OTU_absent_from_env_samples

length(OTU_absent_from_env_samples) #Il y en a 81

#**OTUs dont les reads sont corrigés----
decontaminated$reads.removed %>%
  replace(.,is.na(.),0) %>%
  filter(!OTU %in% otus_to_remove) %>% 
  filter(!OTU %in% OTU_absent_from_env_samples) %>% 
  pull(OTU) -> partially_removed 
length(partially_removed) #Il y en a 19

#**Phyloseq object pour faire graph OTU conta----

otu_table(my_otu,taxa_are_rows = TRUE) -> OTU.decon 
tax_table(as.matrix(my_tax)) -> TAX.decon
sample_data(my_mtd) -> MTD

ps_decon = phyloseq(OTU.decon, TAX.decon, MTD)

#**Visualisation of decontamination----
ps.pa <- transform_sample_counts(ps_decon, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "NEG", ps.pa) 
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$type == "NEG"), ps.pa)

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=case_when(
                      taxa_names(ps.pa) %in% OTU_absent_from_env_samples ~ "Only neg",
                      taxa_names(ps.pa) %in% otus_to_remove ~ "Full contaminant",
                      taxa_names(ps.pa) %in% partially_removed ~ "Partial contaminant", 
                      TRUE ~ "True OTU"
                    ))

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)")  + 
  ggtitle("Prevalence of conserved and removed OTUs following microDecon")

#Graphique avec les nombres de reads 
ps.decon.neg <- prune_samples(sample_data(ps_decon)$type == "NEG", ps_decon) 
ps.decon.pos <- prune_samples(!(sample_data(ps_decon)$type == "NEG"), ps_decon)

df.pa1 <- data.frame(
  OTU = taxa_names(ps.pa),
  pa.pos1 = log10(taxa_sums(ps.decon.pos) + 1),
  pa.neg1 = log10(taxa_sums(ps.decon.neg) + 1),
  contaminant = case_when(
    taxa_names(ps.pa) %in% OTU_absent_from_env_samples ~ "Only neg",
    taxa_names(ps.pa) %in% otus_to_remove ~ "Full contaminant",
    taxa_names(ps.pa) %in% partially_removed ~ "Partial contaminant",
    TRUE ~ "True OTU"
  )
)

ggplot(data=df.pa1, aes(x=pa.neg1, y=pa.pos1, color=contaminant)) + 
  geom_point(position = "jitter")+
  geom_label_repel(aes(label = OTU,  color = contaminant), size = 3)+
  xlab("Number of reads in negative Controls (Log scale))") + 
  ylab("Number of reads in true Samples (Log-scale)")  + 
  ggtitle("Decontamination following microDecon")

my_tax %>% 
  filter(OTU == "5")

#*******************----
#Objets finaux----
#**Table d'OTU décontaminée----
decontaminated$decon.table %>% 
  mutate(total = rowSums(.[,3:ncol(.)])) %>%
  filter(total > 0) %>%
  select(-total,-Mean.blank) %>% 
  column_to_rownames(var = "OTU") -> my_cleaned_otu_table
#cleaned OTU table est bien débarrassée des otus uniquement présentes dans les neg et des otu à enlever complètement (639-95-81=463 OTU restantes)

my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% as.data.frame() %>% 
  pull(1)-> total_reads_per_sample

hist(total_reads_per_sample)
summary(total_reads_per_sample)

# rm(ps.pa, ps.pa.neg, ps.pa.pos, ps.decon, ps.decon.neg, ps.decon.pos)

#**Taxonomie à l'ancetre commun sur la table d'OTU décontaminée----
my_tax %>% 
  filter(OTU %in% rownames(my_cleaned_otu_table)) -> my_cleaned_tax_table

#**Table de références---- 
distrib_OTU <- function(otu_tab){
  read_tsv(illumina_data,
           na = c("", "NA", "No_hit", "0", "0.0"),
           show_col_types = FALSE) %>% 
    filter(!is.na(references)) %>% 
    select(OTU, references) -> references_column
  
  otu_tab %>% 
    mutate(OTU = as.numeric(rownames(.))) %>% 
    left_join(y=references_column, by = "OTU") %>%
    pivot_longer(cols = -c("OTU", "references"), names_to = "sample_code", values_to = "reads") %>%
    filter(reads > 0) %>% 
    separate_rows(references, sep = ",") %>%
    group_by(OTU, sample_code) %>%
    arrange(., sample_code) %>% 
    mutate(corrected_reads = reads/n()) %>%
    ungroup() %>%
    mutate(code_unique = paste(OTU, references, sep="_")) %>% 
    select(-c(references,OTU,reads)) %>% 
    pivot_wider(names_from = sample_code, values_from = corrected_reads) %>% 
    replace(.,is.na(.),0) %>% 
    column_to_rownames(var = "code_unique")
}  

distrib_OTU(my_cleaned_otu_table) -> my_otu_distrib

#**Taxonomie distribuée----
distrib_tax <- function(otu_tab) {
  read_tsv(
    illumina_data,
    na = c("", "NA", "No_hit", "0", "0.0"),
    show_col_types = FALSE
  ) %>%
    select(OTU, taxonomy, references) %>%
    filter(OTU %in% rownames(otu_tab)) %>%
    separate_rows(references, sep = ",") %>%
    left_join(
      read_tsv(file = "databioinfo/SILVA_138.1_SSURef_tax_silva.acc2taxo.gz"),
      by = c("references" = "accession")
    ) %>%
    select(-taxonomy.x) %>%
    separate_rows(taxonomy.y, sep = "\\|") %>%
    left_join(
      x = .,
      y = taxa2ranks,
      by = c("taxonomy.y" = "taxonomic_path")
    ) %>%
    filter(taxonomic_rank %in% good_taxonomic_levels) %>%  # Enlève les NA (possiblement au niveau de l'espèce)
    mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>%
    pivot_wider(names_from = taxonomic_rank, values_from = taxonomy.y) %>%
    mutate(code_unique = paste(OTU, references, sep = "_")) %>%  # On ajoute code unique avec OTU_références
    select(-references) %>%
    column_to_rownames(var = "code_unique") # la colonne "code unique" devient le nom des lignes
}

distrib_tax(my_cleaned_otu_table) -> my_tax_distrib

#Il y a 0 ref assignées à autre chose que des bactéries
my_tax_distrib %>%
  filter(!domain == "Bacteria") %>%
  rownames()

#**Objets phyloseq ----
#Deux objets phyloseq avec uniquement les puis environnementaux : un avec OTU et taxonomie à l'ancêtre commun. L'autre avec references et taxonomie distribuée

my_mtd %>% 
  filter(type=="ENV") -> my_mtd_env 
sample_data(my_mtd_env) -> MTD

otu_table(my_cleaned_otu_table,taxa_are_rows = T) -> OTU
tax_table(as.matrix(my_cleaned_tax_table)) -> TAX

ps = phyloseq(OTU, TAX, MTD)


otu_table(my_otu_distrib,taxa_are_rows = T) -> OTU
tax_table(as.matrix(my_tax_distrib)) -> TAX

ps_multiple = phyloseq(OTU, TAX, MTD)

#Qualité du séquençage----
#**Profondeurs de séquençage----
#On ajoute les nombres de reads totaux de chaque puit dans le tableau des métadonnées
my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("Label") %>% 
  rename(total_reads = 2) %>%
  left_join(my_mtd_env, ., by = "Label") -> my_mtd_env

summary(my_mtd_env$total_reads)

#Plein de graphiques A TRIER
ggplot(my_mtd_env, aes(log10(total_reads), fill=locality))+
  geom_histogram()+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Total reads", y = "Frequency")

# ggsave("res/RUN4_Distribution_of_bacterial_reads.pdf")

my_mtd_env <- my_mtd_env %>% arrange(total_reads)
sample_order<-my_mtd_env$sample_ID
label_order<-my_mtd_env$Label

ggplot(my_mtd_env, aes(x=sample_ID, y=total_reads, col=env))+
  geom_point(size=3)+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)

# ggsave("res/RUN4_Distribution_of_total_reads.pdf")

ggplot(my_mtd_env, aes(x=population, y=total_reads, fill = population))+
  geom_bar(stat = "identity")
  
#Les profondeurs de séquençage n'ont pas l'air de varier selon env, sex ou locality. Peut-être un test pour vérifier ça ? On n'a pas beaucoup de libertés...

#**Proportions of unassigned reads----
columns_to_keep <- c("OTU", "references", only_environmental_samples)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(all_of(columns_to_keep)) %>%  
  mutate_at(., c(3:ncol(.)),~replace(.,is.na(.),0)) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>%
  filter(is.na(references)) -> unassigned 
#il n'y a qu'une seule OTU non assignée. OTU6. Si on la blaste on trouve une mitochondrie de Bactrocera. Pas très anormal. Elle a un total de 139 reads. Pas de quoi perturber l'analyse. On laisse tomber l'analyse des non assignées.

#*******************----
#Courbes de raréfaction----
library("iNEXT")

load("results_fruit/diversities.Rdata")

#Selon le nombre d'ASV
out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:72]) %>%
  filter(Order.q == 0) %>% 
  ggplot(., aes(x = m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

#Selon l'indice de Shannon
out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:72]) %>%
  filter(Order.q == 1) %>%
  filter(m<=50000) %>% 
  ggplot(., aes(x = m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

#Par échantillons, en fonction des nombres de Hill
out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:12]) %>% 
  ggplot(., aes(x = m, y = qD, col = as.factor(Order.q)))+
  geom_line()+
  facet_wrap(~Assemblage, ncol=4)+
  xlab("Number of reads") +
  ylab("Number of ASV") +
  theme_minimal()+
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

#*******************----
#BARPLOT - Composition bactérienne par échantillon----
#Une façon de faire basée sur les abondances relatives
as.data.frame(my_tax_distrib) %>%  #get taxonomy
  rownames_to_column(var = "tax_ID") -> df 
as.data.frame(my_otu_distrib) %>% #get reads 
  rownames_to_column(var = "tax_ID") %>% 
  left_join(df, ., by="tax_ID") %>% 
  select(-tax_ID) %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  left_join(., my_mtd_env, by = c("Sample" = "Label")) -> df

list_of_population <- c("ANS", "ETA","REL_GB","VER","GUI","MAR","REL_GF","PAL", "VID","BER", "BOU", "CIR")

# Composition par population
palette_16 <- c("#3E4A63", "#E56B6F", "#63B69F", "#F6AE2D", "#9B5DE5", "#F96E46", "#5B9BD5", "#E75D6F", "#4BAF7F", "#F5CA5C", "#7C4D79", "#E89644", "#3F6849", "#C98B9D", "#718CA1", "#E2979C")

df %>% 
  count(family,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>%
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  mutate(env_lab = paste0(substr(Sample, 1, 2), "_", env)) %>% 
  mutate(population = factor(population, levels = list_of_population)) %>%
  ggplot(., aes(x=env_lab, y = frq, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))+
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Echantillons environnementaux", y = "Abondances relatives", fill="Familles")+
  scale_fill_manual(values = palette_16) +
  facet_wrap(~population, scales="free_x", ncol = 12)

#Par fruit hôte 
df %>% 
  count(family,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>%
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  ggplot(., aes(x=Sample, y = frq, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2))+
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Echantillons environnementaux", y = "Abondances relatives", fill="Familles")+
  scale_fill_manual(values = palette_16)+
  facet_wrap(~env, scales="free_x", ncol = 5)

#Par population et genre
colors <- palette("set 3")

df %>% 
  count(genus,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  mutate(genus = ifelse(all(frq < 0.05), "Autres", genus)) %>% 
  ungroup() %>%
  mutate(genus = replace(genus, genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium")) %>%
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  mutate(env_lab = paste0(substr(Sample, 1, 2), "_", env)) %>% 
  mutate(population = factor(population, levels = list_of_population)) %>%
  ggplot(., aes(x=env_lab, y = frq, fill = genus)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 11))+
  guides(fill = guide_legend(ncol = 1))+
  labs(x="Echantillons environnementaux", y = "Abondances relatives", fill="Genres")+
  scale_fill_manual(values = colors) +
  facet_wrap(~population, scales="free_x", ncol = 12)

# Abondance des familles par échantillon
df %>% 
  count(family,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>% 
  group_by(family, Sample) %>% 
  summarise(frq = sum(frq)) %>% 
  ungroup() %>% 
  mutate(frq = frq*100) %>% 
  # filter(family == "Enterobacteriaceae") %>%
  pivot_wider(
    id_cols = 1,
    names_from = Sample,
    values_from = frq,
    values_fill = NA
  ) %>%
  view()

# Abondance des genres par échantillon
df %>% 
  count(genus,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(genus) %>% 
  mutate(genus = ifelse(all(frq < 0.05), "Autres", genus)) %>% 
  ungroup() %>% 
  group_by(genus, Sample) %>% 
  summarise(frq = sum(frq)) %>% 
  ungroup() %>% 
  mutate(frq = frq*100) %>% 
  # filter(genus == "Baccilus") %>%
  pivot_wider(
    id_cols = 1,
    names_from = Sample,
    values_from = frq,
    values_fill = NA
  ) %>%
  view()

# Abondance par localité  
df %>% 
  count(family,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>% 
  mutate(Ab = frq * sum(total)) %>% 
  group_by(family, Sample) %>% 
  summarise(Ab = sum(Ab)) %>% 
  ungroup() %>% 
  # mutate(frq = frq*100) %>% 
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  count(family,locality, wt = Ab, name = "loca_tot") %>% 
  # left_join(., my_mtd_env, by = "locality") %>% View()
  group_by(locality) %>% 
  mutate(loca_tot = loca_tot/sum(loca_tot)) %>% 
  # filter(locality == "MAR") %>% 
  # mutate(test = sum(loca_tot)) %>% View()
  ungroup() %>% 
  pivot_wider(
    id_cols = 1,
    names_from = locality,
    values_from = loca_tot,
    values_fill = NA
  ) %>% View()

# Abondance par fruit
df %>% 
  count(family,Sample, wt = Abundance, name = "total") %>% 
  group_by(Sample) %>% 
  mutate(frq = total/sum(total)) %>% 
  ungroup() %>% 
  group_by(family) %>% 
  mutate(family = ifelse(all(frq < 0.05), "Autres", family)) %>% 
  ungroup() %>% 
  mutate(Ab = frq * sum(total)) %>% 
  group_by(family, Sample) %>% 
  summarise(Ab = sum(Ab)) %>% 
  ungroup() %>% 
  # mutate(frq = frq*100) %>% 
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  count(family,locality, wt = Ab, name = "loca_tot") %>% 
  # left_join(., my_mtd_env, by = "locality") %>% View()
  group_by(locality) %>% 
  mutate(loca_tot = loca_tot/sum(loca_tot)) %>% 
  # filter(locality == "MAR") %>% 
  # mutate(test = sum(loca_tot)) %>% View()
  ungroup() %>% 
  pivot_wider(
    id_cols = 1,
    names_from = locality,
    values_from = loca_tot,
    values_fill = NA
  ) %>% View()

# Afficher les abondances par taxon
#Niveau genre
df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  left_join(., my_mtd_env, by = c("Sample" = "Label")) %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() -> df_with_NA

df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i) {
    paste(domain[i], phylum[i], class[i], order[i], family[i], genus[i], sep =
            "_")
  })) %>%
  group_by(tax_ID, Sample) %>%
  summarise(Rel_ab = sum(Rel_ab)) %>%
  pivot_wider(
    id_cols = 1,
    names_from = Sample,
    values_from = Rel_ab,
    values_fill = NA
  ) %>%
  view()

#Niveau famille
df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i) {
    paste(domain[i], phylum[i], class[i], order[i], family[i], sep = "_")
  })) %>%
  group_by(tax_ID, Sample) %>%
  summarise(Rel_ab = sum(Rel_ab)) %>%
  pivot_wider(
    id_cols = 1,
    names_from = Sample,
    values_from = Rel_ab,
    values_fill = NA
  ) %>% view()

#*******************----
#NMDS----
#Avec matrice présence/absence
#Env/sex
library(ade4)
set.seed(1)
ps_rare = rarefy_even_depth(ps, sample.size = 24000)
ps_bin = transform_sample_counts(ps_rare, function(x) 1*(x>0))

jaccard <- ordinate(ps_bin, "NMDS", "jaccard") #stress trop haut
jaccard <- ordinate(ps_bin, "PCoA", "jaccard") 

col = c("purple", "darkgreen", "red","gray","orange")

plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "env",
                shape = "sex") +
  theme_bw() +
  geom_point(size = 5) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm")) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Sexes", values=c(17,19),
                     labels = c("Femelle", "Mâle")) +
  scale_color_manual(name = "Fruits hôtes", values = col) +
  geom_text_repel(aes(label = locality, fontface = 'italic'), size = 5,max.overlaps = 65)

#Env/Locality
col = c("red","Blue","Yellow","Orange","Green")

plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "locality",
                shape = "env") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "host fruit", values=c(1:5)) +
  # scale_color_manual(name = "locality") +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

ggsave("results/NMDS_flies_plot.svg", width=15, height = 10, units = "cm")

#T°/plu
col = c("red","blue","darkgreen")

plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "grp_t",
                shape = "grp_plu") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Pluviosité", values = c(15,16,17),
                     labels = c("Faible","Moyenne","Forte")) +
  scale_color_manual(name = "Température", values = col,
                     labels = c("Elevée","Basse","Moyenne")) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

# ggsave("resultS_fruit/NMDS_flies_plot.svg", width=15, height = 10, units = "cm")

#T°/Alt
col = c("red","blue","darkgreen")

plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "grp_alt",
                shape = "grp_plu") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Pluviosité", values = c(15,16,17),
                     labels = c("Faible","Moyenne","Forte")) +
  scale_color_manual(name = "Altitude", values = col,
                     labels = c("Elevée","Basse","Moyenne")) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

# ggsave("resultS_fruit/NMDS_flies_plot.svg", width=15, height = 10, units = "cm")

#taxons
library(paletteer)
ncol = 3
col_tax <- paletteer_d("ggsci::default_igv")[seq(ncol, ncol + 17, 1)]

unique(my_cleaned_tax_table$class) 

plot_ordination(ps_bin,
                jaccard,
                type = "taxa",
                color = "class") +
  theme_bw() +
  geom_point(size = 2) +
  theme(legend.position = "right") +
  scale_color_manual(values = col_tax) +
  #  scale_shape_manual(values = c(19, 17, 3)) +
  geom_text_repel(aes(label = family, fontface = 'italic'), size = 2.5)

ggsave("results/NMDS_bacteria_plot.svg", width = 15,  height=10, units = "cm")

#Permanova----
# Sur données non métrique 
df <- as.data.frame(t(otu_table(ps_rare)))
my_mtd_env %>%
  filter(Label != c("09_B-ETA-male","24_G-VER-fem")) -> mtd_perm
set.seed(1)
Permanovaresults <- adonis2(df~population*sex, data = mtd_perm,
                           method = "bray", permutation = 9999)
Permanovaresults

Permanovaresults <- adonis2(df~env*grp_t*grp_plu, data = mtd_perm,
                            method = "bray", permutation = 9999)
Permanovaresults <- adonis2(df~env+grp_t+grp_plu, data = mtd_perm,
                            method = "bray", permutation = 9999)
Permanovaresults
#J'ai pris ces résultats pour le moment qui semble mieux car pas de transformation de données

#Permanova
#Transformation des données 
df <- as.data.frame(t(otu_table(ps_rare))) %>% 
  sqrt()
my_mtd_env %>%
  filter(Label != c("09_B-ETA-male","24_G-VER-fem")) -> mtd_perm
set.seed(1)
Permanovaresults <- adonis2(df~population*sex, data = mtd_perm,
                            method = "euclidean", permutation = 9999)
Permanovaresults

Permanovaresults <- adonis2(df~env*grp_t*grp_plu, data = mtd_perm,
                            method = "euclidean", permutation = 9999)
Permanovaresults
Permanovaresults <- adonis2(df~env+grp_t+grp_plu, data = mtd_perm,
                            method = "euclidean", permutation = 9999)
Permanovaresults

#{{extract permanova R2 and p-values}}
perm_R2 <- Permanovaresults$R2
perm_pvalues <- Permanovaresults$`Pr(>F)`

#NMDS abondance
#env/localité

# ps_rare = rarefy_even_depth(ps, sample.size = 24000)

bray <- ordinate(ps_rare, "NMDS", "bray") #stress élevé

bray <- ordinate(ps_rare, "PCoA", "bray") 

plot_ordination(ps_rare,
                bray,
                type = "samples",
                shape = "env", 
                col = "locality") +
  theme_bw() +
  geom_point(size = 2.5)+
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Host fruit", values=c(15:19)) +
  # scale_colour_manual(name = "Locality", values = col_env) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

#env/sex

col = c("purple", "darkgreen", "red","gray","orange")

plot_ordination(ps_rare,
                bray,
                type = "samples",
                color = "env",
                shape = "sex") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Sexes", values=c(17,19),
                     labels = c("Femelle", "Mâle")) +
  scale_color_manual(name = "Fruits hôtes", values = col) +
  geom_text_repel(aes(label = locality, fontface = 'italic'), size = 3,max.overlaps = 65)

#env/sex
eucli <- ordinate(ps_rare, "PCoA", "euclidean") 

col = c("purple", "darkgreen", "red","gray","orange")

plot_ordination(ps_rare,
                eucli,
                type = "samples",
                color = "env",
                shape = "sex") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Sexes", values=c(17,19),
                     labels = c("Femelle", "Mâle")) +
  scale_color_manual(name = "Fruits hôtes", values = col) +
  geom_text_repel(aes(label = locality, fontface = 'italic'), size = 3,max.overlaps = 65)

#T°/plu
col = c("red","blue","darkgreen")

sample_data(ps_rare) %>% View()

plot_ordination(ps_rare,
                bray,
                type = "samples",
                color = "grp_t",
                shape = "grp_plu") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10)) +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "Pluviosité", values = c(15,16,17),
                     labels = c("Faible","Moyenne","Forte")) +
  scale_color_manual(name = "Température", values = col,
                     labels = c("Elevée","Basse","Moyenne")) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

#taxon

plot_ordination(ps_rare,
                bray,
                type = "taxa",
                color = "phylum") +
  theme_bw() +
  geom_point(size = 2) +
  theme(legend.position = "right") +
  #  scale_color_manual(values = c(col_tax[1], col_tax[3], col_tax[2], col_tax[4:23])) +
  #  scale_shape_manual(values = c(19, 17, 3)) +
  geom_text_repel(aes(label = family, fontface = 'italic'), size = 2.5)

#*******************----
#Analyse de diversité----
load("results_fruit/diversities.Rdata")

div_df<-out$AsyEst[out$AsyEst$Diversity=="Shannon diversity",]

div_df<-left_join(div_df, my_mtd, by = c("Assemblage" = "Label"))

#div alpha par site 
div_df %>% 
  group_by(locality) %>% 
  mutate(.before = 4, "mean_loca" = mean(Estimator)) %>% 
  distinct(mean_loca, .keep_all = TRUE) %>% View()

#div alpha par env
div_df %>% 
  group_by(env) %>% 
  mutate(.before = 4, "mean_env" = mean(Estimator)) %>%
  distinct(mean_env, .keep_all = TRUE) %>% View

#div alpha par temp
div_df %>% 
  group_by(grp_t) %>% 
  mutate(.before = 4, "mean_t" = mean(Estimator)) %>%
  distinct(mean_t, .keep_all = TRUE) %>% 
  select(c(mean_t, join_Tmoy, grp_t)) %>% View()

#div alpha par pluv
div_df %>% 
  group_by(grp_plu) %>% 
  mutate(.before = 4, "mean_pluv" = mean(Estimator)) %>%
  distinct(mean_pluv, .keep_all = TRUE) %>%
  select(c(mean_pluv, join_pluie_moy, grp_plu)) %>% View()

#div alpha par populations
div_df %>% 
  group_by(population) %>% 
  mutate(.before = 4, "mean_pop" = mean(Estimator)) %>%
  distinct(mean_pop, .keep_all = TRUE) %>% View()

#div alpha par sexe
div_df %>% 
  group_by(sex) %>% 
  mutate(.before = 4, "mean_sex" = mean(Estimator)) %>%
  distinct(mean_sex, .keep_all = TRUE) %>% View()

#Box plot locality
colors <- c("#FF5F5F", "#FF8C42", "#FFC168", "#FFEE8D", "#D1FFA9", "#A4FFC5", "#77FFD2", "#4AFFD9", "#1DFFFF", "#79E7FF", "#A4B3FF")

ggplot(div_df, aes(
  x = locality,
  y = Estimator,
  fill = locality
)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  ) +
  scale_x_discrete(name = "Echantillons environnementaux") +
  scale_y_continuous(name = "Diversité alpla") +
  scale_fill_manual(
    name = "Localités",
    values = colors,
    labels = c("ANS","BER","BOU","CIR","ETA","GUI","MAR","PAL","REL","VER","VID")
  )+
  xlab("Generations in the laboratory") + ylab("Alpha diversity (in genus equivalents)")

#Box plot population
colors <- c("#FF5F5F", "#FF8C42", "#FFC168", "#FFEE8D", "#D1FFA9", "#A4FFC5", "#77FFD2", "#4AFFD9", "#1DFFFF", "#79E7FF", "#A4B3FF","#F3A0F2")

loca <- ggplot(div_df, aes(
  x = population,
  y = Estimator,
  fill = population
)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  ) +
  scale_x_discrete(name = "Echantillons environnementaux") +
  scale_y_continuous(name = "Diversité alpla") +
  scale_fill_manual(
    name = "Population",
    values = colors
  )

#Box plot sexe
colors <- c("#FF5F5F","#1f77b4")

sexe <- ggplot(div_df, aes(
  x = sex,
  y = Estimator,
  fill = sex
)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  )+
  scale_fill_manual(name = "Sexe",
                    values = colors,
                    labels = c("Femelles","Mâles"))+
  scale_x_discrete(name = "Echantillons environnementaux",
                   labels = c("Femelles","Mâles")) 
scale_y_continuous(name = "Diversité alpla") 

library(cowplot)
plot_grid(loca,sexe,labels = c("Population","Sexe"), ncol = 2)

#Box plot env
env <- ggplot(div_df, aes(
  x = env,
  y = Estimator,
  fill = env
)) +
  geom_boxplot() +
  # geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.5,
  #               position=position_dodge(.9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  ) +
  scale_x_discrete(name = "Echantillons environnementaux") +
  scale_y_continuous(name = "Diversité alpla") +
  scale_fill_manual(
    name = "Fruit hôte",
    values = c("#F15F79", "#FFB347", "#77DD77", "#AEC6CF", "#F3A0F2"),
    labels = c("Badamier", "Goyavier blanc", "Goyavier fraise", "Jamrosat", "Mangue")
  )

#Box plot temp
colors <- c("#FF5F5F","#FFEE8D","#A4FFC5", "#77FFD2", "#4AFFD9", "#1DFFFF", "#79E7FF", "#A4B3FF","#F3A0F2")

temp <- ggplot(div_df, aes(
  x = grp_t,
  y = Estimator,
  fill = grp_t
)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  )+
  scale_fill_manual(name = "Température",
                    values = colors,
                    labels = c("Elevée","Basse","Moyenne")) +
  scale_x_discrete(name = "Echantillons environnementaux",
                   labels = c("Elevée","Basse","Moyenne")) +
  scale_y_continuous(name = "Diversité alpla") 

#Box plot pluv
colors <- c("#79E7FF", "#A4B3FF","#F3A0F2")

pluv <- ggplot(div_df, aes(
  x = grp_plu,
  y = Estimator,
  fill = grp_plu
)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  )+
  scale_fill_manual(name = "Pluviosité",
                    values = colors,
                    labels = c("Faible","Moyenne","Forte")) +
  scale_x_discrete(name = "Echantillons environnementaux",
                   labels = c("Faible","Moyenne","Forte")) +
  scale_y_continuous(name = "Diversité alpla") 


plot_grid(env,pluv,temp, labels =c("Fruit hôte","Pluviosité","Température"), ncol = 3)

#*******************----
#CLUSTERING FROM READ COUNTS----
library("sbm")

set.seed(1)
ps_rare = rarefy_even_depth(ps, sample.size = 24000)
hist(log10(taxa_sums(ps_rare)+1))
Log_FB <- log10(t(otu_table(ps_rare)+1))

plotMyMatrix(Log_FB, dimLabels = c(row = "Samples", col = "ASV"))

#Estimation du LBM
LBM_log_reads_gaussian <-
  as.matrix(Log_FB) %>%
  estimateBipartiteSBM(
    model = 'gaussian')

load("results_fruit/LBM_log_reads_gaussian24000_18_26.Rdata")

#Représentation graphique de la matrice----
memb_spl_obs <- LBM_log_reads_gaussian$memberships$row
memb_bac_obs <- LBM_log_reads_gaussian$memberships$col

df=as.data.frame(t(Log_FB))
df$taxon=rownames(df)
df_long<-pivot_longer(df, cols=names(df)[-ncol(df)], names_to =
                        "sample_ID", values_to = "reads")

bac_list<-colnames(Log_FB)
sample_list<-rownames(Log_FB)

LBM_log_reads_gaussian$nbBlocks #how many clusters
Qspl<-LBM_log_reads_gaussian$nbBlocks["row"][[1]]
Qasv<-LBM_log_reads_gaussian$nbBlocks["col"][[1]]

LBM_log_reads_gaussian$memberships #cluster composition
lapply(1:Qspl,function(q){sample_list[LBM_log_reads_gaussian$memberships$row
                                      == q]})
lapply(1:Qasv,function(q){bac_list[LBM_log_reads_gaussian$memberships$col == q]}) -> compo_clusters

df_tax_cluster=data.frame(taxon=bac_list,tax_cluster=LBM_log_reads_gaussian$memberships$col)
df_long<-left_join(df_long, df_tax_cluster, by="taxon")

df_long %>% arrange(tax_cluster) -> df_long
correct_taxon_order <- unique(df_long$taxon)

df_spl_cluster=data.frame(sample_ID=sample_list,spl_cluster=LBM_log_reads_gaussian$memberships$row)
df_long<-left_join(df_long, df_spl_cluster, by="sample_ID")

df_long %>% arrange(spl_cluster) -> df_long
correct_sample_order <- unique(df_long$sample_ID)

ggplot(df_long,aes(x=taxon, y=sample_ID, fill=reads))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y = element_text(vjust = 0, hjust = 1, size = 6),
        axis.text.x = element_text(angle = 90, size = 6))+
  scale_fill_gradient(low="white", high="Black")+
  scale_x_discrete(limits = correct_taxon_order, label =
                     my_tax$family[match(my_tax$OTU, colnames(FB_bin))])+
  scale_y_discrete(limits = correct_sample_order)+
  geom_hline(yintercept=0.5+cumsum(table(LBM_log_reads_gaussian$memberships$row))[-length(table(LBM_log_reads_gaussian$memberships$row))],
             color = "red")+
  # geom_vline(xintercept=0.5+cumsum(table(LBM_log_reads_gaussian$memberships$col))[-length(table(LBM_log_reads_gaussian$memberships$col))],
  #            color = "red")+
  ylab("Echantillons")+
  xlab("Bactéries")

#Exploration clustering
my_tax_distrib %>% 
  filter(OTU %in% compo_clusters[[1]]) %>% View()

my_tax_distrib %>% 
  filter(OTU %in% compo_clusters[[1]]) %>% 
  filter(family == "Enterobacteriaceae") %>% View()

my_tax_distrib %>% 
  filter(OTU %in% compo_clusters[[1]]) %>% 
  filter(genus == "Klebsiella") %>% View()

my_tax_distrib %>% 
  filter(OTU %in% compo_clusters[[1]]) %>% 
  filter(genus == "Enterobacter") %>% View()

my_cleaned_otu_table["1",]

#--------- Null model 
set.seed(1)
ps_rare = rarefy_even_depth(ps, sample.size = 24000)
hist(log10(taxa_sums(ps_rare)+1))

#Cluster randomize
library(parallel)

set.seed(1)
nbCores <- 1
n_iterations <- 10000
Log_FB <- log10(t(otu_table(ps_rare)+1))

sum_reads <- sum(Log_FB)
sum_reads_per_OTU <- apply(Log_FB, 1, sum)
sum_reads_per_samples <- apply(Log_FB, 2, sum)
nrow_log_FB <- nrow(Log_FB)

sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads %>% str()

## Null model ------

randomize.WEDD <- function(OTU_matrix_t) {
  apply(sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads, 
        c(1, 2), function(x)
          rnorm(1, mean = x, 
                sd = sd(sum_reads_per_OTU %o% sum_reads_per_samples / sum_reads)))
}

simulate_membership <- function(OTU_matrix_t) {
  sbm::estimateBipartiteSBM(netMat = OTU_matrix_t,
                            model = "gaussian",
                            estimOptions = list(
                              verbosity = 0,
                              plot = FALSE,
                              nbCores = nbCores)) -> LBM_sim
  
  LBM_sim$memberships$row
}

randomize_and_simulate <- function(n) {
  #argument n non utilisé mais obligatoire pour mclapply
  Log_FB %>%
    randomize.WEDD(.) %>%
    simulate_membership(.)
}

nrow_log_FB <- nrow(FB_bin)

system.time(
  parallel::mclapply(1:n_iterations,
                     randomize_and_simulate,
                     mc.cores = (parallel::detectCores()/2-1)) %>% 
    unlist(.) %>%
    matrix(., nrow = nrow_log_FB , byrow = TRUE) -> clust
)
rm(sum_reads, sum_reads_per_OTU, sum_reads_per_samples)

# save(clust, file = "results_fruit/clust_NPERM10000_log.Rdata")

#------CCA------
#if not alreay loaded : load("results_fruit/clust_NPERM10000_log.Rdata")
library(lme4)
source('scriptENV/scriptsFRUIT/fruit 2023/functions_modified.R') # Kindly provided by Massol et al. 2021 JAE

# Factors
env_vect <- my_mtd_env$env[match(rownames(Log_FB), my_mtd_env$Label)]
sex_vect <- my_mtd_env$sex[match(rownames(Log_FB), my_mtd_env$Label)]
loc_vect <- my_mtd_env$locality[match(rownames(Log_FB), my_mtd_env$Label)]

analysis.function.alt_wgt <-
  function(obs_memb, var1, var2, var3, configs = NULL) {
    #industrial processing of CCA, two permutations
    n <- length(obs_memb)
    modules <- dummy(as.factor(obs_memb), levelsToKeep = levels(as.factor(obs_memb)))
    one.1 <- cca(modules ~ var1, na.action = na.omit)
    one.2 <- cca(modules ~ var2, na.action = na.omit)
    one.3 <- cca(modules ~ var3, na.action = na.omit)
    two.1 <- cca(modules ~ var1 + var2, na.action = na.omit)
    two.2 <- cca(modules ~ var1 + var3, na.action = na.omit)
    two.3 <- cca(modules ~ var2 + var3, na.action = na.omit)
    three <- cca(modules ~ var1 + var2 + var3, na.action = na.omit)
    one.1.alone <- cca(modules ~ var1 + Condition(var2) + Condition(var3), na.action = na.omit)
    one.2.alone <- cca(modules ~ Condition(var1) + var2 + Condition(var3), na.action = na.omit)
    one.3.alone <- cca(modules ~ Condition(var1) + Condition(var2) + var3, na.action = na.omit)
    one.1.vs.3 <- cca(modules ~ var1 +  Condition(var3), na.action = na.omit)
    one.1.vs.2 <- cca(modules ~ var1 +  Condition(var2), na.action = na.omit)
    one.2.vs.3 <- cca(modules ~ var2 +  Condition(var3), na.action = na.omit)
    one.2.vs.1 <- cca(modules ~ var2 +  Condition(var1), na.action = na.omit)
    one.3.vs.1 <- cca(modules ~ var3 +  Condition(var1), na.action = na.omit)
    one.3.vs.2 <- cca(modules ~ var3 +  Condition(var2), na.action = na.omit)
    xxx <-list(one.1,one.2,one.3,two.1,two.2, two.3,three,
               one.1.alone, one.2.alone, one.3.alone, one.1.vs.3,one.1.vs.2,one.2.vs.3,
               one.2.vs.1, one.3.vs.1, one.3.vs.2)
    chi2 <- sapply(1:16, function(x) chi2.computation(xxx[[x]]))
    Forms <- sapply(1:16, function(x) extractTerms(xxx[[x]]))
    F <- sapply(1:16, function(x) F.computation(xxx[[x]], n))
    P <- sapply(1:16, function(x) anova(xxx[[x]], permutations = how(nperm = NPERM - 1))$P[1])
    if (is.null(configs)) {
      results <- data.frame(
        "Formulas" = Forms,
        "chi2" = chi2,
        "F" = F,
        "P" = P
      )}
    else {
      depth <- dim(configs)[2]
      all.F <-
        sapply(1:depth, function(x)
          only.F.alt(configs[, x], var1, var2, var3))
      all.ecdf <- apply(all.F, 1, ecdf)
      P2 <- sapply(1:16, function(x)
        1 - ((all.ecdf[[x]])(F[x])))
      results <-
        data.frame(
          "Formulas" = Forms,
          "chi2" = chi2,
          "F" = F,
          "P" = P,
          "P2" = P2
        )
    }
    results
  }

cca.data <-
  data.frame(
    env = as.factor(env_vect),
    sex = as.factor(sex_vect),
    loc = as.factor(loc_vect)
  )

wholenetwork.cca.results.NULL <- analysis.function.alt_wgt(memb_spl_obs,
                                                           cca.data$env,
                                                           cca.data$sex,
                                                           cca.data$loc,
                                                           configs = NULL)

wholenetwork.cca.results.CONFIGS <- analysis.function.alt_wgt(memb_spl_obs,
                                                              cca.data$env,
                                                              cca.data$sex,
                                                              cca.data$loc,
                                                              configs = clust)

# save(wholenetwork.cca.results.NULL, wholenetwork.cca.results.CONFIGS, file = "results/CCA_weighted.Rdata")
write.csv2(wholenetwork.cca.results.CONFIGS, file = "cca_weighted.csv")

#CCA sur loc et sexe----
#if not alreay loaded : load("results_fruit/clust_NPERM10000_log.Rdata")

# Factors
sex_vect <- my_mtd_env$sex[match(rownames(FB_bin), my_mtd_env$Label)]

my_mtd_env$loc_code = sapply(1:nrow(my_mtd_env), function(i) {
  fruit_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][1]
  site_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][2]
  paste(fruit_code, site_code, sep = "-")
})
my_mtd_env$loc_code %>%
  str_replace("C-", "GF-") %>% str_replace("G-", "GB-") ->
  my_mtd_env$loc_code

loc_vect_bis <- my_mtd_env$loc_code[match(rownames(FB_bin),
                                          my_mtd_env$Label)]

dummyMod<-dummy(as.factor(memb_spl_obs), levelsToKeep = levels(as.factor(memb_spl_obs)))

#CCA
cca.data <-
  data.frame(
    "sex" = as.factor(sex_vect),
    "loc" = as.factor(loc_vect_bis)
  )

analysis.function.alt_wgt2 <-
  function(obs_memb, var1, var2, configs = NULL) {
    #industrial processing of CCA, two permutations
    n <- length(obs_memb)
    modules <- dummy(as.factor(obs_memb), levelsToKeep =
                       levels(as.factor(obs_memb)))
    one.1 <- cca(modules ~ var1, na.action = na.omit)
    one.2 <- cca(modules ~ var2, na.action = na.omit)
    two.1 <- cca(modules ~ var1 + var2, na.action = na.omit)
    one.1.vs.2 <- cca(modules ~ var1 +  Condition(var2), na.action =
                        na.omit)
    one.2.vs.1 <- cca(modules ~ var2 +  Condition(var1), na.action =
                        na.omit)
    xxx <-list(one.1,one.2,two.1,one.1.vs.2,one.2.vs.1)
    chi2 <- sapply(1:5, function(x) chi2.computation(xxx[[x]]))
    Forms <- sapply(1:5, function(x) extractTerms(xxx[[x]]))
    F <- sapply(1:5, function(x) F.computation(xxx[[x]], n))
    P <- sapply(1:5, function(x) anova(xxx[[x]], permutations =
                                         how(nperm = NPERM - 1))$P[1])
    if (is.null(configs)) {
      results <- data.frame(
        "Formulas" = Forms,
        "chi2" = chi2,
        "F" = F,
        "P" = P
      )}
    else {
      depth <- dim(configs)[2]
      all.F <-
        sapply(1:depth, function(x)
          only.F.alt2(configs[, x], var1, var2))
      all.ecdf <- apply(all.F, 1, ecdf)
      P2 <- sapply(1:5, function(x)
        1 - ((all.ecdf[[x]])(F[x])))
      results <-
        data.frame(
          "Formulas" = Forms,
          "chi2" = chi2,
          "F" = F,
          "P" = P,
          "P2" = P2
        )
    }
    results
  }

wholenetwork.cca.results.NULL <- analysis.function.alt_wgt2(memb_spl_obs,
                                                            cca.data$sex,
                                                            cca.data$loc,
                                                            configs = NULL)

wholenetwork.cca.results.CONFIGS <- analysis.function.alt_wgt2(memb_spl_obs,
                                                               cca.data$sex,
                                                               cca.data$loc,
                                                               configs = clust)
write.csv2(wholenetwork.cca.results.CONFIGS, file = "results_fruit/cca_weighted_loca_sexe_2.csv")

#Venn (CCA)
varp <- varpart(dummyMod, cca.data$sex, cca.data$loc, chisquare = T)

plot(varp, Xnames = list("Sexe", "Population"))

#CCA (env, temp, pluv)----
env_vect <- my_mtd_env$env[match(rownames(Log_FB), my_mtd_env$Label)]
temp_vect <- my_mtd_env$grp_t[match(rownames(Log_FB), my_mtd_env$Label)]
pluv_vect <- my_mtd_env$grp_plu[match(rownames(Log_FB), my_mtd_env$Label)]

cca.data <-
  data.frame(
    "env" = as.factor(env_vect),
    "temp" = as.factor(temp_vect),
    "pluv" = as.factor(pluv_vect)
  )

#Get the corresponding null distribution of clustering
clust -> cle.config

#CCA with line randomization
wholenetwork.cca.results.NULL <- analysis.function.alt(as.data.frame(Log_FB),
                                                       memb_spl_obs,
                                                       cca.data$env,
                                                       cca.data$temp,
                                                       cca.data$pluv,
                                                       configs = NULL)
#CCA with full randomization
wholenetwork.cca.results.CONFIGS <- analysis.function.alt(as.data.frame(Log_FB),
                                                          memb_spl_obs,
                                                          cca.data$env,
                                                          cca.data$temp,
                                                          cca.data$pluv,
                                                          configs = cle.config)
write.csv2(wholenetwork.cca.results.CONFIGS, file = "results_fruit/cca_weighted_tplu2.csv")


#Venn (CCA TPE)
dummyMod<-dummy(as.factor(memb_spl_obs), levelsToKeep = levels(as.factor(memb_spl_obs)))

varp <- varpart(dummyMod, cca.data$env, cca.data$pluv,cca.data$temp,  chisquare = T)

plot(varp, Xnames = list("Fruit hôte", "Pluviosité","Température"))

#------NMI----
library(igraph)
load("results_fruit/clust_NPERM10000_log.Rdata")

# Population
my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("loc_code") %>% 
  as.factor() -> vect

memb_spl_obs <- LBM_log_reads_gaussian$memberships$row

nmi_obs <- compare(as.factor(memb_spl_obs), vect, method = "nmi")

clust %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p1 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Population")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

# température
my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("grp_t") %>% 
  as.factor() -> vect 

memb_spl_obs <- LBM_log_reads_gaussian$memberships$row

nmi_obs <- compare(as.factor(memb_spl_obs), vect, method = "nmi")

clust %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p2 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2)+ 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Température")


print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

# pluvio
my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("grp_plu") %>% 
  as.factor() -> vect 

nmi_obs <- compare(as.factor(memb_spl_obs), vect, method = "nmi")

clust %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p3 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2)+ 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Pluviométrie")


print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

# Env
my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("env") %>% 
  as.factor() -> vect 

nmi_obs <- compare(as.factor(memb_spl_obs), vect, method = "nmi")

clust %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p4 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2)+ 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Fruit hôte")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

# sex
my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("sex") %>% 
  as.factor() -> vect 

nmi_obs <- compare(as.factor(memb_spl_obs), vect, method = "nmi")

clust %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p5 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2)+ 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Sexe")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

library(cowplot)
plot_grid(p1, p5, blankPlot ,p2, p3, p4, labels=c("A", "B"," ","C","D","E"), ncol = 3, nrow = 2)

#*******************----
#CLUSTERING FROM PRESENCE ABSENCE DATA----

#***LBM (Latent Block Model = modèle de détection des commu) LBM pour matrice binaire ???----
library(sbm)

LBM_bernoulli <-
  as.matrix(FB_bin) %>% 
  estimateBipartiteSBM(
    model = 'bernoulli', 
    dimLabels = c(row = "Samples", col = "ASV"))

memb_spl_obs_ber <- LBM_bernoulli$memberships$Samples
memb_bac_obs_ber <- LBM_bernoulli$memberships$ASV

# save(LBM_bernoulli, ps_rare, file = "results_fruit/LBM_bernouilli24000.Rdata")
load("results_fruit/LBM_bernouilli24000.Rdata")

#Graph cluster----
library(igraph)
NPERM = 10000 # Number of configuration matrices for each rarefied matrix

set.seed(1)
ps_rare = rarefy_even_depth(ps, sample.size = 24000) 
ps_bin = transform_sample_counts(ps_rare, function(x) 1*(x>0))

FB_bin <- as.data.frame(t(otu_table(ps_bin)))
df=as.data.frame(t(FB_bin))
df$taxon=rownames(df)
df_long<-pivot_longer(df, cols=names(df)[-ncol(df)], names_to =
                        "sample_ID", values_to = "reads")

bac_list<-colnames(FB_bin)
sample_list<-rownames(FB_bin)

LBM_bernoulli$nbBlocks #how many clusters
Qspl<-LBM_bernoulli$nbBlocks["Samples"][[1]]
Qasv<-LBM_bernoulli$nbBlocks["ASV"][[1]]

LBM_bernoulli$memberships #cluster composition
lapply(1:Qspl,function(q){sample_list[LBM_bernoulli$memberships$Samples
                                      == q]})
lapply(1:Qasv,function(q){bac_list[LBM_bernoulli$memberships$ASV == q]})

df_tax_cluster=data.frame(taxon=bac_list,tax_cluster=LBM_bernoulli$memberships$ASV)
df_long<-left_join(df_long, df_tax_cluster, by="taxon")

df_long %>% arrange(tax_cluster) -> df_long
correct_taxon_order <- unique(df_long$taxon)

df_spl_cluster=data.frame(sample_ID=sample_list,spl_cluster=LBM_bernoulli$memberships$Samples)
df_long<-left_join(df_long, df_spl_cluster, by="sample_ID")

df_long %>% arrange(spl_cluster) -> df_long
correct_sample_order <- unique(df_long$sample_ID)

ggplot(df_long,aes(x=taxon, y=sample_ID, fill=reads))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y = element_text(vjust = 0, hjust = 1, size = 6),
        axis.text.x = element_text(angle = 90, size = 6))+
  scale_fill_gradient(low="white", high="Black")+
  scale_x_discrete(limits = correct_taxon_order, label =
                     my_tax$class[match(my_tax$OTU, colnames(FB_bin))])+
  scale_y_discrete(limits = correct_sample_order)+
  geom_hline(yintercept=0.5+cumsum(table(LBM_bernoulli$memberships$Samples))[-length(table(LBM_bernoulli$memberships$Samples))],
             color = "red")+
  # geom_vline(xintercept=0.5+cumsum(table(LBM_bernoulli$memberships$ASV))[-length(table(LBM_bernoulli$memberships$ASV))],
             # color = "red")+
  ylab("Echantillons")+
  xlab("Bactéries")

# --------- Null model ------
library(parallel)

set.seed(1)
nbCores <- 1
n_iterations <- 10000

#Simulated binary networks----
NECH <- nrow(FB_bin)
nrow_FB_bin = nrow(FB_bin)

FB_bin %>%
  as.data.frame() %>%
  nullmodel(.,"curvball") %>%
  simulate(nsim =1) -> netBin.config #Fred préfère

#produce NPERM simulated graphs
#.Random.seed

simulate_membership <- function(OTU_matrix_t) {
  sbm::estimateBipartiteSBM(netMat = OTU_matrix_t,
                            model = "bernoulli",
                            estimOptions = list(
                              verbosity = 0,
                              plot = FALSE,
                              nbCores = nbCores)) -> LBM_sim

  LBM_sim$memberships$row
}

randomize_and_simulate <- function(n) {
  #argument n non utilisé mais obligatoire pour mclapply
    netBin.config <- simulate(nullmodel(as.data.frame(FB_bin), "curveball"), nsim = 1)
    simulate_membership(netBin.config[,,1])
}

system.time(
  parallel::mclapply(1:n_iterations,
                     randomize_and_simulate,
                     mc.cores = (parallel::detectCores()/2-1)) %>%
    unlist(.) %>%
    matrix(., nrow = nrow_FB_bin , byrow = TRUE) -> clust_bin
)
rm(sum_reads, sum_reads_per_OTU, sum_reads_per_samples)

save(clust_bin, file = "results_fruit/clust_NPERM10000_050623_bin.Rdata.Rdata")

load("results_fruit/clust_NPERM10000_050623_bin.Rdata")
#----NMI----
# Population
my_mtd_env$loc_code = sapply(1:nrow(my_mtd_env), function(i) {
  fruit_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][1]
  site_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][2]
  paste(fruit_code, site_code, sep = "-")
})
my_mtd_env$loc_code %>%
  str_replace("C-", "GF-") %>% str_replace("G-", "GB-") ->
  my_mtd_env$loc_code

my_mtd_env %>%
  filter(Label %in% rownames(Log_FB)) %>%
  pull("loc_code") %>% 
  as.factor() -> vect

memb_spl_obs_ber <- LBM_bernoulli$memberships$Samples

nmi_obs <- igraph::compare(as.factor(memb_spl_obs_ber), vect, method = "nmi")

clust_bin %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p1 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Population")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

#Température
my_mtd_env %>%
  filter(Label %in% rownames(FB_bin)) %>%
  pull("grp_t") %>% 
  as.factor() -> vect

nmi_obs <- compare(as.factor(memb_spl_obs_ber), vect, method = "nmi")

clust_bin %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p2 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Température")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

#Pluvio
my_mtd_env %>%
  filter(Label %in% rownames(FB_bin)) %>%
  pull("grp_plu") %>% 
  as.factor() -> vect

nmi_obs <- compare(as.factor(memb_spl_obs_ber), vect, method = "nmi")

clust_bin %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p3 <-ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Pluviométrie")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

#Env
my_mtd_env %>%
  filter(Label %in% rownames(FB_bin)) %>%
  pull("env") %>% 
  as.factor() -> vect

nmi_obs <- compare(as.factor(memb_spl_obs_ber), vect, method = "nmi")

clust_bin %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p4 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Fruit hôte")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

#sex
my_mtd_env %>%
  filter(Label %in% rownames(FB_bin)) %>%
  pull("sex") %>% 
  as.factor() -> vect

nmi_obs <- compare(as.factor(memb_spl_obs_ber), vect, method = "nmi")

clust_bin %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Simulation", values_to = "Memberships") %>%
  group_by(Simulation) %>% 
  summarise(., nmi_sim = compare(as.factor(Memberships), vect, method = "nmi")) -> simulated_nmi

p5 <- ggplot(simulated_nmi, aes(x = nmi_sim)) +
  geom_histogram(binwidth = 0.003, fill = "lightblue", color = "black") +
  geom_vline(xintercept = nmi_obs,
             color = "red",
             size = 1.2) + 
  labs(x = "valeur de congruence", y = "Effectif cumulé", title = "Sexe")

print(c("NMI =",nmi_obs, "P-value = ", 1 - ecdf(simulated_nmi$nmi_sim)(nmi_obs)))

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

library(cowplot)
plot_grid(p1, p5, blankPlot ,p2, p3, p4, labels=c("A", "B"," ","C","D","E"), ncol = 3, nrow = 2)

#Simulated binary networks----
NECH <- nrow(FB_bin)

#Alluvial plot----
library(alluvial)
library(dplyr)

env_vect = my_mtd_env$env[match(rownames(FB_bin), my_mtd_env$Label)]
sex_vect <- my_mtd_env$sex[match(rownames(FB_bin), my_mtd_env$Label)]
loc_vect <- my_mtd_env$locality[match(rownames(FB_bin), my_mtd_env$Label)]
T_vect <- my_mtd_env$grp_t[match(rownames(FB_bin), my_mtd_env$Label)] %>%
  recode("low" = "Basse", "mid" = "Moyenne", "high" = "Élevée")
Plu_vect <- my_mtd_env$grp_plu[match(rownames(FB_bin), my_mtd_env$Label)] %>%
  recode("1" = "Faible", "2" = "Moyenne", "3" = "Forte")
Alt_vect <- my_mtd_env$grp_alt[match(rownames(FB_bin), my_mtd_env$Label)]

my_mtd_env$loc_code = sapply(1:nrow(my_mtd_env), function(i) {
  fruit_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][1]
  site_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][2]
  paste(fruit_code, site_code, sep = "-")
})
my_mtd_env$loc_code %>% 
  str_replace("C-", "GF-") %>% str_replace("G-", "GB-") -> my_mtd_env$loc_code

loc_vect_bis <- my_mtd_env$loc_code[match(rownames(FB_bin), my_mtd_env$Label)]

B <-
  as.data.frame(table(env_vect,loc_vect_bis, as.factor(memb_spl_obs),sex_vect))
colnames(B) = c( "Fruit","Locality", "Gut bacteria",  "Sex", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
alluvial(B[, c(1, 2, 3, 4)], freq = B$Freq)

alluvial(B[, c(1, 2, 3, 4)],
         freq = B$Freq,
         col = case_when(B$Fruit == "Jamrosat" ~ "lightyellow",
                         B$Fruit == "Badamier" ~ "purple",
                         B$Fruit == "Goyavier fraise" ~ "red",
                         B$Fruit == "Goyavier blanc" ~ "lightgreen",
                         B$Fruit == "Mangue" ~ "orange",
                         TRUE ~ "blue"))

#Alluvial plot avec juste les métadonnées
B <- as.data.frame(table(env_vect,Plu_vect,T_vect,loc_vect))
colnames(B) = c( "Fruit","Pluviosité","Température","Localité", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]

alluvial(B[, c(1,2, 3,4)],
         freq = B$Freq,
         col = case_when(B$Fruit == "Jamrosat" ~ "yellow",
                         B$Fruit == "Badamier" ~ "purple",
                         B$Fruit == "Goyavier fraise" ~ "red",
                         B$Fruit == "Goyavier blanc" ~ "lightgreen",
                         B$Fruit == "Mangue" ~ "orange",
                         TRUE ~ "blue"))


#CCA (env, temp, pluv)----
library(lme4)

env_vect <- my_mtd_env$env[match(rownames(Log_FB), my_mtd_env$Label)]
temp_vect <- my_mtd_env$grp_t[match(rownames(Log_FB), my_mtd_env$Label)]
pluv_vect <- my_mtd_env$grp_plu[match(rownames(Log_FB), my_mtd_env$Label)]

cca.data <-
  data.frame(
    "env" = as.factor(env_vect),
    "temp" = as.factor(temp_vect),
    "pluv" = as.factor(pluv_vect)
  )

#Get the corresponding null distribution of clustering
clust_bin -> cle.config

#CCA with line randomization
wholenetwork.cca.results.NULL <- analysis.function.alt(as.data.frame(FB_bin),
                                                       memb_spl_obs_ber,
                                                       cca.data$env,
                                                       cca.data$temp,
                                                       cca.data$pluv,
                                                       configs = NULL)
#CCA with full randomization
wholenetwork.cca.results.CONFIGS <- analysis.function.alt(as.data.frame(FB_bin),
                                                          memb_spl_obs_ber,
                                                          cca.data$env,
                                                          cca.data$temp,
                                                          cca.data$pluv,
                                                          configs = cle.config)
write.csv2(wholenetwork.cca.results.CONFIGS, file = "results_fruit/cca_bin_TPE2.csv")

#Venn bianaire (CCA TPE)
memb_spl_obs_ber <- LBM_bernoulli$memberships$Samples

dummyMod<-dummy(as.factor(memb_spl_obs_ber), levelsToKeep =
                  levels(as.factor(memb_spl_obs_ber)))

varp <- varpart(dummyMod, cca.data$env, cca.data$pluv, cca.data$temp,  chisquare = T)

plot(varp, Xnames = list("Fruit hôte", "Pluviosité","Température"))

write.csv2(
  rbind(varp$part$fract, varp$part$indfract, varp$part$contr1),
  "results_fruit/SVD_varpart_result3.csv")

#CCA sur loc et sexe----

# Factors
sex_vect <- my_mtd_env$sex[match(rownames(FB_bin), my_mtd_env$Label)]

my_mtd_env$loc_code = sapply(1:nrow(my_mtd_env), function(i) {
  fruit_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][1]
  site_code = str_split(my_mtd_env$code_bennath[i], pattern="-")[[1]][2]
  paste(fruit_code, site_code, sep = "-")
})
my_mtd_env$loc_code %>%
  str_replace("C-", "GF-") %>% str_replace("G-", "GB-") ->
  my_mtd_env$loc_code

loc_vect_bis <- my_mtd_env$loc_code[match(rownames(FB_bin),
                                          my_mtd_env$Label)]

dummyMod<-dummy(as.factor(memb_spl_obs_ber), levelsToKeep =
                  levels(as.factor(memb_spl_obs_ber)))

#CCA
cca.data <-
  data.frame(
    "sex" = as.factor(sex_vect),
    "loc" = as.factor(loc_vect_bis)
  )

#CCA with line randomization
wholenetwork.cca.results.NULL <-
  analysis.function.alt2(as.data.frame(FB_bin),
                         memb_spl_obs_ber,
                         cca.data$sex,
                         cca.data$loc,
                         configs = NULL)

#CCA with full randomization
wholenetwork.cca.results.CONFIGS <-
  analysis.function.alt2(as.data.frame(FB_bin), 
                         memb_spl_obs_ber,
                         cca.data$sex,
                         cca.data$loc, 
                         configs = clust_bin)

write.csv2(wholenetwork.cca.results.CONFIGS, file = "results_fruit/cca_bin_loca_sexe_test.csv")

#Venn (CCA)
varp <- varpart(dummyMod, cca.data$sex, cca.data$loc, chisquare = T)

plot(varp, Xnames = list("Sexe", "Population"))

###SVD------
#install.packages("ROCR")
library(ROCR)
library(lme4)
library(sbm)
source('scriptENV/scriptsFRUIT/fruit 2023/functions_modified.R') # Kindly provided by Massol et al. 2021 JAE
library(igraph)
NPERM = 10000 # Number of configuration matrices for each rarefied matrix

env_vect = my_mtd_env$env[match(rownames(FB_bin), my_mtd_env$Label)]
sex_vect <- my_mtd_env$sex[match(rownames(FB_bin), my_mtd_env$Label)]
loc_vect <- my_mtd_env$locality[match(rownames(FB_bin), my_mtd_env$Label)]
loc_vect_bis <- my_mtd_env$loc_code[match(rownames(FB_bin), my_mtd_env$Label)]
T_vect <- my_mtd_env$grp_t[match(rownames(FB_bin), my_mtd_env$Label)]
Plu_vect <- my_mtd_env$grp_plu[match(rownames(FB_bin), my_mtd_env$Label)]
Alt_vect <- my_mtd_env$grp_alt[match(rownames(FB_bin), my_mtd_env$Label)]

NECH = nrow(FB_bin)

#Decomposition
netBin.svd <- svd(FB_bin)
netBin.U <- netBin.svd$u
netBin.S <- diag(netBin.svd$d)
netBin.Ssqrt <-
  structure(vapply(netBin.S, sqrt, numeric(1)), dim = dim(netBin.S))
netBin.traits <- netBin.U %*% netBin.Ssqrt
netBin.traits2 <- netBin.svd$v %*% netBin.Ssqrt

#choose the right number of vectors/traits
nmi_vect <- sapply(1:NECH, function(x) {
  approx_matrix <-
    approx.from.svd(FB_bin, netBin.traits, netBin.traits2, x)
  LBM_bernoulli_sim <- as.matrix(approx_matrix) %>%
    estimateBipartiteSBM(model = 'bernoulli', dimLabels = c(row = "Samples", col = "ASV"))
  memb_spl_sim <- LBM_bernoulli_sim$memberships$Samples
  compare(memb_spl_obs_ber, memb_spl_sim, method = "nmi")
})

plot(1:NECH,nmi_vect,
     ylim = c(0, 1),
     type = "l",
     xlab = "Nombres de vecteurs retenus",
     ylab = "Valeurs de congruence (NMI)"
)

load("results_fruit/n_retained.Rdata")

#RDA sur loc sexe----
NECH = nrow(FB_bin)

#RDA
dummySex <- dummy(sex_vect, levelsToKeep = levels(as.factor(sex_vect)))
dummyPop <- dummy(loc_vect_bis, levelsToKeep =
                    levels(as.factor(loc_vect_bis)))

N_RETAINED = 11

varp <- varpart(netBin.traits[, 1:N_RETAINED],dummySex, dummyLoc)

plot(varp, Xnames = list("Sexe", "Population"))

varp

write.csv2(
  rbind(varp$part$fract, varp$part$indfract, varp$part$contr1),
  "results_fruit/SVD_varpart_sexe_pop_result2.csv")

####tests with randomization of rows
dummyVar1 <- dummySex
dummyVar2 <- dummyPop

NPERM = 10000

tests.whole <- c(
  
  anova(rda(netBin.traits[, 1:N_RETAINED] ~ dummyVar1, na.action = 
              na.omit),  permutations = how(nperm = NPERM - 1))$P[1],
  anova(rda(netBin.traits[, 1:N_RETAINED] ~ dummyVar2, na.action = 
              na.omit), permutations = how(nperm = NPERM - 1))$P[1],
  anova(rda(netBin.traits[, 1:N_RETAINED] ~ dummyVar2 + dummyVar1, 
            na.action = na.omit), permutations = how(nperm = NPERM - 1))$P[1],
  anova(rda(netBin.traits[, 1:N_RETAINED] ~ dummyVar1 + 
              Condition(dummyVar2),na.action = na.omit), permutations = how(nperm = 
                                                                              NPERM - 1))$P[1],
  anova(rda(netBin.traits[, 1:N_RETAINED] ~ dummyVar2 + 
              Condition(dummyVar1), na.action = na.omit), permutations = how(nperm = 
                                                                               NPERM - 1))$P[1]
)
write.csv2(tests.whole,"results_fruit/wholenet_row_sexe_pop.csv")

####tests with randomization of edges
.Random.seed
netBin.config <-
  simulate(nullmodel(as.data.frame(FB_bin), "curveball"), nsim = NPERM)

varp.config <- sapply(1:NPERM,function(x)
  analysis.function.rdpg2(netBin.config[,,x], dummyVar1, dummyVar2,
                          N_RETAINED))

varp.config.ecdf <- apply(varp.config,1,ecdf)

varp.Rsq <- analysis.function.rdpg2(FB_bin, dummyVar1, dummyVar2,
                                    N_RETAINED)

tests2.whole <- sapply(1:7,function(x)
  1-((varp.config.ecdf[[x]])(varp.Rsq[x])))#all 7 p-values (some are not needed)

write.csv2(tests2.whole,"results_fruit/wholenet_edgeperm_sexe_loca.csv")

#RDA EPT----
load("results_fruit/n_retained.Rdata")

dummyEnv <- dummy(env_vect, levelsToKeep = levels(as.factor(env_vect)))
dummyT <- dummy(T_vect, levelsToKeep = levels(as.factor(T_vect)))
dummyPlu <- dummy(Plu_vect, levelsToKeep = levels(as.factor(Plu_vect)))


#Venn RDA
N_RETAINED = 11

varp <-
  varpart(netBin.traits[, 1:N_RETAINED], dummyEnv, dummyPlu, dummyT)

varp 

plot(varp, Xnames = list("Fruit hôte", "Pluviosité", "Température"))
write.csv2(
  rbind(varp$part$fract, varp$part$indfract, varp$part$contr1),
  "results_fruit/Venn_RDA_EPT.csv")

####tests with randomization of rows
dummyType <- dummyEnv
dummyDiet <- dummyPlu
dummyGenotype <- dummyT

NPERM = 10000

tests.whole <- c(
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyType, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyDiet, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyGenotype, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyType + dummyDiet, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyType + dummyGenotype, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyGenotype + dummyDiet, na.action = na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyType + dummyDiet + dummyGenotype,
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyType + Condition(dummyDiet) + Condition(dummyGenotype),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyDiet + Condition(dummyGenotype) + Condition(dummyType),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyGenotype + Condition(dummyType) + Condition(dummyDiet),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyType + Condition(dummyGenotype),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyType + Condition(dummyDiet), na.action =
          na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyDiet + Condition(dummyGenotype),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(netBin.traits[, 1:N_RETAINED] ~ dummyDiet + Condition(dummyType), na.action =
          na.omit),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyGenotype + Condition(dummyType),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1],
  anova(
    rda(
      netBin.traits[, 1:N_RETAINED] ~ dummyGenotype + Condition(dummyDiet),
      na.action = na.omit
    ),
    permutations = how(nperm = NPERM - 1)
  )$P[1]
)
write.csv2(tests.whole, file = "results_fruit/RDA_rowperm2_TPE.csv")

####tests with randomization of edges
.Random.seed
netBin.config <-
  simulate(nullmodel(as.data.frame(FB_bin), "curveball"), nsim = NPERM)

varp.config<-sapply(1:NPERM,function(x) analysis.function.rdpg(netBin.config[,,x],dummyType,dummyDiet,dummyGenotype,N_RETAINED))
varp.config.ecdf<-apply(varp.config,1,ecdf)
varp.Rsq<-analysis.function.rdpg(FB_bin,dummyType,dummyDiet,dummyGenotype,N_RETAINED)
tests2.whole<-sapply(1:21,function(x) 1-((varp.config.ecdf[[x]])(varp.Rsq[x])))#all 21 p-values (some are not needed)
write.csv2(tests2.whole, file = "results_fruit/RDA_edgeperm2_.csv")

# Test RDA à la main----
library(vegan)

df1 <- as.data.frame(t(otu_table(ps_rare))) %>% 
  decostand(., method = "hellinger")
my_mtd_env %>%
  filter(Label != c("09_B-ETA-male","24_G-VER-fem")) -> mtd_perm

my_rda <- rda(df1 ~ env + grp_plu + grp_t, data = mtd_perm)
summary(my_rda)
my_rda
plot(my_rda)
RsquareAdj(my_rda)
anova.cca(my_rda, permutations = 999)
anova.cca(my_rda, by = "axis")
anova.cca(my_rda, by = "terms")
sqrt(vif.cca(my_rda))

install.packages("packfor", repos = "http://R-Forge.R-project.org")

test <- cca(df)
summary(test, display = "reg") 
plot(test)

data(dune)
cca(dune)

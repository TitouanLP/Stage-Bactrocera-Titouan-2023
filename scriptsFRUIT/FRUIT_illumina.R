#*******************#
#### NGB - FRUIT ####
#*******************#

#Packages
library(tidyverse)
library(phyloseq)

#Load data----
#**Nouvelle taxo----
silva_data = "databioinfo/tax_slv_ssu_138.1.txt.gz"
illumina_data = "dataFRUIT/NGB_16S_Illumina_214_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table3"

## ------------------------------------------------------------------ functions

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

## ----------------------------------------------------------------------- main

## find potential issues
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

#Basic statistics----
#Nombre de séquences : 1167
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%  dim()

#Nombre de séquences non assignées à quoi que ce soit : 80
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  filter(is.na(references)) %>% dim()

#Nombre de séquences avec une assignation multiple : 690 (622 assignations distinctes)
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  filter(!is.na(references)) %>% 
  filter(str_detect(references, ",")) %>% 
  select(references) %>% duplicated() %>% table()


#*******************----
#Building a phyloseq object for Illumina data (based on OTUs, not on taxonomy)----
#toutes les données dans même objet
#**OTU table----
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  filter(!is.na(references)) %>%
  select(-taxonomy,-total, -cloud,-amplicon,-length,-abundance,-chimera,-spread,-quality,-sequence,-identity,-references) %>% 
  replace(is.na(.),0) %>% 
  column_to_rownames(var = "OTU") -> my_otu_table

my_otu_table %>% colnames()

#c'est ici que tu peux essayer de mettre des noms d'échantillons plus sensés (avec match et le fichier MTD)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%
  head(n=2) %>% 
  pivot_longer(cols = -c(1:13), names_to = "sample_code", values_to = "reads") %>% 
  filter(reads > 0) %>%
  separate_rows(references, sep = ",") %>% 
  group_by(OTU, sample_code) %>% 
  mutate(corrected_reads = reads/n()) %>% 
  ungroup() %>% 
  left_join(read_tsv(file = "databioinfo/SILVA_138.1_SSURef_tax_silva.acc2taxo.gz"), by=c("references" = "accession")) %>% 
  group_by(OTU, total, amplicon, length, abundance, quality, spread, chimera, sequence, identity, sample_code,taxonomy.y, reads, corrected_reads) %>%
  summarize(sum_corrected_reads = sum(corrected_reads), .groups = "drop") -> tmp_sum

# Fonction pour match (merci GPT4)
label_by_sample <- my_mtd %>%
  select(sample_ID, Label) %>%
  distinct()

# Remplacer les "sample_code" par les "label" correspondants (merci GPT4)
tmp_sum2 <- tmp_sum %>%
  left_join(label_by_sample, by = c("sample_code" = "sample_ID")) %>%
  mutate(sample_code = if_else(!is.na(Label), as.character(Label), NA_character_)) %>%
  filter(!is.na(sample_code)) %>% 
  select(-Label)


write.table(tmp_sum2, file = paste0("C:/Users/titou/OneDrive/Bureau/Stage/CIRAD Master 1/R stage/Stage-Bactrocera-Titouan-2023/dataFRUIT/tmp_sum2.csv"), row.names = FALSE)

  # write.table(tmp, file = paste0("C:/Users/titou/OneDrive/Bureau/Stage/CIRAD Master 1/R stage/Stage-Bactrocera-Titouan-2023/dataFRUIT/tmp.csv"), row.names = FALSE)

  # library(dplyr)
  # read_tsv(illumina_data,
  #          na = c("", "NA", "No_hit", "0", "0.0"),
  #          show_col_types = FALSE) %>%
  #   head(n=2) %>%
  #   mutate(total_reads = rowSums(select(., -(1:13)), na.rm = TRUE)) %>%
  #   select(total_reads) %>%
  #   view()
  # script GPT4 pour vérifier si la somme des reads de chaque échantillon par OTU est bien égal au total  

#**TAX table----
my_ref_list <- rownames(my_otu_table)

#***Différentes approches pour gérer les assignations multiples----
#Approche "Ancêtre commun"
load_assignments(illumina_data) %>%
  separate_rows(taxonomy, sep = "\\|") %>%
  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy" = "taxonomic_path")) %>%
  filter(taxonomic_rank %in% good_taxonomic_levels) %>%
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>%
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy) %>%
  column_to_rownames(var = "OTU") -> my_tax_table

#Approche "Distribution des reads"
  read_tsv(illumina_data,
           na = c("", "NA", "No_hit", "0", "0.0"),
           show_col_types = FALSE) %>% 
    select(OTU, taxonomy, references) %>% 
    separate_rows(references, sep = ",") %>% 
    #récupérer les chemins complets dans silva
    left_join(read_tsv(file = "databioinfo/SILVA_138.1_SSURef_tax_silva.acc2taxo.gz"), by=c("references" = "accession")) %>% 
    select(-taxonomy.x) %>%  
  separate_rows(taxonomy.y, sep = "\\|") %>%
  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy.y" = "taxonomic_path")) %>%
  filter(taxonomic_rank %in% good_taxonomic_levels) %>%
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>%
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy.y) %>% 
    mutate(code_unique = paste(OTU, references, sep="_")) %>% 
    select(-OTU, -references) %>% 
  column_to_rownames(var = "code_unique") -> my_tax_table

  colnames(tmp)
  
  length(unique(tmp$references))
  view(tmp)
  
#Il y a 4 OTUS assignées à autre chose que des bactéries  
my_tax_table %>% 
  filter(!domain=="Bacteria") %>% 
  dim()

my_otu_table %>% 
  filter(rownames(.) %in% c("495","553","813","1006")) %>% 
  apply(.,1,sum)

my_tax_table %>% 
  filter(domain=="Bacteria") %>% 
  rownames() -> OTU_to_keep

dim(my_otu_table)

head(my_otu_table)

dim(my_tax_table)

head(my_tax_table)

# filter(my_otu_table, rownames(my_otu_table) %in% rownames(my_tax_table)) -> my_otu_table

otu_table(my_otu_table,taxa_are_rows = TRUE) -> OTU

tax_table(as.matrix(my_tax_table)) -> TAX

#**Metadata----
my_mtd <- read.csv2("dataFRUIT/metadata_illumina.csv") %>% 
  subset(., experiment == "RUN4")

rownames(my_mtd) <- my_mtd$sample_ID

summary(my_mtd)

sample_data(my_mtd)->MTD

#**Phyloseq object----
ps_illumina = phyloseq(OTU, TAX, MTD)

ps_illumina = prune_taxa(OTU_to_keep, ps_illumina)

#**Add read depth and percentage of unassigned OTUs----
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  select(-taxonomy,-total, -cloud,-amplicon,-length,-abundance,-chimera,-spread,-quality,-sequence,-identity,-references) %>% 
  replace(is.na(.),0) %>% 
  column_to_rownames(var = "OTU") %>% 
  apply(., 2,sum) -> total_reads

total_reads

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  filter(is.na(references)) %>% 
  select(-taxonomy,-total, -cloud,-amplicon,-length,-abundance,-chimera,-spread,-quality,-sequence,-identity,-references) %>% 
  replace(is.na(.),0) %>% 
  column_to_rownames(var = "OTU") %>% 
  apply(., 2,sum) -> unassigned_reads

unassigned_reads

sample_data(ps_illumina)$total_reads <- total_reads

sample_data(ps_illumina)$unassigned_reads <- unassigned_reads ## A quoi servent ces deux lignes ? 

sample_data(ps_illumina)$unassigned_percentages <- 
  100*sample_data(ps_illumina)$unassigned_reads/sample_data(ps_illumina)$total_reads

sample_data(ps_illumina)$bacterial_reads <- sample_sums(ps_illumina)

# sample_data(ps_illumina) %>% View()

#Add an alert on very low sequencing depths
sample_data(ps_illumina)$failed<-(sample_data(ps_illumina)$total_reads<5000 & sample_data(ps_illumina)$type=="ENV") 

#******************####
#Distribution of unassigned reads----
#**RUN2 and RUN4----
mtd<-data.frame(sample_data(ps_illumina))
mtd<- mtd %>% arrange(unassigned_percentages)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = unassigned_percentages), color=env, shape=experiment)+
  geom_point()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=5))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 et RUN4 - Illumina",
       x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2&4_unassigned_ALL_sorted.pdf")

mtd<- mtd %>% arrange(experiment, type, Label)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = unassigned_percentages, color=env, shape=experiment))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=5))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 et RUN4 - Illumina",
       x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2&4_unassigned_ALL.pdf")

#**RUN4----
ps <- prune_samples(sample_data(ps_illumina)$experiment=="RUN4", ps_illumina)

mtd<-data.frame(sample_data(ps))
mtd<- mtd %>% arrange(env, locality)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = unassigned_percentages, color=env))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN4_unassigned_ALL.pdf")

ggplot(mtd, aes(x = env, y = unassigned_percentages, fill=env))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN4_unassigned_BY_FRUIT.pdf")

ggplot(mtd, aes(x = locality, y = unassigned_percentages,fill=locality))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN4_unassigned_BY_SITE.pdf")

ggplot(mtd, aes(x = sex, y = unassigned_percentages, fill=sex))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN4_unassigned_BY_SEX.pdf")

#**RUN2----
ps <- prune_samples(sample_data(ps_illumina)$experiment=="RUN2", ps_illumina)

mtd<-data.frame(sample_data(ps))
mtd<- mtd %>% arrange(env, host, env)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

summary(mtd)

ggplot(mtd, aes(x = sample_ID, y = unassigned_percentages, color=env))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=6))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 - Illumina",
       x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2_unassigned_ALL.pdf")

ggplot(mtd, aes(x = env, y = unassigned_percentages, fill=env))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2_unassigned_BY_FRUIT.pdf")

ggplot(mtd, aes(x = host, y = unassigned_percentages,fill=host))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2_unassigned_BY_HOST.pdf")

ggplot(mtd, aes(x = sex, y = unassigned_percentages, fill=sex))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=0))+  
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Percentage of unassigned reads")

ggsave("res/RUN2_unassigned_BY_SEX.pdf")

#******************####

#Distribution of the number of reads across samples----
#**RUN2 & RUN4----
mtd<-data.frame(sample_data(ps_illumina))
mtd<- mtd %>% arrange(total_reads)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = total_reads,color=type, shape=experiment))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=6))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 et RUN4 - Illumina",
       x ="Samples", y = "Sequencing depth (total read number)")

ggsave("res/RUN2&4_sequencing_depth_sorted_ALL.pdf")

mtd<- mtd %>% arrange(experiment, type, Label)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = total_reads,color=env, shape=experiment))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4, size=5))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 et RUN4 - Illumina",
       x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN2&4_sequencing_depth_ALL.pdf")

#**RUN4----
ps <- prune_samples(sample_data(ps_illumina)$experiment=="RUN4", ps_illumina)

mtd<-data.frame(sample_data(ps))
mtd<- mtd %>% arrange(total_reads)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x = sample_ID, y = total_reads,color=env))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN4_sequencing_depth_ALL.pdf")

#3 échantillons ont de très faibles profondeurs de séquençage = du même ordre que les témoins négatifs

mtd %>% 
  filter(total_reads<25000) %>% 
  filter(type=="ENV") -> low_read_samples


ggplot(mtd, aes(x = env, y = total_reads, fill=env))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN4_sequencing_depth_BY_FRUIT.pdf")

ggplot(mtd, aes(x = locality, y = total_reads,fill=locality))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN4_sequencing_depth_BY_SITE.pdf")

ggplot(mtd, aes(x = sex, y = total_reads, fill=sex))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+  
  labs(title="Effet fruit - RUN4 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN4_sequencing_depth_BY_SEX.pdf")

#**RUN2----
ps <- prune_samples(sample_data(ps_illumina)$experiment=="RUN2", ps_illumina)

mtd<-data.frame(sample_data(ps))
mtd<- mtd %>% arrange(env, host, env)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

summary(mtd)

ggplot(mtd, aes(x = sample_ID, y = total_reads, color=env))+
  geom_point()+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  labs(title="Effet fruit - RUN2 - Illumina",
       x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN2_sequencing_depth_ALL.pdf")

ggplot(mtd, aes(x = env, y = total_reads, fill=env))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN2_sequencing_depth_BY_FRUIT.pdf")

ggplot(mtd, aes(x = host, y = total_reads, fill=host))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN2_sequencing_depth_BY_HOST.pdf")

ggplot(mtd, aes(x = sex, y = total_reads, fill=sex))+
  geom_boxplot()+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=0))+  
  labs(title="Effet fruit - RUN2 - Illumina", x ="Samples", y = "Sequencing depth (read number)")

ggsave("res/RUN2_sequencing_depth_BY_SEX.pdf")

#******************####
#Decontamination----

ps <- prune_samples(sample_data(ps_illumina)$experiment=="RUN4", ps_illumina)

mtd<-data.frame(sample_data(ps))

#Distribution of bacterial reads
mtd<- mtd %>% arrange(env, host, env)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(log10(bacterial_reads), fill=type))+
  geom_histogram()+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Bacterial reads", y = "Frequency")

ggsave("res/RUN4_Distribution_of_bacterial_reads.pdf")

#Distribution of total reads
mtd<- mtd %>% arrange(total_reads)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

ggplot(mtd, aes(x=sample_ID, y=total_reads, col=type))+
  geom_point(size=3)+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)

ggsave("res/RUN4_Distribution_of_total_reads.pdf")

#Barplot
ps_even=transform_sample_counts(ps,  function(x) 1E6 * x / sum(x) )
plot_bar(physeq = ps_even, fill="Order")+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))
#Controls and environmental samples have clearly different profiles. Empty has a profile close to environmental samples. I will not count it as a control in the following. The "failed" samples (ie, those with very low sequencing depths = 67_S43 and 78_S54) have normal profiles.

#**Decontam----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
install.packages("decontam")
library("decontam")

#DNA concentration method
contamdf.freq <- isContaminant(ps, method="frequency", conc=sample_data(ps)$DNA_quantif)
table(contamdf.freq$contaminant)
#42 contaminants selon la méthode freq

ps_contam1 <- prune_taxa(contamdf.freq$contaminant, ps)
list_decontam_freq=taxa_names(ps_contam1)

#Prevalence method
contamdf.prev <- isContaminant(ps, method="prevalence", neg=(sample_data(ps)$type=="NEG"))
table(contamdf.prev$contaminant)

ps_contam2 <- prune_taxa(contamdf.prev$contaminant, ps)
list_decontam_prev=taxa_names(ps_contam2)

intersect(taxa_names(ps_contam1),taxa_names(ps_contam2))

#Graph illustrating the DNA concentration method
#Contaminant OTUs are expected to have less reads in higher concentration wells 
plot_frequency(ps, list_decontam_freq, conc="DNA_quantif") + 
  xlab("DNA Concentration")
#Non-contaminant OTUs should not have less reads in higher concentration wells 
list_ok_freq=sample(taxa_names(ps)[which(!(taxa_names(ps)%in%list_decontam_freq))],12)
plot_frequency(ps, list_ok_freq, conc="DNA_quantif") + 
  xlab("DNA Concentration")

#Distribution of contaminant and non-contaminant OTUs in environmental samples and controls
#Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "NEG", ps.pa)
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$type == "NEG"), ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.freq$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position="jitter") +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("Prevalence of conserved and removed OTUs following decontam - DNA concentration method")

df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("Prevalence of conserved and removed OTUs following decontam - prevalence method")

#**microDecon----
# install.packages("remotes")
# remotes::install_github("donaldtmcknight/microDecon")
library("microDecon")

df_decon<-data.frame(OTU_ID=taxa_names(ps))
ps_neg<-prune_samples(sample_data(ps)$host %in% c("Library_PCR", "Lysis_buffer","Water"), ps)
df_decon<-cbind.data.frame(df_decon, otu_table(ps_neg))  
ps_env<-prune_samples(sample_data(ps)$type %in% c( "ZYMO"),ps)
df_decon<-cbind.data.frame(df_decon, otu_table(ps_env))

decontaminated <- decon(data = df_decon,numb.blanks=length(sample_data(ps_neg)$sample_ID),numb.ind=length(sample_data(ps_env)$sample_ID), taxa=F)

list_microDecon<-sort(decontaminated$OTUs.removed$OTU_ID)
list_microDecon %>% length()

#View(tax_table(prune_taxa(taxa_names(ps) %in% list_microDecon, ps)))

list_decontam_freq %>% length()
list_decontam_prev %>% length()

intersect(list_decontam_prev, list_microDecon) %>% length()
intersect(list_decontam_freq, list_microDecon) %>% length()
intersect(list_decontam_prev, intersect(list_decontam_freq, list_microDecon)) %>% length()

#Distribution of contaminant and non-contaminant OTUs in environmental samples and controls
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$type == "NEG", ps.pa)
ps.pa.pos <- prune_samples(!(sample_data(ps.pa)$type == "NEG"), ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=(taxa_names(ps.pa)%in%decontaminated$OTUs.removed$OTU_ID))
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)")  + 
  ggtitle("Prevalence of conserved and removed OTUs following microDecon")

#Decontamination per se
cleaned_otu_table<-decontaminated$decon.table
rownames(cleaned_otu_table)<-cleaned_otu_table$OTU_ID
cleaned_otu_table$OTU_ID<-NULL
cleaned_otu_table$Mean.blank<-NULL

ps_cleaned<-prune_samples(sample_data(ps)$sample_ID %in% colnames(cleaned_otu_table), ps)
ps_cleaned<-prune_taxa(taxa_names(ps_cleaned) %in% rownames(cleaned_otu_table), ps_cleaned)
ps_cleaned<-prune_taxa(taxa_sums(ps_cleaned)>0, ps_cleaned)

cleaned_otu_table %>%
  filter(rownames(cleaned_otu_table) %in% taxa_names(ps_cleaned)) -> cleaned_otu_table

otu_table(ps_cleaned)<-otu_table(cleaned_otu_table,taxa_are_rows = TRUE)

#**Mock----
ps_mock<-prune_samples(sample_data(ps)$type == "ZYMO", ps)
ps_mock<-prune_taxa(taxa_sums(ps_mock)>0, ps_mock)

ps_mock_cleaned<-prune_samples(sample_data(ps_cleaned)$type == "ZYMO", ps_cleaned)
ps_mock_cleaned<-prune_taxa(taxa_sums(ps_mock_cleaned)>0, ps_mock_cleaned)

#Expected distribution of species relative frequencies in Zymo mock
Mock.expected=read.csv2(file="data/mock_zymo_expected.csv")
Mock.expected$Rel_ab=as.numeric(Mock.expected$Rel_ab)

#proportion of assignments at each level 
apply(tax_table(ps_mock)[,2:7],2,function(x){1-mean(is.na(x))})

# tax_table(ps_mock) %>% View()
# tax_table(ps_mock_cleaned) %>% View()

table(as.data.frame(tax_table(ps_mock))$Genus %in% Mock.expected$Genus)
table(as.data.frame(tax_table(ps_mock_cleaned))$Genus %in% Mock.expected$Genus)

#xxxxxxxxxxxxx
#Any level
transform_ps_into_df=function(ps_object, level, na_choice){
  as.data.frame(tax_table(ps_object)) %>%  #get taxonomy
    rownames_to_column(var = "tax_ID") -> tmp 
  as.data.frame(otu_table(ps_object)) %>% #get reads 
    rownames_to_column(var = "tax_ID") %>% 
    left_join(tmp, ., by="tax_ID")  %>%  #assemble taxonomy and reads 
    replace(is.na(.), "Others") -> tmp
  if(na_choice) tmp else filter(tmp, !tmp[[level]]=="Others") -> tmp
  as.data.frame(tapply(tmp$Tpos_S78, tmp[[level]], sum)/sum(tmp$Tpos_S78)) %>%
    rownames_to_column(var = level) 
}

compare_decontamination_methods=function(level, na_choice){
  transform_ps_into_df(ps_mock, level, na_choice) %>% 
    rename("Uncleaned"=colnames(.)[2]) -> tab_mock
  
  transform_ps_into_df(ps_mock_cleaned, level, na_choice) %>% 
    rename("microDecon"=colnames(.)[2]) -> tab_mock_microDecon
  
  prune_taxa(!contamdf.prev$contaminant, ps) %>% 
    prune_samples(sample_data(.)$type=="ZYMO",.) %>% 
    prune_taxa(taxa_sums(.)>0, .) %>% 
    transform_ps_into_df(., level, na_choice) %>% 
    rename("Decontam_PREV"=colnames(.)[2]) -> tab_mock_decontam_PREV
  
  prune_taxa(!contamdf.freq$contaminant, ps) %>% 
    prune_samples(sample_data(.)$type=="ZYMO",.) %>% 
    prune_taxa(taxa_sums(.)>0, .) %>% 
    transform_ps_into_df(., level, na_choice) %>% 
    rename("Decontam_FREQ"=colnames(.)[2]) -> tab_mock_decontam_FREQ
  
  full_join(tab_mock,tab_mock_decontam_FREQ, by=level) %>% 
    full_join(.,tab_mock_decontam_PREV, by=level) %>% 
    full_join(.,tab_mock_microDecon, by=level)
}

Mock.expected %>% 
  select("Family", "Rel_ab") %>%
  group_by(Family) %>%
  mutate(Rel_ab=sum(Rel_ab)) %>%
  distinct() %>% 
  full_join(.,compare_decontamination_methods("Family", FALSE)) %>% 
  pivot_longer(cols = !Family, names_to="Treatment", values_to = "Abundance") %>%
  replace(is.na(.), 0) %>% 
  ggplot(aes(x = Treatment, y=Abundance, fill=Family))+
  geom_col()

#LOrsqu'on fait tourner microDecon ou decontam_PREV sur les échantillons environnementaux + zymo, la decontamination supprime 92 OTUs dont des OTUs importantes de la mock (ex. Bacillus). Si on fait tourner microDecon ou decontam_PREV uniquement sur ZYMO, la mock est bcp moins affectée. Seules 10 OTUs sont retirées totalement et aucune n'est un vrai positif de la mock. Cette grande variabilité de résultat selon le poll d'échantillons positifs utilisés est inquiétante par rapport à cette méthode de décontamination. 

plot_bar(ps_mock, fill="Genus")

#on ne peut pas descendre à l'espèce
plot_bar(ps_mock, fill="Species")

#

cond <- is.na(tax_table(ps_mock)[,"Phylum"])|(tax_table(ps_mock)[,"Phylum"] %in% Mock.expected$Phylum & is.na(tax_table(ps_mock)[,"Class"]))|(tax_table(ps_mock)[,"Class"] %in% Mock.expected$Class & is.na(tax_table(ps_mock)[,"Order"]))|(tax_table(ps_mock)[,"Order"] %in% Mock.expected$Order & is.na(tax_table(ps_mock)[,"Family"]))|(tax_table(ps_mock)[,"Family"] %in% Mock.expected$Family & is.na(tax_table(ps_mock)[,"Genus"]))|(tax_table(ps_mock)[,"Genus"] %in% Mock.expected$Genus)

ps_mock_not_ok=prune_taxa(taxa_names(ps_mock)[which(!cond, taxa_names(ps_mock))], ps_mock)

as.data.frame(otu_table(ps_mock_not_ok)) %>% View()

#l'OTU dont on est sur qu'elle est un faux négatif la plus abondante est à 169 reads sur un total de 138398, soit une abondance relative de 0.001221116 en incluant les unassigned

sample_sums(ps_mock)

tax_table(ps_mock_not_ok) %>% View()

ps_mock_ok=prune_taxa(tax_table(ps_mock)[,"Genus"] %in% Mock.expected$Genus, ps_mock)
tax_table(ps_mock_ok) %>% View()
taxa_sums(ps_mock_ok) 

#Il y a des tas d'OTUs à très peu de reads (1 à 5) et qui sont d'un genre présent dans la mock

#
plot_bar(ps_mock_ok, fill="Order")+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina")+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))


ps_mock_even=transform_sample_counts(ps_mock_ok,  function(x) 100*x / sum(x) )
otu_table(ps_mock_even)


#**Galan et al----
ps_neg<-prune_samples(sample_data(ps)$type=="NEG", ps)
ps_neg<-prune_taxa(taxa_sums(ps_neg)>0, ps_neg)
# il y a 177 OTUs présentes dans les témoins négatifs

tax_table(ps_neg) 


data.frame(otu_table(ps_neg)) %>% 
  apply(.,1,max) -> cross_contam_threshold

cross_contam_threshold %>% length()

data.frame(taxa_sums(ps)) %>%
  filter(rownames(.) %in%  taxa_names(ps_neg)) -> total_reads_cont

t(total_reads_cont) %>% length()

cross_contam_threshold/t(total_reads_cont) -> cross_contam_threshold2

plot(log10(cross_contam_threshold), cross_contam_threshold2)

#******************####
#Barplot----
ps_even=transform_sample_counts(ps,  function(x) 1E6 * x / sum(x) )
prune_samples(sample_data(ps)$failed == "TRUE", ps)
plot_bar(physeq = ps_even, fill="Family")+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))


mtd<-data.frame(sample_data(ps))
mtd<- mtd %>% arrange(env, locality)
sample_order<-mtd$sample_ID
label_order<-mtd$Label

sample_order=c("RUN1_barcode02","RUN1_barcode06","RUN1_barcode03","RUN1_barcode07","RUN1_barcode01", "RUN1_barcode05","RUN3_barcode01","RUN2C_barcode04","RUN2C_barcode05","RUN2C_barcode06","RUN3_barcode03","RUN3_barcode04","RUN1_barcode04","RUN3_barcode05","RUN2C_barcode07", "RUN2C_barcode08")

label_order=mtd$label[match(sample_order, mtd$sample_ID)]

ps_even=transform_sample_counts(physeq001,  function(x) 1E6 * x / sum(x) )

col<-paletteer_d("ggsci::default_igv")[c(1,31,33,37,3,45,7,4,43,2,32)]
plot_bar(physeq = ps_even, fill="Order")+
  scale_x_discrete(limits=sample_order, labels=label_order)+
  scale_fill_manual(values=col)+
  scale_y_continuous(name="Relative abundance",labels=c())+
  theme_minimal()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90, hjust=2.5))

ps = tax_glom(ps_cleaned, "Genus")
sample_data(ps)


save(ps, ps_illumina, taxonomy_dict, file="interm.Rdata")

ps2 <- prune_samples(sample_data(ps)$type=="ENV" & sample_data(ps)$experiment=="RUN4", ps)

plot_bar(physeq = ps2, x= "sample_ID", fill="Family")+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=0))

plot(sort(sample_sums(ps2)))

#Les témoins des 2 RUNS
ps <- prune_samples(!(sample_data(ps_illumina)$type=="ENV"), ps_illumina)

plot_bar(physeq = ps, x= "sample_ID", fill="Phylum")+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=0))


ps_even=transform_sample_counts(ps,  function(x) 1E6 * x / sum(x) )

ps4 <- prune_samples(sample_data(ps)$type=="ENV" & sample_data(ps)$experiment=="RUN4", ps)

ps_even=transform_sample_counts(ps4,  function(x) 1E6 * x / sum(x) )

my_plot_bar(physeq = ps_even, x= "sample_ID", fill="Phylum")+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=0))+
  scale_x_discrete(labels=sample_data(ps_even)$host)


mat<-ggplot(df_long,aes(x=label, y=taxon, fill=log10(reads+1)))+
  coord_flip()+
  geom_tile()+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1,face="italic"),
        axis.text.y = element_text())+
  scale_fill_gradient(low="white", high="Black")+
  scale_y_discrete(name="Bacterial genus", limits = correct_taxon_order)+
  scale_x_discrete(name="Sample", limits = correct_sample_order)+
  geom_vline(xintercept=0.5+cumsum(table(SBMobject$memberships$row))[-length(table(SBMobject$memberships$row))], color = "black", size=1)+
  geom_hline(yintercept=0.5+cumsum(table(SBMobject$memberships$col))[-length(table(SBMobject$memberships$col))], color = "black", size=1)
mat

#--------- LBM Clustering ------
LBM_log_reads_gaussian<- 
  as.matrix(Log_FB) %>% 
  estimateBipartiteSBM(
    model = 'gaussian')

memb_spl_obs<-LBM_log_reads_gaussian$memberships$row

LBM_log_reads_gaussian$storedModels %>% arrange(ICL)

myOrderedMatrixPlot=function(SBMobject, BF_mat){
  colnames(BF_mat)=mtd$label[match(colnames(BF_mat), mtd$sample_ID)]
  df=as.data.frame(BF_mat)
  df$taxon=rownames(df)
  df_long<-pivot_longer(df, cols=names(df)[-17], names_to = "label", values_to = "reads")
  df_long<-left_join(df_long,mtd, by="label")
  
  Qrow<-SBMobject$nbBlocks["row"][[1]]
  Qcol<-SBMobject$nbBlocks["col"][[1]]
  
  df_tax_cluster=data.frame(taxon=rownames(BF_mat),tax_cluster=SBMobject$memberships$col)
  df_long<-left_join(df_long, df_tax_cluster, by="taxon")
  
  df_long<-df_long %>% arrange(tax_cluster)
  correct_taxon_order <- unique(df_long$taxon)
  
  df_spl_cluster=data.frame(label=colnames(BF_mat)
                            ,spl_cluster=SBMobject$memberships$row)
  df_long<-left_join(df_long, df_spl_cluster, by="label")
  
  df_long<-df_long %>% arrange(spl_cluster)
  correct_sample_order <- unique(df_long$label)
  
  mat<-ggplot(df_long,aes(x=label, y=taxon, fill=log10(reads+1)))+
    coord_flip()+
    geom_tile()+
    theme_minimal()+
    theme(legend.position="top",
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1,face="italic"),
          axis.text.y = element_text())+
    scale_fill_gradient(low="white", high="Black")+
    scale_y_discrete(name="Bacterial genus", limits = correct_taxon_order)+
    scale_x_discrete(name="Sample", limits = correct_sample_order)+
    geom_vline(xintercept=0.5+cumsum(table(SBMobject$memberships$row))[-length(table(SBMobject$memberships$row))], color = "black", size=1)+
    geom_hline(yintercept=0.5+cumsum(table(SBMobject$memberships$col))[-length(table(SBMobject$memberships$col))], color = "black", size=1)
  mat
}
myOrderedMatrixPlot(LBM_log_reads_gaussian, BF)

ggsave("res/weigthed_incid_mat.svg")



# MOCKS ####
#Expected distribution of species relative frequencies in Zymo mock
# Mock.expected=read.csv2(file="data/mock_zymo_expected.csv")
# Mock.expected$Rel_ab=as.numeric(Mock.expected$Rel_ab)
# 
# #Loading mock samples data
# ps_mock<-subset_samples(ps, type=="ZYMO")
# 
# #Distribution of the number of reads across samples 
# sample_sums(ps_mock) #for each sample
# sum(sample_sums(ps_m)ock)) #total mock reads
# min(sample_sums(ps_mock)) #min number of reads
# max(sample_sums(ps_mock)) #max number of reads
# 
# #proportion of assignments we have made at each level 
# apply(tax_table(ps_mock)[,2:7],2,function(x){1-mean(is.na(x))})

#Working at  level####

#Sequencing depth
my_mtd$total_reads <- sapply(my_mtd$sample_ID, function(spl) sample_sums(ps_illumina)[[spl]])
names(my_mtd)
my_mtd %>% 
  filter(experiment=="RUN4") %>% 
  ggplot(my_mtd, mapping=aes(x=sample_ID, y=total_reads, color=type)) +
  geom_point()+
  theme_minimal()+
  theme(
    legend.position="right",
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(vjust = .5, hjust=1,angle =90))

my_mtd %>% 
  filter(experiment=="RUN4") %>%
  filter(total_reads<25000) %>% View()

my_mtd %>% 
  filter(experiment=="RUN4") %>%
  filter(type%in%c("NEG", "ZYMO")) %>% 
  View()

ps_neg<-subset_samples(ps_illumina, type%in%c("NEG", "ZYMO"))

ps = tax_glom(ps_neg, "Genus")

ps_even=transform_sample_counts(ps,  function(x) 1E6 * x / sum(x) )

my_plot_bar(ps_even, fill="Family")+
  theme(legend.position="none")

ps_mock<-subset_samples(ps_illumina, type=="ZYMO")

ps_mock_GENUS = tax_glom(ps_mock, fill="Genus")

#Observing the number of reads of false positive and false negative taxa
#Setting as dataframe
df_taxo=as.data.frame(tax_table(ps_mock_GENUS)) #get taxonomy
df_taxo$tax_ID=rownames(df_taxo)
df_reads=as.data.frame(otu_table(ps_mock_GENUS)) #get reads
df_reads$tax_ID=rownames(df_reads)
df_wide=left_join(df_taxo, df_reads, by="tax_ID") #assemble taxonomy and reads
df_mock=pivot_longer(df_wide, cols = names(df_wide)[9:length(names(df_wide))], names_to = "Sample", values_to = "Reads")
df_mock <- df_mock %>% filter(Reads>0) #filter all taxa with null abundance
df_wide=pivot_wider(df_mock,names_from = "Sample", values_from="Reads") #reshape for visibility 
View (df_wide)
#There is no false negative. The minimmal number of reads associated with a true mock genus is 135 (Pseudomonas in RUN3_barcode06). The maximal number of reads associated with a false positive is 5 (Enterobacter in RUN3_barcode06).

#Figure S1
exp_mock=aggregate(Rel_ab~Genus, data=Mock.expected, sum)
exp_mock$Sample="Expected"

df_mock_GENUS<-aggregate(Reads~Genus+Sample, data=df_mock, FUN=sum)

df_mock_nona<-df_mock_GENUS %>% filter(!is.na(Genus))
offset_vect_pre = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset_pre = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect_pre[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab_pre  = df_mock_nona$Reads / df_mock_nona$offset_pre

# The false positive with highest relative abundance is ENterobacter with a frequency of 0.00088

df_mock_nona<-subset(df_mock_nona, df_mock_nona$Rel_ab>0.005)

offset_vect = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab  = df_mock_nona$Reads / df_mock_nona$offset

gg_mock<-rbind(df_mock_nona[,c("Sample", "Genus", "Rel_ab")],exp_mock)
gg_mock %>% arrange(Rel_ab)

ggplot(gg_mock, mapping=aes(x=Sample, y=Rel_ab, fill=Genus))+
  geom_bar(width = 0.8, stat = "identity")+
  xlab("Mock samples")+ 
  ylab("Relative abundances")+
  theme_minimal()+
  theme(
    legend.position="right",
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(vjust = .5, hjust=1,angle =90))+
  scale_x_discrete(labels=c("Expected"="Expected", "RUN1_barcode10" = "RUN1", "RUN2B_barcode11" = "RUN 2B","RUN2C_barcode11" = "RUN 2C","RUN3_barcode06" = "RUN3"))

#ggsave("res/mock_barplot.svg")


#Working at Species level####
ps_mock_Species = tax_glom(ps_mock, "Species")

#Observing the number of reads of false positive and false negative taxa
#Setting as dataframe
df_taxo=as.data.frame(tax_table(ps_mock_Species)) #get taxonomy
df_taxo$tax_ID=rownames(df_taxo)
df_reads=as.data.frame(otu_table(ps_mock_Species)) #get reads
df_reads$tax_ID=rownames(df_reads)
df_wide=left_join(df_taxo, df_reads, by="tax_ID") #assemble taxonomy and reads
df_mock=pivot_longer(df_wide, cols = names(df_wide)[9:length(names(df_wide))], names_to = "Sample", values_to = "Reads")
df_mock <- df_mock %>% filter(Reads>0) #filter all taxa with null abundance
df_wide=pivot_wider(df_mock,names_from = "Sample", values_from="Reads") #reshape for visibility 
View (df_wide)
#There is no false negative. The minimmal number of reads associated with a true mock Species is 135 (Pseudomonas in RUN3_barcode06). The maximal number of reads associated with a false positive is 5 (Enterobacter in RUN3_barcode06).

#Figure S1
exp_mock=aggregate(Rel_ab~Species, data=Mock.expected, sum)
exp_mock$Sample="Expected"

df_mock_Species<-aggregate(Reads~Species+Sample, data=df_mock, FUN=sum)

df_mock_nona<-df_mock_Species %>% filter(!is.na(Species))
offset_vect_pre = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset_pre = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect_pre[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab_pre  = df_mock_nona$Reads / df_mock_nona$offset_pre

# The false positive with highest relative abundance is ENterobacter with a frequency of 0.00088

df_mock_nona<-subset(df_mock_nona, df_mock_nona$Rel_ab>0.005)

offset_vect = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab  = df_mock_nona$Reads / df_mock_nona$offset

gg_mock<-rbind(df_mock_nona[,c("Sample", "Species", "Rel_ab")],exp_mock)
gg_mock %>% arrange(Rel_ab)

ggplot(gg_mock, mapping=aes(x=Sample, y=Rel_ab, fill=Species))+
  geom_bar(width = 0.8, stat = "identity")+
  xlab("Mock samples")+ 
  ylab("Relative abundances")+
  theme_minimal()+
  theme(
    legend.position="right",
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(vjust = .5, hjust=1,angle =90))+
  scale_x_discrete(labels=c("Expected"="Expected", "RUN1_barcode10" = "RUN1", "RUN2B_barcode11" = "RUN 2B","RUN2C_barcode11" = "RUN 2C","RUN3_barcode06" = "RUN3"))

#ggsave("res/mock_barplot.svg")

#Working at Family level####
ps_mock_Family = tax_glom(ps_mock, "Family")

#Observing the number of reads of false positive and false negative taxa
#Setting as dataframe
df_taxo=as.data.frame(tax_table(ps_mock_Family)) #get taxonomy
df_taxo$tax_ID=rownames(df_taxo)
df_reads=as.data.frame(otu_table(ps_mock_Family)) #get reads
df_reads$tax_ID=rownames(df_reads)
df_wide=left_join(df_taxo, df_reads, by="tax_ID") #assemble taxonomy and reads
df_mock=pivot_longer(df_wide, cols = names(df_wide)[9:length(names(df_wide))], names_to = "Sample", values_to = "Reads")
df_mock <- df_mock %>% filter(Reads>0) #filter all taxa with null abundance
df_wide=pivot_wider(df_mock,names_from = "Sample", values_from="Reads") #reshape for visibility 
View (df_wide)
#There is no false negative. The minimmal number of reads associated with a true mock Family is 135 (Pseudomonas in RUN3_barcode06). The maximal number of reads associated with a false positive is 5 (Enterobacter in RUN3_barcode06).

#Figure S1
exp_mock=aggregate(Rel_ab~Family, data=Mock.expected, sum)
exp_mock$Sample="Expected"

df_mock_Family<-aggregate(Reads~Family+Sample, data=df_mock, FUN=sum)

df_mock_nona<-df_mock_Family %>% filter(!is.na(Family))
offset_vect_pre = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset_pre = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect_pre[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab_pre  = df_mock_nona$Reads / df_mock_nona$offset_pre

# The false positive with highest relative abundance is ENterobacter with a frequency of 0.00088

df_mock_nona<-subset(df_mock_nona, df_mock_nona$Rel_ab>0.005)

offset_vect = tapply(df_mock_nona$Reads, df_mock_nona$Sample, sum)
df_mock_nona$offset = unlist(sapply(1:length(df_mock_nona$Sample), function(i) offset_vect[df_mock_nona$Sample[[i]]]))
df_mock_nona$Rel_ab  = df_mock_nona$Reads / df_mock_nona$offset

gg_mock<-rbind(df_mock_nona[,c("Sample", "Family", "Rel_ab")],exp_mock)
gg_mock %>% arrange(Rel_ab)

ggplot(gg_mock, mapping=aes(x=Sample, y=Rel_ab, fill=Family))+
  geom_bar(width = 0.8, stat = "identity")+
  xlab("Mock samples")+ 
  ylab("Relative abundances")+
  theme_minimal()+
  theme(
    legend.position="right",
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(vjust = .5, hjust=1,angle =90))+
  scale_x_discrete(labels=c("Expected"="Expected", "RUN1_barcode10" = "RUN1", "RUN2B_barcode11" = "RUN 2B","RUN2C_barcode11" = "RUN 2C","RUN3_barcode06" = "RUN3"))

#ggsave("res/mock_barplot.svg")

#GARBAGE####
# OTU table #
# read_tsv(illumina_data,
#          na = c("", "NA", "No_hit", "0", "0.0"),
#          show_col_types = FALSE) %>%
#   pivot_longer(cols = contains("_S"),
#                names_to = "samples",
#                values_to = "reads",
#                values_drop_na = TRUE) %>% 
#   select(OTU, references, samples, reads) %>%
#   filter(!is.na(references)) %>%
#   separate_rows(references, sep = ",") %>%
#   group_by(OTU, samples) %>%
#   mutate(weight = 1 / n()) %>%
#   ungroup() %>%
#   mutate(reads = reads * weight) %>%
#   select(-OTU, -weight) %>% 
#   count(references, samples, wt=reads, name = "reads") %>%
#   pivot_wider(names_from = samples, values_from = reads, values_fill = 0) %>%
#   column_to_rownames(var = "references") -> my_otu_table


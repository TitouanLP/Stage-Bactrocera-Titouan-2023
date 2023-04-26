#*******************#
#### NGB - FRUIT ####
#*******************#

#Packages----
library(tidyverse)
library(phyloseq)
library(dplyr)
library(ggrepel)
library (vegan)

#Load data----
#**Nouvelle taxo----
silva_data = "databioinfo/tax_slv_ssu_138.1.txt.gz"
illumina_data = "dataFRUIT/NGB_16S_Illumina_214_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table3"

# illum_bis = "datatransmission/Transmision_MiSeq_16S_20230307_864_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2.gz"


##functions----
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
##main----

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
#**Metadata----
read.csv2("dataFRUIT/metadata_illumina.csv") %>%
  filter(experiment == "RUN4") %>%
  arrange(env, locality) %>% 
  add_count(Label) %>%
  mutate(Label = paste0(str_pad(row_number(),width=2, side="left",pad = "0"), "_", Label)) %>%
  select(-n) -> my_mtd
rownames(my_mtd) <- my_mtd$Label

# read.csv2("dataFRUIT/metadata_illumina.csv") %>% 
#   filter(experiment == "RUN4") %>% 
#   add_count(Label) %>% 
#   mutate(Label = paste0(str_pad(row_number(), width = 2, side = "left", pad = "0"), "_", Label)) %>% 
#   select(-n) %>% 
#   column_to_rownames("Label") %>% 
#   mutate(.before = 4, Label = rownames(.)) -> my_mtd


#Basic statistics----
my_mtd %>%
  select(sample_ID) %>%
  pull(sample_ID) -> all_RUN4_samples

my_mtd %>%
  filter(type=="ENV") %>% 
  select(sample_ID) %>%
  pull(sample_ID) -> only_environmental_samples 

#Nombre de séquences
columns_to_keep <- c("OTU","references",only_environmental_samples)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%  
  select(all_of(columns_to_keep)) %>%
  mutate(.before = 1, total = rowSums(.[,3:ncol(.)],na.rm=T)) %>% 
  filter(total > 0) %>% 
  select(-total) -> illumina_data_RUN4 

illumina_data_RUN4 %>% dim()
#Dans le RUN4, il y a 559 OTU

#Nombre de séquences non assignées à quoi que ce soit : 
illumina_data_RUN4 %>%
  filter(is.na(references)) 
#Dans le RUN4, il y a 1 seule OTU associée à aucune référence. Elle a très peu de reads.

#*******************----
#Building a phyloseq object for Illumina data (based on OTUs, not on taxonomy)----
#toutes les données dans même objet

#**OTU pour tax ancetre commun en vue de decontamination----
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

#**Tax pour décontamination (sans séparer les ref)----
read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% View()
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

#OTUs à retirer selon microDecon
decontaminated$OTUs.removed %>% 
  pull(OTU)-> otus_to_remove
length(otus_to_remove) #il y en a 95

my_tax %>% 
  filter(OTU %in% otus_to_remove) 

#OTUs seulement présentes dans les puits négatifs
my_otu %>% 
  mutate(.before = 1, total_env = rowSums(.[,env_samples])) %>%
  filter(total_env==0) %>% 
  rownames(.)-> OTU_absent_from_env_samples
length(OTU_absent_from_env_samples) #Il y en a 81

#OTUs dont les reads sont corrigés
decontaminated$reads.removed %>%
  replace(.,is.na(.),0) %>%
  filter(!OTU %in% otus_to_remove) %>% 
  filter(!OTU %in% OTU_absent_from_env_samples) %>% 
  pull(OTU) -> partially_removed 
length(partially_removed) #Il y en a 19

#**Phyloseq object pour faire graph OTU conta
otu_table(my_otu,taxa_are_rows = TRUE) -> OTU.decon 

tax_table(as.matrix(my_tax)) -> TAX.decon

sample_data(my_mtd) -> MTD

ps_decon = phyloseq(OTU.decon, TAX.decon, MTD)

#Distribution of contaminant and non-contaminant OTUs in environmental samples and controls
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
  filter(OTU == "44")

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

#clearNames(ps.pa, ps.pa.neg, ps.pa.pos, ps.decon, ps.decon.neg, ps.decon.pos)

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

#**Observation des données----
#Profondeurs de séquençage par échantillon (plein de graphes)
my_cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("Label") %>% 
  rename(total_reads = 2) %>% 
  left_join(my_mtd_env, ., by = "Label") -> my_mtd_env

summary(my_mtd_env$total_reads)

ggplot(my_mtd_env, aes(log10(total_reads), fill=type))+
  geom_histogram()+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Total reads", y = "Frequency")

# ggsave("res/RUN4_Distribution_of_bacterial_reads.pdf")

my_mtd_env <- my_mtd_env %>% arrange(total_reads)
sample_order<-my_mtd_env$sample_ID
label_order<-my_mtd_env$Label

ggplot(my_mtd_env, aes(x=sample_ID, y=total_reads, col=sex))+
  geom_point(size=3)+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order, labels=label_order)

# ggsave("res/RUN4_Distribution_of_total_reads.pdf")

ggplot(my_mtd_env, aes(x=env, y=total_reads))+
  geom_violin(aes(fill=env))

ggplot(my_mtd_env, aes(x=env, y=total_reads))+
  geom_point(aes(shape=sex, col = locality), size = 3)

ggplot(my_mtd_env, aes(x=locality, y=total_reads))+
  geom_violin(aes(fill=env))+
  scale_x_discrete(limits = c("ANS","ETA","VER","REL","GUI","MAR","PAL","VID","BER","BOU","CIR"))

ggplot(my_mtd_env, aes(x=locality, y=total_reads))+
  geom_point(size = 3, aes(col=sex))+
  scale_x_discrete(limits = c("ANS","ETA","VER","REL","GUI","MAR","PAL","VID","BER","BOU","CIR"))

#Proportions of unassigned reads
columns_to_keep <- c("OTU", "references", only_environmental_samples)

read_tsv(illumina_data,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(all_of(columns_to_keep)) %>% View() 
  mutate_at(., c(3:ncol(.)),~replace(.,is.na(.),0)) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>%
  filter(is.na(references)) -> unassigned
#il n'y a qu'une seule OTU non assignée. OTU6. Si on la blaste on trouve une mitochondrie de Bactrocera. Pas très anormal. Elle a un total de 139 reads. Pas de quoi perturber l'analyse. On laisse tomber l'analyse des non assignées.

#**Barplot - Composition bactérienne par sample----
#Une façon de faire basée sur les abondances relatives
as.data.frame(my_tax_distrib) %>%  #get taxonomy
  rownames_to_column(var = "tax_ID") -> df 
as.data.frame(my_otu_distrib) %>% #get reads 
  rownames_to_column(var = "tax_ID") %>% 
  left_join(df, ., by="tax_ID") %>% 
  select(-tax_ID) -> df
df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  left_join(., my_mtd_env, by = c("Sample" = "Label")) -> df_with_NA

# my_mtd_env <- my_mtd_env %>% 
#   arrange(env, locality, sex)
sample_order<-my_mtd_env$Label
label_order<-my_mtd_env$Label

# my_mtd_env %>% 
#   arrange(env, locality, sex) 

ggplot(df_with_NA, aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  # scale_x_discrete(labels = label_order)+
  facet_wrap(~locality, scales="free_x")

my_mtd_env %>% 
  arrange(locality) %>% 
  pull(Label)->samples_bon_ordre 

my_mtd_env %>% 
  arrange(locality) %>%
  pull(locality)->etiquettes_bon_ordre

ggplot(df_with_NA, aes(x = Sample, y = Rel_ab, fill = family))+
  geom_bar(stat = "identity", aes(x = Sample, y = Rel_ab) )+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  scale_x_discrete(limits = samples_bon_ordre, label = etiquettes_bon_ordre)

#Afficher les abondances par taxon
#Niveau genre
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
  ) %>% View
#ça devient les nouveau OTU non ? 

#Niveau famille
df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i){
    paste(domain[i], phylum[i],class[i], order[i], family[i], sep="_") 
  })) %>% 
  group_by(tax_ID, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>%  
  pivot_wider(id_cols = 1, names_from = Sample, values_from = Rel_ab, values_fill = NA) %>% view() 

#STATS#####
#NMDS----

#Avec matrice présence/abnsence
#Env/sex
ps_rare = rarefy_even_depth(ps, sample.size = 24000)  
ps_bin = transform_sample_counts(ps_rare, function(x) 1*(x>0))
otu_table(ps) %>% 
  as.data.frame() %>% 
  summarise(across(everything(), sum)) %>% 
  select_if(~ sum(.x) < 25000) %>% View()


jaccard <- ordinate(ps_bin, "PCoA", "jaccard") 

col = c("red","Blue","Yellow","Orange","Green")
plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "env",
                shape = "sex") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "sex", values=c(17,19)) +
  scale_color_manual(name = "env", values = col) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3,max.overlaps = 65)

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

#ggsave("res/NMDS_flies_plot.svg", width=15, height = 10, units = "cm")

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

#ggsave("res/NMDS_bacteria_plot.svg", width = 15,  height=10, units = "cm")

#NMDS abondance
#env/localité
ncol = 3
col_env <- paletteer_d("ggsci::default_igv")[seq(ncol, ncol + 10, 1)]

my_mtd_env$locality %>% unique() %>% length()

ps_rare = rarefy_even_depth(ps, sample.size = 24000)

bray  <- ordinate(ps_rare, "NMDS", "bray")
  
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
    scale_colour_manual(name = "Locality", values = col_env) +
    geom_text_repel(aes(label = env, fontface = 'italic'), size = 3, max.overlaps = 100)

#env/sex

plot_ordination(ps_rare,
                bray,
                type = "samples",
                color = "env",
                shape = "sex") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) +
  scale_shape_manual(name = "sex", values=c(17,19)) +
  scale_color_manual(name = "env", values = col_env) +
  geom_text_repel(aes(label = env, fontface = 'italic'), size = 3,max.overlaps = 100)

#taxon
  
plot_ordination(ps_rare,
                bray,
                type = "taxa",
                color = "class") +
  theme_bw() +
  geom_point(size = 2) +
  theme(legend.position = "right") +
#  scale_color_manual(values = c(col_tax[1], col_tax[3], col_tax[2], col_tax[4:23])) +
#  scale_shape_manual(values = c(19, 17, 3)) +
geom_text_repel(aes(label = family, fontface = 'italic'), size = 2.5)

#----4. DIVERSITY ANALYSES (Figure 2) ------
#a. Rarefaction curves (Figure S2)------
library("iNEXT")
load("datafruit/diversities.Rdata")

str(out)

my_mtd_env$Label

intersect(out$iNextEst$size_based$Assemblage, my_mtd_env$Label)

out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:72]) %>%
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage))+
  geom_line()+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()
  
out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:15]) %>%
  filter(!Method == "Observed") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_point(size = 3)+
  geom_line()+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

p1 <- out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:15]) %>%
  filter(Method == "Rarefaction") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_point(size = 3)+
  geom_line()+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

p2 <- out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:15]) %>%
  filter(Method == "Extrapolation") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_line(linetype = 3)+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

p1+p2

out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_env$Label[1:8]) %>%
  ggplot(., aes(x= m, y = qD, col = as.factor(Order.q)))+
  geom_line()+
  facet_wrap(~Assemblage,  ncol = 4)+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

plot_richness(ps, measures="Shannon", color = "env")+
  facet_wrap(~env, scales ="free_x")

plot_richness(ps, measures="Shannon", color = "locality")+
  facet_wrap(~locality, scales ="free_x")

data.frame(
  Sample=sample_names(ps),
  Observed = estimate_richness(ps, measures = "Observed"),
  Shannon = exp(estimate_richness(ps, measures = "Shannon")),
  Simpson = 1 / (1 - estimate_richness(ps, measures = "Simpson")),
  Chao1 = estimate_richness(ps, measures = "Chao1"),
  Total_reads = sample_sums(ps)
) %>% 
  left_join(., my_mtd_env, by=c("Sample" = "Label")) %>%  
arrange(env) -> df

my_mtd_env %>% 
  arrange(env) %>% 
  pull(Label)->samples_bon_ordre 

ggplot(df, aes(x=Sample, y = Shannon, fill=env))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  scale_x_discrete(limits = samples_bon_ordre)

my_mtd_env %>% 
  arrange(locality) %>% 
  pull(Label)->samples_bon_ordre 

ggplot(df, aes(x=Sample, y = Shannon, fill=locality))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  scale_x_discrete(limits = samples_bon_ordre)

ggplot(df, aes(x=Sample, y = Shannon, fill=env))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  #scale_x_discrete(limits = samples_bon_ordre)+
  facet_grid(.~env, scales="free_x", shrink = F)

ggplot(df, aes(x=env, y = Shannon, fill=env))+
  geom_boxplot()+
  geom_point(size = 4, aes(shape=df$locality))+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))
#scale_x_discrete(limits = samples_bon_ordre)+
#facet_grid(.~env, scales="free_x", shrink = F)

ggiNEXT(out, type=1, facet.var="None", color.var = "Assemblage")+
  # facet_wrap(env, scales="free", ncol=4)+
  xlab("Number of reads") +
  ylab("Number of genera") +
  theme_minimal()+
  theme(legend.position = "right")+
  guides(shape = "none") +
  scale_shape_manual()

unique_colors <- viridis_pal()(length(unique(out$Data$Assemblage)))
out$DataInfo %>% View()
out$AsyEst %>% View()
out$iNextEst$coverage_based

ggiNEXT(out, type=1, facet.var="None", color.var = "Assemblage") +
  xlab("Number of reads") +
  ylab("Number of genera") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(values = unique_colors)


#b. Diversities------
#Gamma diversity
library(vegan)

my_ps <- ps_genus_nmds
df_taxo = as.data.frame(tax_table(my_ps))
df_taxo$tax_ID = rownames(df_taxo)
df_reads = as.data.frame(otu_table(my_ps))
df_reads$tax_ID = rownames(df_reads)
df_wide = left_join(df_taxo, df_reads, by = "tax_ID") 
df_wide$Species <- NULL 
# View(df_wide)
BF <- df_wide[, c(7, 8:80)]
if (length(which(BF$genus == "uncultured")) == 2)
  BF$genus[which(BF$genus == "uncultured")] <-
  c("uncultured_1", "uncultured_2")
BF$genus[which(BF$genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")] <-
  "Rhizobium"
rownames(BF) <- BF$tax_ID
BF$genus <- NULL
BF$tax_ID <- NULL
# colnames(BF) <- mtd$label[match(colnames(BF), mtd$sample_ID)]
FB <- t(BF)
Log_FB = log10(FB + 1)

FB <- round(FB, digits = 0)
res_spl<-multipart(y=FB, scales = 1, global = TRUE, relative = FALSE, nsimul=99)
res_spl$statistic

#Alpha diversities by sample (Figure 2)
library(iNEXT)

load("diversities.Rdata")

div_df<-out$AsyEst[out$AsyEst$Diversity=="Shannon diversity",]

div_df<-left_join(div_df, mtd, by = c("Assemblage" = "Label"))

names(div_df)

plot_heatmap(ps_genus_nmds, "RDA", "none", sample.label="env", taxa.label="family")

#--------- LBM Clustering ------
library(sbm)
LBM_log_reads_gaussian <-
  as.matrix(Log_FB) %>%
  estimateBipartiteSBM(model = 'gaussian') #associe les bactéries qui s'associent de façon similaire avec les échantillons

save(LBM_log_reads_gaussian, )

load("datafruit/LBM_log_reads_gaussian.Rdata")

str(LBM_log_reads_gaussian)

memb_spl_obs <- LBM_log_reads_gaussian$memberships$row

LBM_log_reads_gaussian$nbBlocks

Qrow <- LBM_log_reads_gaussian$nbBlocks["row"][[1]]
Qcol <- LBM_log_reads_gaussian$nbBlocks["col"][[1]]

memb_spl_obs <- LBM_log_reads_gaussian$memberships$row

LBM_log_reads_gaussian$storedModels %>% arrange(ICL)

SBMobject <- LBM_log_reads_gaussian

BF_mat <- as.matrix(as.data.frame(otu_table(ps)))

BF_mat %>% view()

df = as.data.frame(BF_mat)
df$taxon = rownames(df)
df_long <-
  pivot_longer(df,
               cols = names(df)[-ncol(df)],
               names_to = "Label",
               values_to = "reads")
df_long <- left_join(df_long, my_mtd_env, by = "Label")

Qrow <- SBMobject$nbBlocks["row"][[1]]
Qcol <- SBMobject$nbBlocks["col"][[1]]

df_tax_cluster = data.frame(taxon = rownames(BF_mat),
                            tax_cluster = SBMobject$memberships$col)

df_long <- left_join(df_long, df_tax_cluster, by = "taxon")

df_long <- df_long %>% arrange(tax_cluster)

correct_taxon_order <- unique(df_long$taxon)

df_spl_cluster = data.frame(label = colnames(BF_mat)
                            ,
                            spl_cluster = SBMobject$memberships$row)
df_long <- left_join(df_long, df_spl_cluster, by = "label")

df_long <- df_long %>% arrange(spl_cluster)
correct_sample_order <- unique(df_long$label)

mat <- ggplot(df_long, aes(
  x = label,
  y = taxon,
  fill = log10(reads + 1)
)) +
  coord_flip() +
  geom_tile() +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(
      angle = 90,
      vjust = 1,
      hjust = 1,
      face = "italic"
    ),
    axis.text.y = element_text()
  ) +
  scale_fill_gradient(low = "white", high = "Black") +
  scale_y_discrete(name = "Bacterial genus", limits = correct_taxon_order) +
  scale_x_discrete(name = "Sample", limits = correct_sample_order) +
  geom_vline(
    xintercept = 0.5 + cumsum(table(SBMobject$memberships$row))[-length(table(SBMobject$memberships$row))],
    color = "red",
    size = 1
  ) +
  geom_hline(
    yintercept = 0.5 + cumsum(table(SBMobject$memberships$col))[-length(table(SBMobject$memberships$col))],
    color = "red",
    size = 1
  )
mat

myOrderedMatrixPlot(LBM_log_reads_gaussian, BF)

ggsave("res/weigthed_incid_mat.svg")

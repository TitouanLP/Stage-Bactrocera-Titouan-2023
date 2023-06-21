#*******************#
#### NGB - TRANSMISSION ####
#*******************#

#Packages
library(tidyverse)
library(phyloseq)
library("iNEXT")
library("microDecon")
library("ggrepel")

#Load data----
#**Nouvelle taxo----
silva_data = "databioinfo/tax_slv_ssu_138.1.txt.gz"
illumina_data_trans = "datatransmission/Transmision_MiSeq_16S_20230307_864_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2.gz"

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
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%  dim()

#Nombre de séquences non assignées à quoi que ce soit : 80
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  filter(is.na(references)) %>% dim()

#Nombre de séquences avec une assignation multiple : 690 (622 assignations distinctes)
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  filter(!is.na(references)) %>% 
  filter(str_detect(references, ",")) %>% 
  select(references) %>% duplicated() %>% table()

#*******************----
#Building a phyloseq object for Illumina data (based on OTUs, not on taxonomy)----
#toutes les données dans même objet

#**Metadata----
#Formatage de my_mtd pour phyloseq
read.csv2("datatransmission/mtd_transmission_cor.csv") %>% 
  mutate(row = ID) %>% 
  column_to_rownames(., var = "row") -> my_mtd_trans

#**OTU pour décontamination----

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% #7275 OTU
  filter(!is.na(references)) %>% #7244 OTU (donc que 31 mal assignés)
  mutate_all(~replace(., is.na(.), 0)) %>% #enlève les TRUE FALSE par 0/1
  select(-c(742:805)) %>%
  select(-c(744:797)) %>% 
  select(-c(total,cloud,amplicon,length,abundance,chimera,spread,quality,sequence,identity,taxonomy)) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>% 
  select(-total,-references) %>% 
  column_to_rownames(var = "OTU") -> my_otu_trans

#**Tax pour décontamination (sans séparer les ref)----

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(OTU,taxonomy) %>% 
  filter(!is.na(taxonomy)) %>% 
  filter(OTU %in% rownames(my_otu_trans)) %>% 
  separate_rows(taxonomy, sep = "\\|") %>% 
  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy" = "taxonomic_path")) %>% 
  filter(taxonomic_rank %in% good_taxonomic_levels) %>% #enlève les mauvais niveaux taxo dont les NA
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>% 
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy) %>%  
  column_to_rownames(var = "OTU") %>% 
  cbind(.,OTU = rownames(my_otu_trans)) -> my_tax_trans

# options(max.print = 10000)
# rownames(my_otu_trans) == rownames(my_tax_trans) 

#*******************----
#MicroDecon----
# install.packages("remotes")
# remotes::install_github("donaldtmcknight/microDecon")
library("microDecon")

neg_samples <- my_otu_trans %>% #enlever mock 
  select(731:738) %>% 
  summarise(across(everything(), sum)) 

neg_samples <- my_otu_trans %>%
  select(731:738) %>%
  as.data.frame()

env_samples <- my_otu_trans %>% #enlever mock 
  select(1:728) %>% 
  summarise(across(4:ncol(.), sum)) %>% 
  select_if(~ sum(.) < 800) 

env_samples <- my_otu_trans %>% #enlever mock 
  select(1:728) %>% 
  as.data.frame() 

data.frame(OTU = rownames(my_otu_trans)) %>% 
  cbind(., neg_samples) %>% 
  cbind.data.frame(., env_samples) -> df_decon

decontaminated <- decon(data = df_decon,numb.blanks=length(colnames(neg_samples)),numb.ind=length(colnames(env_samples)), taxa=F) # je suis obligé d'enlever les NEG de PCR car en TRUE FALSE ou si on change c'est 1/0 donc trop faible 

# save(decontaminated, file="deconta_trans.Rdata")
load("datatransmission/deconta_trans.Rdata")

#OTUs à retirer selon microDecon
decontaminated$OTUs.removed %>% 
  pull(OTU)-> otus_to_remove
length(otus_to_remove) #il y en a 433

my_tax_trans %>% 
  filter(OTU %in% otus_to_remove) 

#OTUs seulement présentes dans les puits négatifs
my_otu_trans %>% 
  mutate(.before = 1, total_env = rowSums(.[, 1:728])) %>%
  filter(total_env == 0) %>% 
  rownames(.) -> OTU_absent_from_env_samples
length(OTU_absent_from_env_samples) #Il y en a 58

#OTUs dont les reads sont corrigés
decontaminated$reads.removed %>%
  replace(.,is.na(.),0) %>%
  filter(!OTU %in% otus_to_remove) %>% 
  filter(!OTU %in% OTU_absent_from_env_samples) %>% 
  pull(OTU) -> partially_removed 
length(partially_removed) #Il y en a 66

#**Phyloseq object pour faire graph OTU conta
my_tax_trans %>%
  filter(!domain=="Bacteria") %>%
  rownames() #29 OTU qui ne sont pas des bactéries

my_tax_trans %>% 
  filter(domain=="Bacteria") %>% 
  rownames() -> OTU_to_keep

filter(my_otu_trans, rownames(my_otu_trans) %in% OTU_to_keep) -> my_otu_trans
filter(my_tax_trans, rownames(my_tax_trans) %in% OTU_to_keep) -> my_tax_trans

otu_table(my_otu_trans,taxa_are_rows = TRUE) -> OTU.decon 

tax_table(as.matrix(my_tax_trans)) -> TAX.decon

sample_data(my_mtd_trans) -> MTD

ps_decon = phyloseq(OTU.decon, TAX.decon, MTD)

sample_names(OTU.decon) == sample_names(MTD)

#Distribution of contaminant and non-contaminant OTUs in environmental samples and controls----
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

#**Decontamination de la table d'OTU----
decontaminated$decon.table %>% 
  mutate(OTU = as.numeric(OTU)) %>% 
  arrange(OTU) %>% 
  mutate(total = rowSums(.[,3:ncol(.)])) %>%
  filter(total > 0) %>%
  select(-total,-Mean.blank) %>% 
  column_to_rownames(var = "OTU") %>% 
  summarise(across(everything(), sum)) %>% 
  select_if(~ sum(.) == 0) %>% 
  names() -> cols_to_remove # 4 samples ou il n'y a plus de reads 

decontaminated$decon.table %>%  
  mutate(OTU = as.numeric(OTU)) %>% 
  arrange(OTU) %>% 
  select(-all_of(cols_to_remove)) %>% 
  mutate(.before = 1 ,total = rowSums(.[,3:ncol(.)])) %>%
  filter(total > 0) %>% 
  select(-total,-Mean.blank) %>%
  column_to_rownames(var = "OTU") -> cleaned_otu_table

cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% as.data.frame() %>% 
  pull(1)-> total_reads_per_sample

hist(total_reads_per_sample)
summary(total_reads_per_sample)

#**Table de références---- 
my_tax_trans %>% 
  filter(OTU %in% rownames(cleaned_otu_table)) -> cleaned_tax_table #POURQUOI C PAS LE MEME NOMBRE QUE OTU TABLE 

distrib_OTU <- function(otu_tab){
  read_tsv(illumina_data_trans,
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

distrib_OTU(cleaned_otu_table) -> my_otu_distrib

#**TAX table après décontamination----
distrib_tax <- function(otu_tab) {
  read_tsv(
    illumina_data_trans,
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

distrib_tax(cleaned_otu_table) -> my_tax_distrib

#150 ref à autres choses que bactéries 
my_tax_distrib %>%
  filter(!domain == "Bacteria") %>%
  rownames()

my_tax_distrib %>% 
  filter(domain=="Bacteria") %>% 
  rownames() -> OTU_to_keep

filter(my_otu_distrib, rownames(my_otu_distrib) %in% OTU_to_keep) -> my_otu_distrib
filter(my_tax_distrib, rownames(my_tax_distrib) %in% OTU_to_keep) -> my_tax_distrib

#**phyloseq objet----
my_mtd_trans %>%
  filter(ID %in% colnames(cleaned_otu_table)) %>% 
  filter(type == "ENV") -> my_mtd_trans 
sample_data(my_mtd_trans) -> MTD
  
otu_table(cleaned_otu_table,taxa_are_rows = T) -> OTU
tax_table(as.matrix(cleaned_tax_table)) -> TAX

ps = phyloseq(OTU, TAX, MTD)
  
otu_table(my_otu_distrib,taxa_are_rows = T) -> OTU
apply(my_tax_distrib, c(1,2), as.character) %>%
  tax_table(my_tax_distrib) -> TAX

sample_names(OTU)==sample_names(MTD)

ps_multiple_trans = phyloseq(OTU, TAX, MTD)

#**Observation des données----
#Total Reads 
cleaned_otu_table %>% 
  summarise(across(everything(), sum)) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("ID") %>% 
  rename(total_reads = 2) %>% 
  left_join(my_mtd_trans, ., by = "ID") -> my_mtd_trans

t(cleaned_otu_table) %>% 
  as.data.frame() %>% 
  mutate(.before = 1, total_reads = rowSums(.[,1:ncol(.)],na.rm=T)) %>% View()

summary(my_mtd_trans$total_reads)
which(is.na(my_mtd_trans$total_reads))

ggplot(my_mtd_trans, aes(log10(total_reads), fill=type))+
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black")+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Total reads", y = "Frequency")

ggplot(my_mtd_trans, aes(total_reads, fill=type))+
  geom_histogram(fill = "lightblue", color = "black")+
  theme_minimal()+
  labs(title="Effet fruit - RUN4 - Illumina",
       x ="Total reads", y = "Frequency")

# ggsave("res/RUN4_Distribution_of_bacterial_reads.pdf")

my_mtd_trans <- my_mtd_trans %>% arrange(total_reads)
sample_order <- my_mtd_trans$Sample_ID

ggplot(my_mtd_trans, aes(x=Sample_ID, y=total_reads, col=Stage))+
  geom_point()+
  theme_minimal()+
  labs(title="Effet transmission - Illumina",
       x ="Bacterial reads", y = "Frequency")+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.4))+
  scale_x_discrete(limits=sample_order)

# ggsave("res/RUN4_Distribution_of_total_reads.pdf")

Stage_order <- c("egg","larva","pupa","adult","mother")

ggplot(my_mtd_trans, aes(x=Stage, y=total_reads))+
  geom_violin(aes(fill=Stage)) + 
  scale_x_discrete(limits=Stage_order) +
  scale_fill_manual(values = c("egg" = "#E6194B", "larva" = "#3CB44B", "pupa" = "#FFE119", "adult" = "#4363D8", "mother" = "#F58231"),
                    labels = c("egg", "larva", "pupa", "adult", "mother")) + 
  guides(fill = guide_legend(override.aes = list(fill = c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231"))))

ggplot(my_mtd_trans, aes(x=Stage, y=total_reads))+
  geom_violin(aes(fill=protein)) + 
  scale_x_discrete(limits=Stage_order)

ggplot(my_mtd_trans, aes(x=Stage, y=total_reads))+
  geom_violin(aes(fill=culture)) + 
  scale_x_discrete(limits=Stage_order)

ggplot(my_mtd_trans, aes(x=Stage, y=total_reads))+
  geom_violin(aes(fill=replicate_culture)) + 
  scale_x_discrete(limits=Stage_order)

ggplot(my_mtd_trans, aes(x=Stage, y=total_reads))+
  geom_point(aes(shape=Stage, col = protein), size = 3)

#Unassigned Reads
my_mtd_trans %>%
  pull(ID) -> only_environmental_samples

columns_to_keep <- c("OTU", "references", only_environmental_samples)

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(all_of(columns_to_keep)) %>% 
  mutate_at(., c(3:ncol(.)),~replace(.,is.na(.),0)) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>%
  filter(is.na(references)) -> unassigned

#**Barplot - Composition bactérienne par sample----
#Une façon de faire basée sur les abondances relatives
as.data.frame(my_tax_distrib) %>%  #get taxonomy
  rownames_to_column(var = "tax_ID") -> df 
as.data.frame(my_otu_distrib) %>% #get reads 
  rownames_to_column(var = "tax_ID") %>% 
  left_join(df, ., by="tax_ID") %>% 
  select(-tax_ID) -> df
df %>% 
  select(-`CTL-DNA-Moc-B_S761`,-`CTL-DNA-Moc-A_S377`) %>% 
  mutate(family = as.character(family)) -> df

df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  left_join(., my_mtd_trans, by = c("Sample" = "ID")) %>% 
  group_by(family, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>% 
  pivot_wider(
    id_cols = 1,
    names_from = Sample,
    values_from = Rel_ab,
    values_fill = NA
  ) %>% View() #-> df_with_NA_family

df_with_NA %>% head(n=100000) %>% View()

df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  left_join(., my_mtd_trans, by = c("Sample" = "ID")) %>% 
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

df_with_NA %>% 
  arrange(Stage) %>% 
  pull(Sample) -> label_order

colors <- c("#FF5F5F", "#FF8C42", "#FFC168", "#FFEE8D", "#D1FFA9", "#A4FFC5", "#77FFD2", "#4AFFD9", "#1DFFFF", "#79E7FF", "#A4B3FF", "#C4A6FF", "#FFA6FF", "#FFD6FF", "#FFC6A9")

df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(family, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>% 
  ungroup() %>%
  group_by(family) %>% 
  mutate(family = ifelse(all(Rel_ab < 0.1), "Autres", as.character(family))) %>% 
  ungroup() %>% 
  left_join(., my_mtd_trans, by = c("Sample" = "ID")) %>%
  mutate(Stage = factor(Stage, levels = c("egg", "larva", "pupa", "adult", "mother"),
                        labels = c("Oeufs",
                                   "Larves", "Nymphes", "Adultes", "Mères"))) %>% 
  ggplot(., aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(aes(fill = as.character(family)), stat = "identity")+
  theme_minimal()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box = "vertical",
        legend.key.size = unit(0.6, "cm"),
        axis.text.x = element_text(angle=90, vjust=0.4, hjust=1.2),
        legend.text = element_text(size = 10))+
  guides(fill = guide_legend(ncol = 1))+
  scale_fill_manual(values = colors) +
  # scale_x_discrete(labels = label_order)+
  labs(
    x = "Echantillons environnementaux",
    y = "Abondances réelles",
    fill = "Familles"
  ) +
  facet_wrap(~Stage, scales="free_x", ncol = 5)


palette <- colorRampPalette(colors = c("blue", "green", "yellow", "red"))(num_colors)

ggplot(df_with_NA, aes(x = Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1.2)) +
  facet_wrap(~Stage, scales = "free_x")

ggplot(df_with_NA, aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  facet_wrap(~Stage, scales="free_x")

# Vérification des types de données uniques après la conversion
unique_types <- unique(sapply(df_with_NA$family, class))
print(unique_types)


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

#STATS#####
#NMDS----

#Avec matrice présence/abnsence
#Stage
ps_rare = rarefy_even_depth(ps, sample.size = 1000)  
ps_bin = transform_sample_counts(ps_rare, function(x) 1*(x>0))
otu_table(ps) %>% 
  as.data.frame() %>% 
  summarise(across(everything(), sum)) %>% 
  select_if(~ sum(.x) < 1000) %>% View() # max read est à 47000 et 3 en dessous de 1000 et 4 entre 0 et 1 mais bizzare 
sample_data(ps_bin) %>% View()


jaccard <- ordinate(ps_bin, "PCoA", "jaccard") 

col = c("red","Blue","#FFD6FF","Orange","#A4FFC5")
colors <- c("#FF5F5F", "#FF8C42", "#FFC168", "#FFEE8D", "#D1FFA9", "#A4FFC5", "#77FFD2", "#4AFFD9", "#1DFFFF", "#79E7FF", "#A4B3FF", "#C4A6FF", "#FFA6FF", "#FFD6FF", "#FFC6A9")
col <- c("#FF5F5F","#FFEE8D", "#D1FFA9","#79E7FF", "#FFA6FF")
col = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#F0E442", "#00BA38", "#619CFF", "#D55E00", "#CC79A7", "#999999")
col = c("purple", "darkgreen", "red","gray","orange")
ps_bin <- reorder(ps_bin, "Stage")
plot_ordination(ps_bin,
                jaccard,
                type = "samples",
                color = "Stage") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",
  legend.key.size = unit(1, "cm"))+
  theme(legend.text = element_text(face = "italic", size = 14),
        legend.title = element_text(size = 14)) +
  # scale_shape_manual(name = "protein", values=c(17,19,21)) +
  scale_color_manual(name = "Stade de développement", values = col,
                     labels = c("Adulte","Oeuf","Larve","Mère","Pupe"),
                     guide = guide_legend(override.aes = list(size = 4.5)))# +
  geom_text_repel(aes(label = Stage, fontface = 'italic'), size = 3,max.overlaps = 407)


bray <- ordinate(ps_rare, "PCoA", "bray") 
plot_ordination(ps_rare,
                bray,
                type = "samples",
                color = "Stage") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right",
        legend.key.size = unit(1, "cm"))+
  theme(legend.text = element_text(face = "italic", size = 14),
        legend.title = element_text(size = 14)) +
  # scale_shape_manual(name = "protein", values=c(17,19,21)) +
  scale_color_manual(name = "Stade de développement", values = col,
                     labels = c("Adulte","Oeuf","Larve","Mère","Pupe"),
                     guide = guide_legend(override.aes = list(size = 4.5)))# +
geom_text_repel(aes(label = Stage, fontface = 'italic'), size = 3,max.overlaps = 407)


tax_table(ps_bin) %>% View()

plot_ordination(ps_bin,
                jaccard,
                type = "taxa",
                color = "phylum") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) #+
  # scale_shape_manual(name = "protein", values=c(17,19,21)) +
  # scale_color_manual(name = "family", values = col) +
  #geom_text_repel(aes(label = class, fontface = 'italic'), size = 3,max.overlaps = 2809)

ps_bin_family = tax_glom(ps_bin, taxrank = "family") 

jaccard <- ordinate(ps_bin_family, "NMDS", "jaccard") 

plot_ordination(ps_bin_family,
                jaccard,
                type = "samples",
                color = "Stage") +
  theme_bw() +
  geom_point(size = 3) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(face = "italic")) #+
# scale_shape_manual(name = "protein", values=c(17,19,21)) +
# scale_color_manual(name = "family", values = col) +
#geom_text_repel(aes(label = class, fontface = 'italic'), size = 3,max.overlaps = 2809)

plot_ordination(ps_bin_family,
                jaccard,
                type = "taxa",
                color = "family") +
  theme_bw() +
  geom_point(size = 3) +
theme(legend.position = "bottom")
  
# comparer attendu avec mock avec observé ou faire à ma manière 
# dans full join : prendre objet phylo final et recup échantillon mock, agglo au niveau famille (ou genre) et faire pourcentage 
# rajouter en nombre de reads

#----4. DIVERSITY ANALYSES (Figure 2) ------
#a. Rarefaction curves (Figure S2)------
library("iNEXT")

# tab <- round(as.data.frame(otu_table(ps)), 0)

# out <- iNEXT(tab, q=c(0, 1, 2), datatype="abundance", size=c(500, 1000, 5000, 10000, 25000, 50000))
# 
# save(out, file ="datatransmission/out.Rdata")
load("datatransmission/out.Rdata")

str(out)

intersect(out$iNextEst$size_based$Assemblage, my_mtd_trans$ID)

out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_trans$ID[1:724]) %>%
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage))+
  geom_line()+
  theme(legend.position="none")+ 
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")

out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_trans$ID[1:15]) %>%
  filter(!Method == "Observed") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_point(size = 3)+
  geom_line()+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme(legend.position="none")

p1 <- out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_trans$ID[1:15]) %>%
  filter(Method == "Rarefaction") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_point(size = 3)+
  geom_line()+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

p2 <- out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_trans$ID[1:15]) %>%
  filter(Method == "Extrapolation") %>% 
  filter(Order.q == 0) %>% 
  ggplot(., aes(x= m, y = qD, col = Assemblage, shape = Method))+
  geom_line(linetype = 3)+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

out$iNextEst$size_based %>% 
  filter(Assemblage %in% my_mtd_trans$ID[20:35]) %>%
  ggplot(., aes(x= m, y = qD, col = as.factor(Order.q)))+
  geom_line()+
  facet_wrap(~Assemblage,  ncol = 4)+
  xlab("Nomdre de reads échantillonnés")+
  ylab("Number of OTUs")+
  theme_minimal()

plot_richness(ps, measures = "Shannon", color = "Stage") +
  facet_grid(~reorder(Stage, match(Stage, c("egg", "pupa", "larva", "adult", "mother"))), 
             scales = "free_x")

data.frame(
  Sample=sample_names(ps),
  Observed = estimate_richness(ps, measures = "Observed"),
  Shannon = exp(estimate_richness(ps, measures = "Shannon")),
  Simpson = 1 / (1 - estimate_richness(ps, measures = "Simpson")),
  Chao1 = estimate_richness(ps, measures = "Chao1"),
  Total_reads = sample_sums(ps)
) %>% 
  left_join(., my_mtd_trans, by=c("Sample" = "ID")) -> df

my_mtd_trans %>% 
  mutate(stage_ind = case_when(
    Stage == "egg" ~ 1,
    Stage == "larva" ~ 2,
    Stage == "pupa" ~ 3,
    Stage == "adult" ~ 4,
    Stage == "mother" ~ 5,
    TRUE ~ NA_integer_
  )) %>% 
  arrange(stage_ind) %>% 
  pull(ID) -> samples_bon_ordre

my_mtd_trans %>% 
  mutate(stage_ind = case_when(Stage == "egg" ~ 1,Stage == "larva" ~ 2,Stage == "pupa" ~ 3,Stage == "adult" ~ 4,Stage == "mother" ~ 5,TRUE ~ NA_integer_)) %>% 
  arrange(stage_ind) %>% 
  distinct(Stage) %>%
  pull(Stage) %>% 
  as.character()-> etiquettes_bon_ordre

ggplot(df, aes(x=Sample, y = Shannon, fill=Stage))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  scale_x_discrete(limits = samples_bon_ordre)

ggplot(df, aes(x=Sample, y = Shannon, fill=factor(Stage, levels=etiquettes_bon_ordre)))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  labs(fill ="Stage")+
  scale_x_discrete(limits = samples_bon_ordre)

ggplot(df, aes(x=Sample, y = Shannon, fill=Stage))+
  geom_col()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))+
  facet_grid(.~reorder(Stage, match(Stage, c("egg", "larva", "pupa", "adult", "mother"))), 
             scales = "free_x")

# ggiNEXT(out, type=1, facet.var="None", color.var = "Assemblage")+
#   # facet_wrap(env, scales="free", ncol=4)+
#   xlab("Number of reads") +
#   ylab("Number of genera") +
#   theme_minimal()+
#   theme(legend.position = "right")+
#   guides(shape = "none") +
#   scale_shape_manual()

# unique_colors <- viridis_pal()(length(unique(out$Data$Assemblage)))
# out$DataInfo %>% View()
# out$AsyEst %>% View()
# out$iNextEst$coverage_based
# 
# ggiNEXT(out, type=1, facet.var="None", color.var = "Assemblage") +
#   xlab("Number of reads") +
#   ylab("Number of genera") +
#   theme_minimal() +
#   theme(legend.position = "right") +
#   scale_color_manual(values = unique_colors)


#b. Diversities------
#Gamma diversity
library(vegan)

otu_table(ps_bin_family) %>% View()

my_ps <- ps
df_taxo = as.data.frame(tax_table(my_ps))
df_taxo$tax_ID = rownames(df_taxo)
df_reads = as.data.frame(otu_table(my_ps))
df_reads$tax_ID = rownames(df_reads)
df_wide = left_join(df_taxo, df_reads, by = "tax_ID") 
BF <- df_wide[, c(9:ncol(df_wide))]
if (length(which(BF$genus == "uncultured")) == 2)
  BF$genus[which(BF$genus == "uncultured")] <-
  c("uncultured_1", "uncultured_2")
BF$genus[which(BF$genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")] <-
  "Rhizobium"
rownames(BF) <- BF$tax_ID

# colnames(BF) <- mtd$label[match(colnames(BF), mtd$sample_ID)]
FB <- t(BF)
Log_FB = log10(FB + 1)

FB <- round(FB, digits = 0)
res_spl<-multipart(y=FB, scales = 1, global = TRUE, relative = FALSE, nsimul=99)
res_spl$statistic

#Alpha diversities by sample (Figure 2)
library(iNEXT)

load("datatransmission/out.Rdata")

div_df<-out$AsyEst[out$AsyEst$Diversity=="Shannon diversity",]

div_df<-left_join(div_df, my_mtd_trans, by = c("Assemblage" = "ID"))

div_df$Family <- as.character(div_df$Family)
div_df$replicate_ind <- as.character(div_df$replicate_ind)

div_df %>% 
  mutate(Stage = factor(Stage, levels = c("egg", "larva", "pupa", "adult", "mother"),
                        labels = c("Oeufs",
                                   "Larves", "Nymphes", "Adultes", "Mères"))) -> div_df

plot_heatmap(ps_rare, "RDA", "none", sample.label="Stage", taxa.label="family")

div_df %>%  View()

library(lme4)
library(nlme)

ggplot(div_df, aes(x=Estimator, fill=Family))+
  geom_histogram(binwidth = 0.3)+
  theme_minimal()+
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#F0E442", "#00BA38", "#619CFF", "#D55E00", "#CC79A7", "#999999"))

ggplot(div_df, aes(x= Stage, y=Estimator, fill=Stage))+
  geom_boxplot()+
  theme_minimal()+
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#F0E442", "#00BA38", "#619CFF", "#D55E00", "#CC79A7", "#999999"))+
  labs(x = "Stade de développement", y= "Indice de Shannon")

ggplot(div_df, aes(x= Stage, y=Estimator, fill=Stage))+
  geom_histogram(binwidth = 0.01)+
  theme_minimal()+
  scale_fill_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#F0E442", "#00BA38", "#619CFF", "#D55E00", "#CC79A7", "#999999"))+
  labs(x = "Stade de développement", y= "Indice de Shannon")

mod <- glm(Estimator ~ Stage*culture*replicate_ind*protein*replicate_culture*Family, data=div_df, family = "quasigaussian")
mod <- glm(Estimator ~  Stage, data=div_df, family = Gamma(link = "inverse"))
summary(mod)
anova(mod)

mod <- lmer(Estimator ~ Stage + (1|Family), data=div_df)
summary(mod)
anova(mod)

mod <- glmer(Estimator ~ 0 + Stage + (1|Family), data=div_df, family = Gamma(link = "inverse"))
summary(mod)
library(multcomp)
glht(mod, linfct = mcp(Stage = "Tukey"))

# Modèle GLMER avec l'effet aléatoire de Family 
## On a testé le glmer mais ça a pas fonctionné 
mod <- glmer(Estimator ~ Stage*replicate_ind*replicate_culture + (1|Family), data=div_df, family = poisson)
# Modèle LMER avec l'effet aléatoire de Family
mod <- lmer(Estimator ~ Stage*protein + (1|Family), data=div_df)
aov1 <- aov(Estimator ~ Stage+replicate_ind+replicate_culture+Family, data = div_df)
summary(aov1)
#des variables semblent trop corrélées et sont enlevées 
cor(, method = "spearman")
mod
anova(mod)
summary(mod)

#ATTention je mets ça là
div_df %>% 
  ggplot(., aes(x=Estimator))+
  geom_histogram()
div_df %>% 
  filter(!Stage=="mother")-> div_petit

mod <- glmer(Estimator ~ Stage + (1|Family), data=div_petit, family="inverse.gaussian")
warnings(mod)

mod <- glmer(Estimator ~ Stage + (1|Family), data = div_df, family = poisson(link = "log"), 
             glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


modlm <- lm(Estimator ~ Stage*Family*replicate_ind*replicate_culture, data=div_df)
anova(modlm)
summary(modlm)

modror <- aov(Estimator~Stage+protein+Error(Family),data=div_df)
summary(modror)

library(lmtest)
model.res <- residuals(modlm)
dwtest(modlm)
shapiro.test(model.res)
bptest(modlm)

modlm2 <- aov(Estimator ~ Stage*Family*replicate_ind*replicate_culture, data=div_df)
anova(modlm2)
#spécifié le tyoe d'anova (1, 2 ou 3)

library(lmtest)
model.res1 <- residuals(modlm2)
dwtest(modlm2)
shapiro.test(model.res1)
bptest(modlm2)

#aov et lm ont meme resultats partout

mod1 <- aov(Estimator ~ Stage*culture*replicate_ind*protein*replicate_culture*Family, data = div_df)
anova(mod1)

mod2 <- aov(Estimator ~ Stage, data = div_df)
plot(TukeyHSD(mod2))

#Bar plot Stage
div_df %>% 
  arrange(Stage) %>% 
  pull(Assemblage) -> ordre

ggplot(div_df, aes(
  x = Assemblage,
  y = Estimator,
  fill = Stage
)) +
  geom_col() +
  geom_errorbar(aes(ymin=LCL, ymax=UCL), width=.5,
                position=position_dodge(.9)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(face = "italic")
  ) +
  scale_x_discrete(name = "Echantillons environnementaux", limits = ordre) +
  scale_y_continuous(name = "Diversité alpha") +
  scale_fill_manual(
    name = "Stade de développement",
    values = c("#F15F79", "#FFB347", "#77DD77", "#AEC6CF", "#F3A0F2")
  )
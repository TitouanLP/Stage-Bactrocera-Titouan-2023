#*******************#
#### NGB - TRANSMISSION ####
#*******************#

#Packages
library(tidyverse)
library(phyloseq)
library(dplyr)

#Load data----
#**Nouvelle taxo----
silva_data = "databioinfo/tax_slv_ssu_138.1.txt.gz"
illumina_data_trans = "datatransmission/Transmision_MiSeq_16S_20230307_864_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2.gz"

#Functions----
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

#Issues
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
good_taxonomic_levels <- c("domain", "phylum", "class", "order", "family", "genus")

#*******************----
#Basic statistics----
#Nombre de séquences : 
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% dim() #1167 séquences sur RUN 2 et RUN 4

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(c(806,807)) %>% 
  filter(.[,1]> 0) %>% dim () #33 séquences dans le tube MOCK

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(c(806,807)) %>% 
  filter(.[,2]> 0) %>% dim () #44 séquences dans le tube MOCK

#Observer la distribution des abondances dans les mocks 
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select("OTU", "total", "length", "spread", "quality", "identity","taxonomy","CTL-DNA-Moc-A_S377","CTL-DNA-Moc-B_S761","references") %>% 
  arrange(desc(.[,8])) %>% view() 
#Il semble que les OTUS de la mocks sont les numéros 2, 3, 4, 39, 62, 66, 77, 94, 105. On les trouve bel et bien avec des nombres de reads >1E3 . Mais l'assignation ne trouve pas de nom fiable pour tous 

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select("CTL-DNA-Moc-A_S377","CTL-DNA-Moc-B_S761") %>% 
  mutate_all(., ~replace(., is.na(.), 0)) %>% 
  summarise_all(.,sum) #22000 et 16600 reads

#Nombre de séquences non assignées à quoi que ce soit : 80
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%   
  filter(is.na(references)) %>% 
  mutate(.before = 1, total = rowSums(.[,14:ncol(.)],na.rm=T)) %>% View() 
#OTU 9 et 10 avec plus de 150 000 reads sont non assignés. Ce sont les OTU avec les 9 et 10 eme nombre de reads. Le 8 eme étant à 272 000. Copier séquences et blaster
#OTU 9 : Apodemus sylvaticus (Mulot sylvstre)
#OTU 10: Myodes glareolus (Campagnole roussâtre)

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%
  filter(!is.na(references)) %>%   
  filter(str_detect(taxonomy, "Mitochondria")) %>% 
  view() #ici on tombe sur des mitochondries de champignons
#SOmmer colonne total pour donner info + valeur abondance comparer à nombre de reads totals 

#**Mock expected----
#On charge la Mock théorique et on l'arrange (attention il faut corriger le fichier à la main pour deux ou trois noms qui ont changé dans silva...)
read.csv2(file = "datatransmission/mock_zymo_expected.csv") %>% 
  mutate(Rel_ab = as.numeric(.$Rel_ab)) %>% 
  rename(genus = Genus, family = Family, order = Order, class = Class, phylum = Phylum) %>% 
  mutate(domain = "Bacteria", 
         OTU = NA, 
         Sample = "Expected", 
         Abundance = Rel_ab) %>% 
  relocate(domain, .before = phylum) %>% 
  select(-Species) -> Mock.expected

Mock.expected %>% view()

#**Metadata----
#Formatage de my_mtd pour phyloseq
# read.csv2("datatransmission/metadata_illumina.csv") %>% 
#   add_count(Label) %>%
#   mutate(Label = paste0(str_pad(row_number(),width=2, side="left",pad = "0"), "_", Label)) %>%
#   select(-n) %>% 
#   filter(type=="ZYMO" | type =="NEG") -> my_mtd
read.csv2("datatransmission/mtd_transmission_cor.csv") %>% 
  mutate(row = ID) %>% 
  column_to_rownames(., var = "row") %>% 
  filter(type == "ZYMO") -> my_mtd_mock   

my_mtd_mock %>%
  select(ID) %>%
  pull(ID) -> samples_to_keep 

#**OTU pour tax ancetre commun ----
columns_to_keep <- c("OTU", "references", samples_to_keep)

read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  filter(!is.na(references)) %>% 
  select(all_of(columns_to_keep)) %>% 
  replace(is.na(.),0) %>% 
  mutate(.before = 1,total = rowSums(.[,3:ncol(.)])) %>% 
  filter(total > 0) %>% 
  select(-total,-references) %>% 
  column_to_rownames(var = "OTU") -> my_otu

#**Tax ancetre commun (sans séparer les ref)----
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(OTU,taxonomy) %>% # rajouter les OTU 9 et 10
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

otu_table(my_otu,taxa_are_rows = TRUE) -> OTU

tax_table(as.matrix(my_tax)) -> TAX

sample_data(my_mtd_mock) -> MTD

ps_full = phyloseq(OTU, TAX, MTD)

# sample_names(OTU) == sample_names(MTD)

ps_zymo = prune_samples(sample_data(ps_full)$type == "ZYMO", ps_full)
ps_zymo = prune_taxa(taxa_sums(ps_zymo)>0, ps_zymo)

#**Visualisation----
ps_object = ps_zymo
as.data.frame(tax_table(ps_object)) %>%  #get taxonomy
  rownames_to_column(var = "tax_ID") -> df 
as.data.frame(otu_table(ps_object)) %>% #get reads 
  rownames_to_column(var = "tax_ID") %>% 
  left_join(df, ., by="tax_ID") %>% 
  select(-tax_ID) -> df

names(df)

df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>%
  filter(total>0) %>%
  select(-total) -> df

df %>% view()

#avec des NA
df %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>%
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>%
  rbind.data.frame(., Mock.expected) -> df_with_NA

ggplot(df_with_NA, aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i){
    paste(domain[i], phylum[i],class[i], order[i], family[i], genus[i], OTU[i], sep="_")})) %>% 
  group_by(tax_ID) %>% 
  mutate(Rel_ab_new = sum(Rel_ab)) %>% 
  pivot_wider(id_cols = 1:6, names_from = Sample, values_from = Rel_ab_new, values_fill = NA) %>% view()

#sans NA
df %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>%
  filter(!is.na(family)) %>% 
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>%
  rbind.data.frame(., Mock.expected) %>% 
  ggplot(., aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

#**OTU distribuées----
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>%
  filter(!is.na(references)) %>%
  select(OTU, references) -> references_column

my_otu %>% 
  mutate(OTU = as.numeric(rownames(.))) %>% 
  left_join(y=references_column, by = "OTU") %>% 
  pivot_longer(cols = -c(OTU, references), names_to = "sample_code", values_to = "reads") %>% 
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
  column_to_rownames(var = "code_unique") -> my_otu_distrib

#**Tax distribuée----
read_tsv(illumina_data_trans,
         na = c("", "NA", "No_hit", "0", "0.0"),
         show_col_types = FALSE) %>% 
  select(OTU, taxonomy, references) %>% 
  filter(OTU %in% rownames(my_otu)) %>%
  separate_rows(references, sep = ",") %>% 
  left_join(read_tsv(file = "databioinfo/SILVA_138.1_SSURef_tax_silva.acc2taxo.gz"), by=c("references" = "accession")) %>% 
  select(-taxonomy.x) %>%  
  separate_rows(taxonomy.y, sep = "\\|") %>%
  left_join(x = .,
            y = taxa2ranks,
            by = c("taxonomy.y" = "taxonomic_path")) %>% 
  filter(taxonomic_rank %in% good_taxonomic_levels) %>%  # Enlève les NA (possiblement au niveau de l'espèce)
  mutate(taxonomic_rank = fct_relevel(taxonomic_rank, good_taxonomic_levels)) %>% 
  pivot_wider(names_from = taxonomic_rank, values_from = taxonomy.y) %>% 
  mutate(code_unique = paste(OTU, references, sep="_")) %>%  # On ajoute code unique avec OTU_références
  select(-references) %>% 
  column_to_rownames(var = "code_unique") -> my_tax_distrib

otu_table(my_otu_distrib,taxa_are_rows = TRUE) -> OTU

tax_table(as.matrix(my_tax_distrib)) -> TAX

ps_full_distrib = phyloseq(OTU, TAX, MTD)

ps_zymo = prune_samples(sample_data(ps_full_distrib)$type == "ZYMO", ps_full_distrib)
ps_zymo = prune_taxa(taxa_sums(ps_zymo)>0, ps_zymo) 

#**Visualisation----
ps_object = ps_zymo
as.data.frame(tax_table(ps_object)) %>%  #get taxonomy
  rownames_to_column(var = "tax_ID") -> df 
as.data.frame(otu_table(ps_object)) %>% #get reads 
  rownames_to_column(var = "tax_ID") %>% 
  left_join(df, ., by="tax_ID") %>% 
  select(-tax_ID) -> df

names(df)

df %>% 
  mutate(total = rowSums(.[,8:ncol(.)])) %>% 
  filter(total>0) %>%
  select(-total) -> df

df %>% view()

#avec des NA
df %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>%
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  rbind.data.frame(., Mock.expected) -> df_with_NA

library(paletteer)
ncol = 2
col_tax <- paletteer_d("palettesForR::Windows")[seq(ncol, ncol + 23, 1)]

colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#FFFAC8",
            "#46F0F0", "#FABEBE", "#AAFFC3", "#008080", "#0082C8",
            "#FFD8B1", "#9A6324", "#800000",  "#808000", "#000075", "#FF6600")

df %>% 
  pivot_longer(., cols = 8:ncol(.), names_to = "Sample", values_to = "Abundance") %>%
  group_by(Sample) %>% 
  mutate(Rel_ab = Abundance/sum(Abundance)) %>% 
  ungroup() %>% 
  filter(family %in% c("Enterobacteriaceae", "Lactobacillaceae", "Bacillaceae", "Staphylococcaceae", "Listeriaceae", "Enterococcaceae", "Pseudomonadaceae", "Hungateiclostridiaceae", "Erwiniaceae")) %>% 
  rbind.data.frame(., Mock.expected) %>% 
  ggplot(., aes(x=Sample, y = Rel_ab, fill = family)) +
  geom_bar(stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = colors) +
  theme(legend.position="right",
        axis.text.x = element_text(angle=90,vjust=0.4,hjust=1.2))

#Afficher les abondances par taxon
#Niveau genre
df_with_NA %>%
  group_by(genus, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>%  
  pivot_wider(id_cols = 1, names_from = Sample, values_from = Rel_ab, values_fill = NA) %>%
  view()
#1 enterobacter à 10-2, Klebselia idem et limocilacobacillus à 10-1. 

df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i){
    paste(domain[i], phylum[i],class[i], order[i], family[i], genus[i], sep="_") 
  })) %>% 
  group_by(tax_ID, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>%  
  pivot_wider(id_cols = 1, names_from = Sample, values_from = Rel_ab, values_fill = NA) %>%
  view()
#Faux positifs est l'OTU le plus abondant : Limosilactobacillus. Il n'est pas assigné au bon genre 
#Trop de mauvaises assignations au niveau genre -> Il faut travailler au niveau famille. 

#Niveau famille
df_with_NA %>%
  mutate(tax_ID = sapply(1:nrow(.), function(i){
    paste(domain[i], phylum[i],class[i], order[i], family[i], sep="_") 
  })) %>% 
  group_by(tax_ID, Sample) %>% 
  summarise(Rel_ab = sum(Rel_ab)) %>%  
  pivot_wider(id_cols = 1, names_from = Sample, values_from = Rel_ab, values_fill = NA) %>%
  view() 
#La famille fausse positive la plus abondante est à 1.2E-3 à 1.5E-3. Sinon c'est parfait !

#Niveau OTU (tableau montrant les confusions des OTUs de la Mock avec d'autres genres/familles)
df_with_NA %>%
  mutate(OTU = str_remove_all(OTU, " ")) %>% 
  filter(OTU %in% c()) %>%
  pivot_wider(id_cols = 1:7, names_from = Sample, values_from = Rel_ab, values_fill = NA, values_fn = sum) %>%
  View()

library(tidyverse)


download.file(url = "https://github.com/carpentries-incubator/bioc-intro/raw/main/episodes/data/rnaseq.csv",destfile = "data/rnaseq.csv")

#Reading in a csv file using readr
rna <- read_csv("data/rnaseq.csv")
rna

#using select() function to select columns from rna data frame
select(rna, gene, sample, tissue, expression)
select(rna, -tissue, -organism)

#using filter() to select specific rows of the rna data frame
filter(rna, sex == "Male")
filter(rna, sex == "Male" & infection == "NonInfected")

#Selected columns representing mouse genes and their human homolog
genes <- select(rna, gene, hsapiens_homolog_associated_gene_name)

#Filtering mouse genes that have no human homolog
filter(genes, is.na(hsapiens_homolog_associated_gene_name))

#Filtering mouse genes that only have a human homolog
filter(genes, !is.na(hsapiens_homolog_associated_gene_name))

#Using pipes %>% (ctrl + shift + m) to manipulate the rna data frame
rna %>% 
  filter(sex == "Male") %>% 
  select(gene, sample, tissue, expression)

#Creating a new data frame from rna
rna2 <- rna %>% 
  filter(sex == "Male") %>% 
  select(gene, sample, tissue, expression)

rna2

rna3 <- rna %>% 
  filter(time == 0 & sex == "Female" & expression > 50000) %>% 
  select(gene, sample, time, expression, age)

rna3


#Using the mutate() function to create a new column
rna %>% 
  mutate(time_hours = time * 24, time_min = time_hours *60) %>% 
  select(time,time_hours, time_min)


rna4 <- rna %>% 
  mutate(expression_log_transformed = log(expression)) %>% 
  filter(chromosome_name == "X" | chromosome_name == "Y") %>% 
  filter(!is.na(phenotype_description)) %>% 
  filter(expression_log_transformed > 5) %>% 
  select(gene, chromosome_name, phenotype_description, sample, expression_log_transformed)
  
rna4

#Grouping rna by gene
rna %>% 
  group_by(gene)

#Grouping rna by sample
rna %>% 
  group_by(sample)

#Utilizing summary() function with group_by() to find the mean expression level for each gene
rna %>% 
  group_by(gene) %>% 
  summarise(mean_expression = mean(expression))

#Finding the mean expression level of all genes in each sample
rna %>% 
  group_by(sample) %>% 
  summarise(mean_expression = mean(expression))

#Finding the mean expression for the time and infection level groups for each gene
rna %>% 
  group_by(gene, infection, time) %>% 
  summarise(mean_expression = mean(expression))

rna %>% 
  group_by(gene, infection, time) %>% 
  summarise(mean_expression = mean(expression), median_expression = median(expression))

#Mean expression level of gene “Dok3” by timepoints.
rna %>% 
  filter(gene == "Dok3") %>% 
  group_by(time) %>% 
  summarise(mean_expression = mean(expression))

#Counting the number of observations for each infection status
rna %>% 
  count(infection)
  
#Arranging a count table of infection status and time by time
rna %>% 
  count(infection, time) %>% 
  arrange(time)

#Arranging a count table of infection status and time by counts
rna %>% 
  count(infection, time) %>% 
  arrange(n)

#Arranging a count table of infection status and time by counts in descending order
rna %>% 
  count(infection, time) %>% 
  arrange(desc(n))

#1. Number of genes in each sample
rna %>% 
  count(sample)

#2.Sequencing depth for each sample 
rna %>% 
  group_by(sample) %>% 
  summarise(sequencing_depth = sum(expression)) %>% 
  arrange(desc(sequencing_depth))

#3. Selecting a sample and evaluating the number of genes by biotype
rna %>% 
  filter(sample == "GSM2545350") %>% 
  count(gene_biotype) %>% 
  arrange(desc(n))

#4.Identifying genes associated with the “abnormal DNA methylation” phenotype description, 
#and calculating their mean expression (in log) at time 0, time 4 and time 8.
rna %>% 
  filter(phenotype_description == "abnormal DNA methylation") %>% 
  mutate(log_expression = log(expression)) %>% 
  group_by(gene, time) %>% 
  summarise(mean_log_expression = mean(log_expression))

#Creating a new data frame from rna containing gene, sample, and expression columns
rna_exp <- rna %>%
  select(gene, sample, expression)

rna_exp


#Creating a new data frame of rna_exp in wide format
rna_wide <- rna_exp %>% 
  pivot_wider(names_from = sample, values_from = expression)

rna_wide

#Creating a new data frame for an example. %in% for filtering for multiple values
#instead of gene == x | gene == y | gene == z
rna_with_missing_values <- rna %>% 
  select(gene, sample, expression) %>% 
  filter(gene %in% c("Asl", "Apod", "Cyp2d22")) %>% 
  filter(sample %in% c("GSM2545336", "GSM2545337", "GSM2545338")) %>% 
  arrange(sample) %>% 
  filter(!(gene == "Cyp2d22" & sample != "GSM2545338"))

rna_with_missing_values

#converting rna_with_missing_values data frame into wider format. Filling missing
#information with 0 instead of "NA"
rna_with_missing_values %>% 
  pivot_wider(names_from = sample, values_from = expression, values_fill = 0)

#Convert rna_wide into long format. -gene in col argument so use all columns except gene
rna_long <- rna_wide %>% 
  pivot_longer(names_to = "sample", values_to = "expression", -gene)

rna_long

rna_wide %>% 
  pivot_longer(names_to = "sample", values_to = "expression", cols = starts_with("GSM"))

##1
rna_mouse <- rna %>% 
  select(gene, expression, mouse) %>% 
  pivot_wider(names_from = mouse, values_from = expression) %>% 
  pivot_longer(names_to = "mouse#", values_to = "expression", cols = -gene)

##2.
rna %>% 
  select(sex, chromosome_name, expression) %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  group_by(sex, chromosome_name) %>% 
  summarise(mean_expression = mean(expression)) %>% 
  pivot_wider(names_from = sex, values_from = mean_expression)

##3.
rna %>% 
  group_by(gene, time) %>% 
  summarise(mean_expression = mean(expression)) %>% 
  pivot_wider(names_from = time, values_from = mean_expression)

##4.
rna %>% 
  group_by(gene, time) %>% 
  summarise(mean_expression = mean(expression)) %>% 
  pivot_wider(names_from = time, values_from = mean_expression) %>% 
  mutate(time_8_vs_0 = `8`/`0`, time_8_vs_4 = `8` /`4`) %>% 
  pivot_longer(names_to = "Time comparisons", values_to = "Fold Change", cols = starts_with("time"))

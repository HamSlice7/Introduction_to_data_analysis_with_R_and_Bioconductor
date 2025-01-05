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

  
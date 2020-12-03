library(ggplot2)
#install.packages("refGenome")
library(refGenome)
#install.packages("wesanderson")
library(wesanderson)

library(dplyr)
#library(tidyr)
#library(data.table)
#library(rlist)
#library(stringr)
main_dir="/home/jl/Dropbox/university/UPF_msc/PGB_2020/project_summary"
raw_dir=paste0(main_dir,"/PGB2020_raw_files")

summary_table <- read.table(paste0(main_dir,'/results/basic_statistics.tsv'), sep='\t', header=TRUE)

summary_table$binomial <- paste0(summary_table$genus,' ',summary_table$species)

ggplot(data=summary_table, aes(x=reorder(binomial,num_transcripts), y=num_transcripts, fill=dataset)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=num_transcripts), vjust=0.5, hjust=-0.1, color="black",
            position = position_dodge(0.9), size=3.5) +
  scale_fill_brewer(palette="Paired") +
  theme_minimal() +
  xlab('Number of reconstructed transcripts')+
  ylab('Species')+
  labs(fill='Dataset') +
  coord_flip()

ggplot(data=summary_table, aes(x=reorder(binomial,total_transcript_length), y=total_transcript_length/1000000, fill=dataset)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=paste0(total_transcript_length/1000000,' Mb')), vjust=0.5, hjust=-0.1, color="black",
              position = position_dodge(0.9), size=3.5) +
    scale_fill_brewer(palette="Set1") +
    theme_minimal() +
    xlab('Number of reconstructed transcripts')+
    ylab('Species')+
    labs(fill='Dataset') +
    coord_flip()

##Construct our big dataset from student's independent files (known and novel gff)
raw_data <- data.frame() #Start empty
setwd(raw_dir)

##Known
#Search all files that meet the pattern "know"
file.list<- list.files(path=raw_dir, pattern='*know*', recursive = TRUE, full.names = FALSE)


#Concatenate all datasets together. File name as column
for (f in file.list) {
    dataset <- ensemblGenome()
    print(f)
    read.gtf(dataset,f)
    df <- as.data.frame(dataset@ev$gtf)
    #set filename and dataset
    df$fname <- f
    df$dataset <- 'known'

    #Exception for kudria, error in format (stringtie version?)
    if (grepl('kudria',f, fixed=TRUE)) {
    df$ref_gene_id <- 'NA' }
    raw_data <- rbind(raw_data, df)
}

##Novel
#Search all files that meet the pattern "novel"
file.list<- list.files(path=raw_dir, pattern='*nove*', recursive = TRUE, full.names = FALSE)

#Concatenate all datasets together. File name as column
for (f in file.list) {
    dataset <- ensemblGenome()
    print(f)
    read.gtf(dataset,f)
    df <- as.data.frame(dataset@ev$gtf)
    #set filename and dataset
    df$fname <- f
    df$dataset <- 'novel'
    #Add missing fields in novel
    df$ref_gene_id <- 'NA'
    df$ref_gene_name <- 'NA'
    df$reference_id <- 'NA'#}
    raw_data <- rbind(raw_data, df)
}

#Convert variables to numeric
raw_data[, c(10,16,17)] <- sapply(raw_data[, c(10,16,17)], as.numeric)

#Get info from species by joining the summary dataset
fname2species <- summary_table %>%
select(binomial, fname)

data <- left_join(raw_data,fname2species, by="fname")


ggplot(data, aes(log(TPM), fill = dataset, colour = dataset)) +
  geom_density(alpha = 0.1)+
  scale_fill_manual(values = wes_palette("Royal1"))+
  scale_colour_manual(values = wes_palette("Royal1"))+
  theme_minimal() +
      xlab('log(TPM)')+
      ylab('Density')+
      labs(colour='Dataset')+
      guides(fill=FALSE)

ggplot(data, aes(log(TPM), fill = dataset, colour = dataset)) +
  geom_density(alpha = 0.1)+
  scale_fill_manual(values = wes_palette("Royal1"))+
  scale_colour_manual(values = wes_palette("Royal1"))+
  theme_minimal() +
      xlab('log(TPM)')+
      ylab('Density')+
      labs(colour='Dataset')+
      guides(fill=FALSE)+
      facet_wrap(~binomial)

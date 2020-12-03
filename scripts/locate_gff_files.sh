#!/bin/bash
main_dir=/home/jl/Dropbox/university/UPF_msc/PGB_2020/project_summary
mkdir $main_dir/results

#Create basic statistics tab separated empty file (only header)
echo -e 'genus\tspecies\tdataset\tnum_transcripts\ttotal_transcript_length\tnum_chr\tfname'
echo -e 'genus\tspecies\tdataset\tnum_transcripts\ttotal_transcript_length\tnum_chr\tfname' > $main_dir/results/basic_statistics.tsv

#Read file with species names
input=${main_dir}/metadata/species_sorted.txt
datasets='known novel'
while IFS= read -r line
do
  for dataset in $datasets; do
    binomial=(${line// / }) #Use space as separator to get information in an array
    genus=${binomial[0]} #Position one is genus
    species=${binomial[1]} #Position two is species
    #Find files based on species name and dataset (minus last position)
    file=`find $main_dir/PGB2020_raw_files/ -iname "*$species*" -type f | grep ${dataset::4}`

    #Obtain total length based on substraction of start and end positions for each transcript
    total_length=`cat $file | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$5-$4 }END{print SUM}'`

    #Get number of transcripts per file
    num_transcripts=`wc -l $file | cut -f1 -d ' '`

    #Get number of contigs ~ chr
    num_chr=`cat $file | cut -f1 | sort | uniq | wc -l`

    #File without directory
    onlyfile=`echo $file | sed 's:.*/::'`

    #Print results to screen and file
    echo -e $genus'\t'$species'\t'$dataset'\t'$num_transcripts'\t'$total_length'\t'$num_chr'\t'$onlyfile
    echo -e $genus'\t'$species'\t'$dataset'\t'$num_transcripts'\t'$total_length'\t'$num_chr'\t'$onlyfile >> $main_dir/results/basic_statistics.tsv
  done
done < "$input"

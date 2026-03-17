#!/bin/bash

#Get the list of fasta files

while getopts ":i:" opt; do
  case $opt in
    i)
      infile="$OPTARG"
      ;;
  esac
done

wdir=$(pwd)
prefix=$(dirname $infile) #prefix to where the file is located
filename=$(basename $prefix) #name of the sample

#Create a file with all spacers that come from repeat regions of >9 repeats
if [ ! -e "./$prefix/filtered_spacers.fasta" ]; then
    while read -r line;
    do
        cat ./$prefix/$line.fa >> ./$prefix/filtered_spacers.fasta
    done < $infile
fi

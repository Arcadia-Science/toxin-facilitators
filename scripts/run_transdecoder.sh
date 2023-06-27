#! /usr/bin/bash

# not going to spend the time making this a better command-line script/workflow because not sure if this is even going to work

for file in *.fasta; do
    name=$(basename $file .fasta);
    pep=${name}.fasta.transdecoder.pep;
    newpep=${name}.transdecoder.pep.fasta;
    mkdir $name;
    mv $file $name/;
    cd $name/;
    TransDecoder.LongOrfs -t $file;
    TransDecoder.Predict -t $file;
    mv $pep $newpep;
    cd ../
done

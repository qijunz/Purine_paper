#!/bin/bash

#unpack RSEM, Bowtie
tar -xzvf RSEM-1.3.1.tar.gz
unzip bowtie2-2.3.4-linux-x86_64.zip

# set PATH
export PATH=$(pwd)/bowtie2-2.3.4-linux-x86_64:$PATH
export PATH=$(pwd)/RSEM-1.3.1/bin:$PATH
export PATH=$(pwd)/RSEM-1.3.1/bin/samtools-1.3:$PATH

#transfer and unpack reads
cp /staging/qzhang333/Ath_Trans_reads/BHHTKYBCXX.$1.R1.paired.cat.qualtrimmed.mouseRemoved.fq .
cp /staging/qzhang333/Ath_Trans_reads/BHHTKYBCXX.$1.R2.paired.cat.qualtrimmed.mouseRemoved.fq .

# transfer RSEM bowtie2 indexed DO metagenes reference
cp /staging/qzhang333/operon_gene.tar.gz .

tar -xzvf operon_gene.tar.gz
rm *.tar.gz

# make directory for expression calculating
mkdir rsem-UAgene-$1

# calculate expression
rsem-calculate-expression -p 1 \
                          --bowtie2 \
                          --paired-end \
                          --no-bam-output \
                          BHHTKYBCXX.$1.R1.paired.cat.qualtrimmed.mouseRemoved.fq \
                          BHHTKYBCXX.$1.R2.paired.cat.qualtrimmed.mouseRemoved.fq \
                          operon_gene/operon_gene \
                          rsem-UAgene-$1/$1
                    

FILE=rsem-UAgene-$1/$1.genes.results
if [ -f "$FILE" ]; then
    mv rsem-UAgene-$1/$1.genes.results /staging/qzhang333
fi

rm -r operon_gene
rm *.fq

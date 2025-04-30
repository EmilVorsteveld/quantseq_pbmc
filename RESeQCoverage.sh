#!/bin/bash
for fn in ./*.sorted.bam;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"

geneBody_coverage.py -r ./coverage/hg38.HouseKeepingGenes.nochr.bed -i ${samp} -o coverage

mv coverage.geneBodyCoverage.curves.pdf ./coverage/figures/${samp}.pdf
mv coverage.geneBodyCoverage.txt ./coverage/data/${samp}.txt

done

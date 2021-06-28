while read line
do 

cd ${line}
cut -f1,9 ${line}._stringtie_hisat.tab | sed "s/^/${line}\t/" |  sed '1 i\Name\tGene_ID\tTPM' | sed '2d' > ${line}.tpm.txt
cp ${line}.tpm.txt ../
cd ..

done<lista

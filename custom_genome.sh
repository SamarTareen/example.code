

###call this script in the same directory as the genome files
# copyright Samar H. K. Tareen

if true; then
##download the gtf file
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
gunzip Mus_musculus.GRCm38.101.gtf.gz

##download the fasta file
wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

#filter the gtf file
cellranger mkgtf Mus_musculus.GRCm38.101.gtf Mus_musculus.GRCm38.101.filtered.gtf --attribute=gene_biotype:protein_coding
fi

if true; then
#get sequence lengths of the custom sequences to be added
#cat eGFP.fa | grep -v "^>" | tr -d "\n" | wc -c
#cat IL2.fa | grep -v "^>" | tr -d "\n" | wc -c
#cat IL2_subseq.fa | grep -v "^>" | tr -d "\n" | wc -c
#cat IL2_ENSMUSG00000027720.fa | grep -v "^>" | tr -d "\n" | wc -c
cat AAV-GFAP2.2-EGFP.fa | grep -v "^>" | tr -d "\n" | wc -c
cat pAAV-GFAP2.2-mIL2-WPRE.fa | grep -v "^>" | tr -d "\n" | wc -c

#creating gtf entries for the custom sequences
#echo -e 'eGFP\tunknown\texon\t1\t729\t.\t+\t.\tgene_id "eGFP"; transcript_id "eGFP"; gene_name "eGFP"; gene_biotype "protein_coding";' > eGFP.gtf
#echo -e 'IL2\tunknown\texon\t1\t955\t.\t+\t.\tgene_id "IL2"; transcript_id "IL2"; gene_name "IL2"; gene_biotype "protein_coding";' > IL2.gtf
#echo -e 'IL2_subseq\tunknown\texon\t1\t510\t.\t+\t.\tgene_id "IL2_subseq"; transcript_id "IL2_subseq"; gene_name "IL2_subseq"; gene_biotype "protein_coding";' > IL2_subseq.gtf
#echo -e 'IL2_ENSMUSG00000027720\tunknown\texon\t1\t6637\t.\t+\t.\tgene_id "IL2_ENSMUSG00000027720"; transcript_id "IL2_ENSMUSG00000027720"; gene_name "IL2_ENSMUSG00000027720"; gene_biotype "protein_coding";' > IL2_ENSMUSG00000027720.gtf
echo -e 'AAV-GFAP2.2-EGFP\tunknown\texon\t1\t6546\t.\t+\t.\tgene_id "AAV-GFAP2.2-EGFP"; transcript_id "AAV-GFAP2.2-EGFP"; gene_name "AAV-GFAP2.2-EGFP"; gene_biotype "protein_coding";' > AAV-GFAP2.2-EGFP.gtf
echo -e 'pAAV-GFAP2.2-mIL2-WPRE\tunknown\texon\t1\t7019\t.\t+\t.\tgene_id "pAAV-GFAP2.2-mIL2-WPRE"; transcript_id "pAAV-GFAP2.2-mIL2-WPRE"; gene_name "pAAV-GFAP2.2-mIL2-WPRE"; gene_biotype "protein_coding";' > pAAV-GFAP2.2-mIL2-WPRE.gtf

#add custom seqs to fasta
cp Mus_musculus.GRCm38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
sed -i -e '$a\' Mus_musculus.GRCm38.dna.primary_assembly_custom.fa #adding end of line to the end of file
#cat eGFP.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
#cat IL2.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
#cat IL2_subseq.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
#cat IL2_ENSMUSG00000027720.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
cat AAV-GFAP2.2-EGFP.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
sed -i -e '$a\' Mus_musculus.GRCm38.dna.primary_assembly_custom.fa #adding end of line to the end of file
cat pAAV-GFAP2.2-mIL2-WPRE.fa >> Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
#grep ">" Mus_musculus.GRCm38.dna.primary_assembly_custom.fa
grep -c "^>" Mus_musculus.GRCm38.dna.primary_assembly_custom.fa

#add custom gtf entries
cp Mus_musculus.GRCm38.101.filtered.gtf Mus_musculus.GRCm38.101.filtered_custom.gtf
#cat eGFP.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf
#cat IL2.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf
#cat IL2_subseq.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf
#cat IL2_ENSMUSG00000027720.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf
cat AAV-GFAP2.2-EGFP.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf
cat pAAV-GFAP2.2-mIL2-WPRE.gtf >> Mus_musculus.GRCm38.101.filtered_custom.gtf


#make the ref genome
cellranger mkref --genome=Mus_musculus.GRCm38.101.filtered_custom --fasta=Mus_musculus.GRCm38.dna.primary_assembly_custom.fa --genes=Mus_musculus.GRCm38.101.filtered_custom.gtf --nthreads=10 --memgb=20

fi


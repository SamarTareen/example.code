
#command script to submit job to queue using ssub
# copyright Samar H. K. Tareen

printf "Submitting one cellranger count job per sample\n\n"

#load ssub and cellranger module
module load ssub
module load cellranger

#exp 3515
#S1_10X_4mos_03082020,S2_10X_4mos_03082020,S3_10X_4mos_03082020,S4_10X_4mos_03082020,S5_10X_4mos_03082020,S6_10X_4mos_03082020

#exp 3527
#S1_10X_4mos_03082020,S2_10X_4mos_03082020,S3_10X_4mos_03082020,S4_10X_4mos_03082020,S5_10X_4mos_03082020,S6_10X_4mos_03082020


##submit 4 months jobs
#S01
ssub -o output_4m_S01.log -c16 --mem=64GB --email \
cellranger count --id=4m_S01 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S01_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S01_3527 \
--sample='S1_10X_4mos_03082020,S1_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S01\n"
#squeue

#S02
ssub -o output_4m_S02.log -c16 --mem=64GB --email \
cellranger count --id=4m_S02 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S02_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S02_3527 \
--sample='S2_10X_4mos_03082020,S2_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S02\n"
#squeue

#S03
ssub -o output_4m_S03.log -c16 --mem=64GB --email \
cellranger count --id=4m_S03 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S03_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S03_3527 \
--sample='S3_10X_4mos_03082020,S3_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S03\n"
#squeue

#S04
ssub -o output_4m_S04.log -c16 --mem=64GB --email \
cellranger count --id=4m_S04 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S04_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S04_3527 \
--sample='S4_10X_4mos_03082020,S4_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S04\n"
#squeue

#S05
ssub -o output_4m_S05.log -c16 --mem=64GB --email \
cellranger count --id=4m_S05 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S05_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S05_3527 \
--sample='S5_10X_4mos_03082020,S5_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S05\n"
#squeue

#S06
ssub -o output_4m_S06.log -c16 --mem=64GB --email \
cellranger count --id=4m_S06 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S06_3515,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/4m_S06_3527 \
--sample='S6_10X_4mos_03082020,S6_10X_4mos_03082020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 4m_S06\n"
#squeue


#exp 3466
#S1-2Y,S2-2Y,S3-2Y,S4-2Y,S5-2Y,S6-2Y

#exp 3497
#S1_2Y_10X_RNA_28052020,S2_2Y_10X_RNA_28052020,S3_2Y_10X_RNA_28052020,S4_2Y_10X_RNA_28052020,S5_2Y_10X_RNA_28052020,S6_2Y_10X_RNA_28052020


##submit 2 year jobs
#S01
ssub -o output_2y_S01.log -c16 --mem=64GB --email \
cellranger count --id=2y_S01 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S01_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S01_3497 \
--sample='S1-2Y,S1_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S01\n"
#squeue

#S02
ssub -o output_2y_S02.log -c16 --mem=64GB --email \
cellranger count --id=2y_S02 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S02_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S02_3497 \
--sample='S2-2Y,S2_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S02\n"
#squeue

#S03
ssub -o output_2y_S03.log -c16 --mem=64GB --email \
cellranger count --id=2y_S03 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S03_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S03_3497 \
--sample='S3-2Y,S3_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S03\n"
#squeue

#S04
ssub -o output_2y_S04.log -c16 --mem=64GB --email \
cellranger count --id=2y_S04 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S04_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S04_3497 \
--sample='S4-2Y,S4_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S04\n"
#squeue

#S05
ssub -o output_2y_S05.log -c16 --mem=64GB --email \
cellranger count --id=2y_S05 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S05_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S05_3497 \
--sample='S5-2Y,S5_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S05\n"
#squeue

#S06
ssub -o output_2y_S06.log -c16 --mem=64GB --email \
cellranger count --id=2y_S06 \
--fastqs=/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S06_3466,/bi/group/liston/samar/AAV_IL2/00_raw_data/fastqs/2y_S06_3497 \
--sample='S6-2Y,S6_2Y_10X_RNA_28052020' \
--transcriptome=/bi/group/liston/samar/AAV_IL2/00_raw_data/genomes/Mus_musculus.GRCm38.101.filtered_custom \
--localcores=16 \
--localmem=64
#print job queue
printf "Submitted 2y_S06\n"
#squeue


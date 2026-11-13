# Convert each vcf file to bgen format
# Completed: All
ChNum=6
qctool -g chr${ChNum}_hrc_fullsamp.vcf.gz -og BGEN/chr${ChNum}_hrc.bgen

VCF_PATH=/FHS_HRC
ChNum=1
# 8,468 samples
bcftools query -l ${VCF_PATH}/chr${ChNum}_hrc_fullsamp.vcf.gz | wc -l
bcftools query -l ${VCF_PATH}/chr${ChNum}_hrc_fullsamp.vcf.gz > ${VCF_PATH}/hrc_fullsamp_list.txt

bgen_path=/FHS_HRC/BGEN
ChNum=19
qctool -g ${bgen_path}/Full_Samp/chr${ChNum}_hrc.bgen \
    -og ${bgen_path}/QC_Samp/chr${ChNum}_hrc_qc.bgen \
    -incl-samples ${bgen_path}/hrc_fin_samp_list.txt

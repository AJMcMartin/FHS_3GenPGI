conda activate SNIPar_GP2_py38
IBD_PATH=unipar_ibd
BGEN_PATH=~/fhs_dat/FHS_HRC_Subs/hrc_bgen_final
IMP_PED=sniparV2_unipar_plink/Affy550K_QC_Sets/snipar_ibd_ped.txt
OUT_PATH=unipar_hdf5

impute.py \
  --ibd ${IBD_PATH}/chr_@.ibd \
  --bgen ${BGEN_PATH}/chr@_hrc_gp \
  --chr_range 22 \
  --pedigree ${IMP_PED} \
  --out ${OUT_PATH}/chr@_hrc_gt \
  --threads 30

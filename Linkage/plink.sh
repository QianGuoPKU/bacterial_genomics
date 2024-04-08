cd $1
plink --vcf merge_for_plink.vcf --maf 0.05 --geno 0.05 --r2 --ld-window 999999  --ld-window-r2 0 --out plink_res --allow-extra-chr --const-fid



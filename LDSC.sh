# run from ~/container.home/gwas_eqtl_colocalization/scripts
# Get GWAS data
for i in AllOA KneeOA TKR; do echo "snpid chr bp a1 a2 beta pval Ncases Ncontrols N" > /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/$i.for-ldscr; done
for i in AllOA KneeOA TKR; do zcat /storage/hmgu/projects/OA_T2D/data_original/GO.FILTER.GW.$i.FULL.09052019.txt.gz | grep -v POS | awk 'OFS="\t"{print $19,$20,$21,$2,$3,$8,$10,$16,$17,$18}' >> /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/$i.for-ldscr ; done

echo "snpid chr bp a1 a2 beta pval" > /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/T2D.for-ldscr
zcat /storage/hmgu/projects/OA_T2D/data_original/Mahajan.NatGenet2018b.T2D.European.txt.gz | grep -v Pos | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$7,$9}' >> /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/T2D.for-ldscr

# Get eQTL data
for i in HighGradeCartilage LowGradeCartilage Synovium; do echo "snpid a1 a2 beta pval N INFO" > /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/$i.for-ldscr; done
for i in HighGradeCartilage LowGradeCartilage Synovium; do zcat /storage/hmgu/projects/OA_T2D/data_original/FunGen_eQTL_$i.FastQTL_perm_nom_info.ForMSK-KP_16Jan2021.txt.gz | grep -v genotype_id | awk 'OFS="\t"{print $1,$27,$26,$5,$4,$7,$28}' >> /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/$i.for-ldscr ; done

echo "snpid chr bp a1 a2 beta pval N" > /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/PancreaticIslets.for-ldscr
zcat /storage/hmgu/projects/OA_T2D/data_original/InsPIRE_PancreaticIslets_Gene_eQTLs_Nominal_Pvalues.txt.gz | grep -v SNPchr | awl 'OFS="\t" {print $6,$7,$8,$10,$9,14,$13,$4}' >> /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/PancreaticIslets.for-ldscr

# Format data for LDSC software
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/AllOA.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/AllOA.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist # --N-cas 68621 --N-con 247846
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/KneeOA.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/KneeOA.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist # --N-cas 19501 --N-con 45354
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/TKR.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/TKR.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist # --N-cas 9062 --N-con 36515
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/T2D.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/T2D.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist --N 824000 --N-cas 74124
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/HighGradeCartilage.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/HighGradeCartilage.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/LowGradeCartilage.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/LowGradeCartilage.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/Synovium.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/Synovium.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist
/storage/hmgu/software/ldsc/munge_sumstats.py --sumstats /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/PancreaticIslets.for-ldscr --out /storage/hmgu/projects/OA_T2D/processed_data/TraitCorr/PancreaticIslets.ldsc --merge-alleles  ~/gwas_eqtl_colocalization/scripts/w_hm3.snplist

# Run LDSC software
# GO only
# for i in AllOA KneeOA TKR; do /storage/hmgu/software/ldsc/ldsc.py --rg $i.ldsc.sumstats.gz,AllOA.ldsc.sumstats.gz,KneeOA.ldsc.sumstats.gz,TKR.ldsc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out $i.correlation --pop-prev 0.0873,0.0873,0.0873 --samp-prev 0.20,0.20,0.20; done

# GO + T2D
for i in AllOA KneeOA TKR T2D; do /storage/hmgu/software/ldsc/ldsc.py --rg $i.ldsc.sumstats.gz,AllOA.ldsc.sumstats.gz,KneeOA.ldsc.sumstats.gz,TKR.ldsc.sumstats.gz,T2D.ldsc.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out $i.correlation --pop-prev 0.0873,0.0873,0.0873,0.1 --samp-prev 0.20,0.20,0.20,0.089956; done

# GWAS + eQTL
# for i in AllOA KneeOA TKR T2D HighGradeCartilage LowGradeCartilage Synovium PancreaticIslets; do /storage/hmgu/software/ldsc/ldsc.py --rg $i.ldsc.sumstats.gz,AllOA.ldsc.sumstats.gz,KneeOA.ldsc.sumstats.gz,TKR.ldsc.sumstats.gz,T2D.ldsc.sumstats.gz, HighGradeCartilage.ldsc.sumstats.gz, LowGradeCartilage.ldsc.sumstats.gz, Synovium.ldsc.sumstats.gz, PancreaticIslets.ldsc.sumstats.gz  --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out $i.correlation --pop-prev 0.0873,0.0873,0.0873,0.1 --samp-prev 0.20,0.20,0.20,0.089956; done

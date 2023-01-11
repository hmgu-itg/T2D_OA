#!/bin/bash
DATA_DIR=/lustre/groups/itg/teams/zeggini/projects/OA_T2D
# DATA_DIR=/project_data

# Format data for LDSC software
#for i in AllOA KneeOA HipOA TKR THR T2D; do
#    if [ $i == "T2D" ]; then
#        echo "snpid chr bp a1 a2 beta pval" > ${DATA_DIR}/LDscore/$i.for.ldscr
#        zcat ${DATA_DIR}/data/Mahajan.NatGenet2018b.T2D.European.txt.gz | grep -v Pos | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$7,$9}' >> ${DATA_DIR}/LDscore/T2D.for.ldscr
#        /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/LDscore/$i.rsid.for.ldscr --out ${DATA_DIR}/LDscore/$i.ldsc --merge-alleles  ${DATA_DIR}/LDscore/w_hm3.snplist --N 824000 --N-cas 74124
#    else
#        echo "snpid chr bp a1 a2 beta pval Ncases Ncontrols N" > ${DATA_DIR}/LDscore/$i.for.ldscr
#        zcat ${DATA_DIR}/data/GO.FILTER.GW.$i.FULL.09052019.txt.gz | grep -v POS | awk 'OFS="\t"{print $19,$20,$21,$2,$3,$8,$10,$16,$17,$18}' >> ${DATA_DIR}/LDscore/$i.for.ldscr
#        /lustre/groups/itg/shared/software/ldsc/munge_sumstats.py --sumstats ${DATA_DIR}/LDscore/$i.for.ldscr --out ${DATA_DIR}/LDscore/$i.ldsc --merge-alleles  ${DATA_DIR}/LDscore/w_hm3.snplist
#    fi
#done

# Run LDSC software
# DATA_DIR=/lustre/groups/itg/teams/zeggini/projects/OA_T2D
traits=(AllOA KneeHipOA TJR) #AllOA KneeOA HipOA TKR THR)
sample_prev=(0.27 0.224 0.125) #0.27 0.187 0.115 0.078 0.078)

#for i in {0..4}; do
#    /lustre/groups/itg/shared/software/ldsc/ldsc.py --rg ${traits[$i]}.ldsc.sumstats.gz,T2D.ldsc.sumstats.gz --ref-ld-chr ${DATA_DIR}/LDscore/eur_w_ld_chr/ --w-ld-chr ${DATA_DIR}/LDscore/eur_w_ld_chr/ --out ${DATA_DIR}/LDscore/$i.correlation --pop-prev 0.0873,0.1 --samp-prev ${sample_prev[$i]},0.089956
#done

my.func () {
    local t=$1
    local s=$2
    local n=$3

    # sample snpid column
    # Rscript --vanilla ${DATA_DIR}/LDscore/permutation_LDscore.R ${DATA_DIR}/LDscore/T2D.ldsc.sumstats.gz ${DATA_DIR}/LDscore/permuted.T2D.ldsc.sumstats.gz
    Rscript --vanilla ${DATA_DIR}/LDscore/permutation_LDscore.R ${DATA_DIR}/LDscore/${t}.ldsc.sumstats.gz ${DATA_DIR}/LDscore/permuted.${t}.${n}.ldsc.sumstats.gz

    # run ldsc.py
    /lustre/groups/itg/shared/software/ldsc/ldsc.py --rg permuted.${t}.${n}.ldsc.sumstats.gz,T2D.ldsc.sumstats.gz --ref-ld-chr ${DATA_DIR}/LDscore/eur_w_ld_chr/ --w-ld-chr ${DATA_DIR}/LDscore/eur_w_ld_chr/ --out ${DATA_DIR}/LDscore/out/${t}.perm${n}.correlation --pop-prev 0.0873,0.1 --samp-prev ${s},0.089956
    echo "$(grep "permuted.${t}.${n}.ldsc.sumstats.gz  T2D.ldsc.sumstats.gz" ${DATA_DIR}/LDscore/out/${t}.perm${n}.correlation.log)" >> ${DATA_DIR}/LDscore/out/${t}.correlation.perm.out
    rm ${DATA_DIR}/LDscore/out/${t}.perm${n}.correlation.log
    rm ${DATA_DIR}/LDscore/permuted.${t}.${n}.ldsc.sumstats.gz
}

for trait in AllOA,0.27 KneeHipOA,0.224 TJR,0.125; do #HipOA,0.115 TKR,0.078 THR,0.078; do
    IFS="," read t s <<< "${trait}"
    echo "p1 p2 rg se z p h2_liab h2_liab_se h2_int h2_int_se gcov_int gcov_int_se" > ${DATA_DIR}/LDscore/out/${t}.correlation.perm.out
    for n1 in {1..40}; do
        for n2 in {1..250}; do 
            n=$(((${n1}-1)*250+${n2}))
            my.func "${t}" "${s}" "${n}" &
        done
        wait
    done
done




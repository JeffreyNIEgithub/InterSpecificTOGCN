[ $# -eq 0 ] && echo -e "\n# One '.RBH' and two '.SBH' files are needed.\n\nbash TOgenepair.sh \\ \n\tsp1_sp2.RHB \\ \n\tsp1_sp2.SBH \\ \n\tsp2_sp1.SBH\n" && exit

awk 'BEGIN{OFS="\t"}ARGIND==1{print $1,$2}ARGIND==2{print $1,$2}ARGIND==3{print $2,$1}' $1 $2 $3 | sort -u -k1 > sp1_sp2.geneid.txt

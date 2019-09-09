#!/bin/bash
touch atac.csv
echo 'library_id,fragments,cells' >> atac.csv
echo 'TILB1,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/TIL-B1/outs/fragments.tsv.gz,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/TIL-B1/outs/singlecell.csv' >> atac.csv
echo 'TILP1,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/TIL-P1/outs/fragments.tsv.gz,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/TIL-P1/outs/singlecell.csv' >> atac.csv
echo 'PTLB2,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/PTLB2/outs/fragments.tsv.gz,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/PTLB2/outs/singlecell.csv' >> atac.csv
echo 'PTLP2,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/PTLP2/outs/fragments.tsv.gz,/home/share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/out/PTLP2/outs/singlecell.csv' >> atac.csv

cellranger-atac aggr --id=atac \
                --csv=atac.csv \
                --normalize=depth \
                --reference=/home/swap/refdata-cellranger-hg19-1.2.0

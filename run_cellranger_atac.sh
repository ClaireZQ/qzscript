#!/bin/bash

for i in `echo PTLP2 TILB1 TILP1`;do
	./cellranger atac.sh $i /share/data3/SH_hjjyd/JYSH1903LHT29_Gaol_4scATAC/ $i 10 
done


  
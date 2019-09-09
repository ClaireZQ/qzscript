#!/bin/bash


for i in `echo B2 P1 P2`;do
	./Cellranger.sh $i /share/data3/SH_hjjyd/hjjydu_premature/JYSH1903LHT29_Gaol_4scRNA/Data/ $i 20 28 91
done

touch all.csv
echo 'library_id,molecule_h5'
echo 'B1,home/zhangquan/project/singlecell/uterus/shhjjy/B1/outs/molecule_info.h5'>>all.csv
echo 'B2,home/zhangquan/project/singlecell/uterus/shhjjy/B2/outs/molecule_info.h5'>>all.csv

cellranger aggr --id=aggop \
                --csv=all.csv \
                --normalize=mapped
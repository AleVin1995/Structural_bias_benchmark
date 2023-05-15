#!/usr/bin/env bash

source ~/.bashrc
conda activate CN_bench_r

ROOT=$1
libs=("Avana" "KY")

if [[ -z "$ROOT" ]]
then
	echo "Empty root path"
	exit 1
fi

# run CCR
for l in ${libs[@]}
do
	Rscript $ROOT/src/exec/run_CCR.r $ROOT/data/raw/"$l"_sgrna_raw_LFC.csv $ROOT/data/"$l"GuideMap.csv data/corrected/ "$l"
done

conda deactivate

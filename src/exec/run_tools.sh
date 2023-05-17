#!/usr/bin/env bash

source ~/.bashrc

ROOT=$1
libs=("Avana" "KY")

if [[ -z "$ROOT" ]]
then
	echo "Empty root path"
	exit 1
fi


# run CCR
run_CCR(){
	conda activate CN_bench_r

	for l in ${libs[@]}
	do
		Rscript $ROOT/src/exec/run_CCR.r $ROOT/data/raw/"$l"_sgrna_raw_LFC.csv $ROOT/data/"$l"GuideMap.csv data/corrected/ "$l"
	done

	conda deactivate
}


# run Chronos
run_Chronos(){
	conda activate CN_bench

	for l in ${libs[@]}
	do
		python3 $ROOT/src/exec/run_Chronos.py --lfc $ROOT/data/raw/"$l"_gene_raw_LFC.csv --cn $ROOT/data/OmicsCNGene.csv -o $ROOT/data/corrected/"$l"_gene_Chronos.csv
	done

	conda deactivate
}


# execute the functions
#run_CCR
run_Chronos
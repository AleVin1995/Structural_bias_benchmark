#!/usr/bin/env bash

#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/group/iorio/Alessandro/CN_benchmark
#SBATCH --output=/group/iorio/Alessandro/CN_benchmark/%x.out
#SBATCH --error=/group/iorio/Alessandro/CN_benchmark/%x.err

source ~/.bashrc

ROOT=$1
ALGO=$2
LIB=$3

if [[ -z "$ROOT" || -z "$ALGO" || -z "$LIB" ]]
then
	echo "One or more arguments are missing"
	echo "Usage: bash src/exec/run_tools.sh <root_path> <algorithm> <library>"
	exit 1
fi


# run CCR
run_CCR(){
	conda activate CN_bench_r

	Rscript $ROOT/src/exec/run_CCR.r \
		$ROOT/data/raw/"$LIB"_sgrna_raw_LFC.csv \
		$ROOT/data/"$LIB"GuideMap.csv \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run Chronos
run_Chronos(){
	conda activate CN_bench

	python3 $ROOT/src/exec/run_Chronos.py \
		--lfc $ROOT/data/raw/"$LIB"_gene_raw_LFC.csv \
		--cn $ROOT/data/OmicsCNGene.csv \
		-o $ROOT/data/corrected/"$LIB"_gene_Chronos.csv

	conda deactivate
}


# run Crispy
run_Crispy(){
	conda activate CN_bench

	python3 $ROOT/src/exec/run_Crispy.py \
		--lfc $ROOT/data/raw/"$LIB"_sgrna_raw_LFC.csv \
		--cn $ROOT/data/OmicsCNSegmentsProfile.csv \
		--map $ROOT/data/OmicsProfiles.csv \
		--lib $ROOT/data/"$LIB"GuideMap.csv \
		-o $ROOT/data/corrected/"$LIB"_gene_Crispy.csv

	conda deactivate
}


# run LDO
run_LDO(){
	conda activate CN_bench_r

	Rscript $ROOT/src/exec/run_LDO_and_GAM.r \
		$ROOT/data/raw/"$LIB"_lfc_exp_cn.rds \
		"LDO" \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run GAM
run_GAM(){
	conda activate CN_bench_r

	Rscript $ROOT/src/exec/run_LDO_and_GAM.r \
		$ROOT/data/"$LIB"_lfc_exp_cn.rds \
		"GAM" \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# execute the functions
if [[ "$ALGO" == "CCR" ]]
then
	run_CCR
elif [[ "$ALGO" == "Chronos" ]]
then
	run_Chronos
elif [[ "$ALGO" == "Crispy" ]]
then
	run_Crispy
elif [[ "$ALGO" == "LDO" ]]
then
	run_LDO
elif [[ "$ALGO" == "GAM" ]]
then
	run_GAM
else
	echo "Algorithm not found"
	echo "Usage: bash src/exec/run_tools.sh <root_path> <algorithm> <library>"
	exit 1
fi

# clean up
# rm *.out
# rm *.err
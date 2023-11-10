#!/usr/bin/env bash

#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/group/iorio/Alessandro_Vinceti/Structural_bias_benchmark
#SBATCH --output=/group/iorio/Alessandro_Vinceti/Structural_bias_benchmark/output/%j.out
#SBATCH --error=/group/iorio/Alessandro_Vinceti/Structural_bias_benchmark/error/%j.err

source ~/.bashrc

ROOT=$1
ALGO=$2
LIB=$3
JOB=${SLURM_ARRAY_JOB_ID}
ID=${SLURM_ARRAY_TASK_ID}

if [[ -z "$ROOT" || -z "$ALGO" || -z "$LIB" ]]
then
	echo "One or more arguments are missing"
	echo "Usage: bash src/method/run_tools.sh <root_path> <algorithm> <library>"
	exit 1
fi


# run CCR
run_CCR(){
	conda activate CN_bench_r

	Rscript $ROOT/src/method/run_CCR.r \
		$ROOT/data/raw/"$LIB"_sgrna_raw_LFC.csv \
		$ROOT/data/"$LIB"GuideMap.csv \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run Chronos
run_Chronos(){
	conda activate CN_bench

	python3 $ROOT/src/method/run_Chronos.py \
		--lfc $ROOT/data/raw/"$LIB"_screen_gene_effect.csv \
		--cn $ROOT/data/OmicsCNGene.csv \
		-o $ROOT/data/corrected/"$LIB"_gene_Chronos.csv

	conda deactivate
}


# run Crispy
run_Crispy(){
	conda activate CN_bench

	python3 $ROOT/src/method/run_Crispy.py \
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

	Rscript $ROOT/src/method/run_LDO_and_GAM.r \
		$ROOT/data/"$LIB"_lfc_exp_cn.rds \
		"LDO" \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run GAM
run_GAM(){
	conda activate CN_bench_r

	Rscript $ROOT/src/method/run_LDO_and_GAM.r \
		$ROOT/data/"$LIB"_lfc_exp_cn.rds \
		"GAM" \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run geometric
run_geometric(){
	conda activate CN_bench_r

	Rscript $ROOT/src/method/run_geometric.r \
		$ROOT/data/"$LIB"_lfc_exp_cn.rds \
		$ROOT/data/cytoband_mapping.csv \
		data/corrected/ \
		"$LIB"

	conda deactivate
}


# run MAGeCK
run_MAGeCK(){
	conda activate CN_bench

	python3 $ROOT/src/method/run_MAGeCK.py \
		mle \
		-k $ROOT/data/raw/"$LIB"_sgrna_raw_readcounts.csv \
		-s $ROOT/data/ScreenSequenceMap.csv \
		-n $ROOT/data/corrected/MAGeCK/"$LIB"_gene_MAGeCK_"$ID" \
		--seed $ID \
		--n-cells 50 \
		--cnv-norm $ROOT/data/OmicsCNGene.csv \
		--permutation-round 0 \
		--no-permutation-by-group

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
elif [[ "$ALGO" == "geometric" ]]
then
	run_geometric
elif [[ "$ALGO" == "CERES" ]]
then
	run_CERES
elif [[ "$ALGO" == "MAGeCK" ]]
then
	run_MAGeCK
else
	echo "Algorithm not found"
	echo "Usage: bash src/method/run_tools.sh <root_path> <algorithm> <library>"
	exit 1
fi

# clean up output and error files/directories
n_jobs=$(squeue -u Alessandro_Vinceti.vinceti | grep -w R | wc -l)

if [[ $n_jobs -eq 1 ]]
then
	rm -r output
	rm -r error
fi
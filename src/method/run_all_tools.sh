#!/usr/bin/env bash

#SBATCH --job-name=master
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/group/iorio/Alessandro/CN_benchmark
#SBATCH --output=/group/iorio/Alessandro/CN_benchmark/master.out
#SBATCH --error=/group/iorio/Alessandro/CN_benchmark/master.err

mkdir -p "output"
mkdir -p "error"
mkdir -p "data/corrected/MAGeCK"

for LIB in "Avana" "KY"
do
    # sbatch --job-name=CCR_"$LIB" --mem-per-cpu=10GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ CCR $LIB
    # sbatch --job-name=Chronos_"$LIB" --mem-per-cpu=100GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ Chronos $LIB
    # sbatch --job-name=Crispy_"$LIB" --mem-per-cpu=50GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ Crispy $LIB
    # sbatch --job-name=LDO_"$LIB" --mem-per-cpu=50GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ LDO $LIB
    # sbatch --job-name=GAM_"$LIB" --mem-per-cpu=10GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ GAM $LIB
    # sbatch --job-name=geom_"$LIB" --mem-per-cpu=50GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ geometric $LIB
    # sbatch --job-name=CERES_"$LIB" --mem-per-cpu=100GB src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ CERES $LIB
    sbatch --job-name=MAGeCK_"$LIB" --mem-per-cpu=100GB --array=1-100 src/method/run_tools.sh /group/iorio/Alessandro/CN_benchmark/ MAGeCK $LIB
done


# clean up
rm master*

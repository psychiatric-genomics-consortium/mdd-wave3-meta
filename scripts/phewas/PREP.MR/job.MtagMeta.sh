#!/bin/sh
#$ -cwd
#$ -N mtag
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -l h_rt=48:00:00
. /etc/profile.d/modules.sh

#module load module load igmm/apps/R/3.6.1
#module load anaconda/5.3.1
#conda config --add envs_dirs /exports/igmm/eddie/GenScotDepression/shen/anaconda3/envs/
#conda config --add pkgs_dirs /exports/igmm/eddie/GenScotDepression/shen/anaconda3/pkgs

conda activate mtag

touch data/MR/running_mtag

while read -r a b c || [ -n "$a" ]; do
  if (( $(bc <<< "$(grep $c data/MR/running_mtag | wc -l)==0")  )) && [ ! -f "${c}_mtag_meta.txt" ]; then
	echo $c >> data/MR/running_mtag

   	fname="$(basename $c)"
	echo $fname
	bash scripts/phewas/PREP.MR/func.mtag_meta.sh ${a} ${b} ${c} 
  fi
done < data/MR/mtag_list


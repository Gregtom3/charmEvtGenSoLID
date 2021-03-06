#!/bin/bash

workdir="/work/halla/solid/gam62/Research/charmEvtGenSoLID/psi2S-eventgen"

beamE=(11 12 15 17 24)
nuc=(0 1)
nuc_str=("p" "D")
type=("photo" "electro")
n=(5000000 100000000)

runJobs="${workdir}/slurm/runJobs.sh"
touch $runJobs
chmod +x $runJobs
echo " " > $runJobs

for energy in "${beamE[@]}"
do
    for i in 0 1
    do
	for j in 0 1
	do
	    file="${workdir}/slurm/shells/${nuc_str[${i}]}_${type[${j}]}_$energy.sh"
	    touch $file
	    echo "#!/bin/bash" > $file
	    echo "#SBATCH --account=halla" >> $file
	    echo "#SBATCH --partition=production" >> $file
	    echo "#SBATCH --mem-per-cpu=5000" >> $file
	    echo "#SBATCH --job-name=psi2S_${nuc_str[${i}]}_${type[${j}]}_${energy}" >> $file
	    echo "#SBATCH --time=01:00:00" >> $file
	    echo "#SBATCH --chdir=${workdir}" >> $file
	    echo "#SBATCH --output=${workdir}/slurm/output/%x-%j-%N.out" >> $file
	    echo "#SBATCH --error=${workdir}/slurm/output/%x-%j-%N.err" >> $file
	    echo "root -l ${type[${j}]}_solid_study_psi2S.C\\(${energy},${nuc[${i}]},${n[${j}]}\\)" >> $file

	    echo "sbatch shells/${nuc_str[${i}]}_${type[${j}]}_$energy.sh" >> $runJobs
	done 
    done
done
    

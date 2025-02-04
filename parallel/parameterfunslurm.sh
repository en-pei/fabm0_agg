#!/bin/bash
#SBATCH --job-name=testser
#SBATCH --partition=pNode,pCluster
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks=48 
#SBATCH --time=08:00:00
#SBATCH --mail-user=enpei.li@hzg.de
#SBATCH --mail-type=END,FAIL
#SBATCH --output=job.o%j
#SBATCH --error=job.e%j

module load applications/utils/cmake3-20.0
module load compilers/intel/2020.1.217
module load intelmpi
module load hdf5/1.10.5
module load netcdf/4.7.0
module load pnetcdf
module load applications/python/3.8


cd /gpfs/work/lie/bale2002

paramvar="breakup_factor" 
paramrange=$(seq 900 100 900) #60000 40000 180000 220000 40000 220000  14000 6000 20000 122000 50000 222000# #1500000 #1500  #1400 $(seq 13000 1 13000) #$(seq 13000 0.1 13000.1) #(seq 50 100 300) $(seq 11500 100 12000) 

paramvar1="coagulation_rate"
paramrange1=$(seq 51 1 51) #260 1 260 110 50 160 (seq 14 1 14)  #1300 #14, 5  #40 #|tr "," ".")   #needs to be changed manually in parameterfun.sh

paramvar2="ks" #"dens_lpm"
paramrange2=$(seq 1.6 0.01 1.6) #2000 100 3000

parampath="/gpfs/work/lie/bale2002/par/"  #paramstudy
unmodpath="/gpfs/work/lie/bale2002/" 

rm -rf $parampath
mkdir $parampath
cd $parampath


#paramval=$1
#parampath=$2
#unmodpath=$3
#paramvar=$4
#paramvar1=$5
#paramval1=$6
#paramrange1=$6
#paramvar2="dens_lpm" #$7
#paramrange2=$8


#currrun=$(echo $(date) $paramvar $paramval |sha512sum|cut -d ' ' -f1) #
for k in $paramrange; do
	for i in $paramrange1; do #260 1 260  #coagulation
		for j in $paramrange2; do  #dens_lpm	
			#currrun=$paramvar$k$paramvar1$i$paramvar2$j
			#pathn=$k"_"$i"_"$j"_"$currrun

                        #currrun=$k"_"$i"_"$j
			#pathn=$paramvar$k"_"$paramvar1$i"_"$paramvar2$j"_"$currrun

		        currrun=$paramvar"-"$k"-"$paramvar1"-"$i"-"$paramvar2"-"$j
		        pathn=$paramvar$k"-"$paramvar1$i"-"$paramvar2$j"-"$currrun
			echo $pathn
			echo $paramvar $k $paramvar1 $i $paramvar2 $j $parampath $umodpath $currrun 
				mkdir $parampath$currrun
				cd  $parampath$currrun
				cp  "$HOME/build/gotm" .   #changed path for gotm.exe
				echo "copying gotm"  >> $(echo $parampath$currrun".log")
			#cp $unmodpath"airsea.nml" .
			#cp $unmodpath"fabm_pelagic.nml" .
			#cp $unmodpath"gotm_fabm.nml" .
			#cp $unmodpath"gotmmean.nml" .
			#cp $unmodpath"gotmrun.nml" .
			#cp $unmodpath"gotmturb.nml" .
			#cp $unmodpath"obs.nml" .
			#cp $unmodpath"gotm_meteo.dat" .
				cp $unmodpath"fabm.yaml" .
				cp $unmodpath"outputall.yaml" output.yaml
		        	cp $unmodpath"gotm.yaml" .
				cp $unmodpath"gotm_meteo.dat" .

				echo "copying namelists"  >> $(echo $parampath$currrun".log")

				sed -i '/!/d' fabm.yaml 
		  		echo $k
		  		echo "modifing agg for variable "$paramvar" with value "$k" and "$paramvar1" with value "$i" " and "$paramvar2" with value "$j"  >> $(echo $parampath$currrun".log")
		  
				sed -i '/'$paramvar'/ c\''      '$paramvar': '$k fabm.yaml #'      '
				sed -i '/'$paramvar1'/ c\''      '$paramvar1': '$i fabm.yaml
				sed -i '/'$paramvar2'/ c\''      '$paramvar2': '$j fabm.yaml #added
		  		echo "starting fabm with "$paramvar" = "$paramval  >> $(echo $parampath$currrun".log")
		  		srun -n 1 ./gotm &>> $(echo $parampath$currrun"/"$paramvar$k$paramvar1$i$paramvar2$j"run.log") 
		  		echo "stopping fabm"  >> $(echo $parampath$currrun".log") 
		  #mv output.nc $paramvar$paramval"out.nc"
		  		echo "renameing output"  >> $(echo $parampath$currrun".log")
		  #ncpdq -O -d time,,,5 -a z,time bale2002.nc test.nc
		  	#ncap2 -A -s "totlpm$paramval=spm_spm+agg_agglpm" -v bale2002.nc ../test.nc
		  	#ncap2 -A -s "esd$paramval=agg_esd" -v bale2002.nc ../test.nc
				cd  ../
		done
	done
done
wait
cd $unmodpath
python3 Plot/postparameter.py 16 9
#cd Plot/evaluation
#cp ../../relativeerror.txt .
#python3 evaluation.py


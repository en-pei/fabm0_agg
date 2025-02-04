#!/bin/bash
####read parameter file from output of slurmeval to choose to plot#############
parampath="/gpfs/work/lie/run/par/test/" ###!!!changed #"/gpfs/work/lie/run/par/6oszmay/"  
balepath="/gpfs/work/lie/bale2002/" #"$HOME/setups/bale2002/" 
unmodpath="/gpfs/work/lie/oszmay_newphysics_minErrRange/" #/gpfs/work/lie/tmzjul/0coagulation_rate50.0/ 
postpath="/gpfs/work/lie/bale2002/Plot/" #"$HOME/setups/bale2002/Plot/"
#rm -rf $parampath
#mkdir $parampath

cd $parampath

while IFS= read -r f; do # IFS= read -r l <&3;
	echo "$f"
done < /gpfs/work/lie/results/paramtochangeoszmay.txt #< l.txt
#f=/gpfs/work/lie/results/paramtochange.txt
f=/gpfs/work/lie/results/paramtochangeoszmay.txt #therun.txt
#IN="40661 breakup_factor 7000.000000000001 coagulation_rate 30.0 ks 1.0 fractal_dimension1 2.1 pws 1.5 kws 350.0 0.4857305760713674 0.24531147541306306 0.7310420514844305 0.3398569469506252 8 1.3911851045338053"

#readarray -t a < /gpfs/work/lie/results/paramtochange.txt

IFS=$'\n' read -d '' -r -a lines < /gpfs/work/lie/results/paramtochangeoszmay.txt
printf '%s\n' "${lines[@]}"
echo "mark here"
echo "${lines[1]}"




n=1
while IFS= read -r "variable$n"; do
  n=$((n + 1))
done < /gpfs/work/lie/results/paramtochange6.txt
for ((i=0; i<4; i++)); do #0, 3
#for ((i=0; i<n-1; i++)); do 
   read -r ADDR0 ADDR1 ADDR2 ADDR3 ADDR4 ADDR5 ADDR6 ADDR7 ADDR8 ADDR9 ADDR10 ADDR11 ADDR12 <<<$(ifs=" ";echo ${lines[$i]}) #$IN
   echo "----> "$ADDR0 "line number:"$i
   name=("$ADDR1" "$ADDR3" "$ADDR5" "$ADDR7" "$ADDR9" "$ADDR11")
   value=("$ADDR2" "$ADDR4" "$ADDR6" "$ADDR8" "$ADDR10" "$ADDR12")
   dv=$i$ADDR3$ADDR4  #$ADDR0  #$i #
   mkdir $parampath$dv
   cd $parampath$dv
   echo "task nr." $dv; 
   cp $unmodpath"fabm.yaml" .
   #cp $unmodpath"output.yaml" .
   cp $balepath"output_more.yaml" output.yaml #outputall.yaml for all variables
   cp $unmodpath"gotm.yaml" .
   cp $unmodpath"gotm_meteo.dat" .
   for j in $(seq 0 1 5); do  #6 parameters from 0 
      sed -i '/'${name[$j]}'/ c\''      '${name[$j]}': '${value[$j]} fabm.yaml #'      '
      echo ${name[$j]}': '${value[$j]}
   done
#   ../../gotm
  ../../../gotm #run in current $d folder 
#  python3 $postpath"postparameter.py" 16 9 #$i  #try postprocessing 
  cd .. #back to test folder
#  echo $i>> $(echo "current.log") #give the run number
#  rm -rf $i  #testing removing folders
done
python3 /gpfs/work/lie/bale2002/Plot/plotoszmay.py 6 6 6
exit


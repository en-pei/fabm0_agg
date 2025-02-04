parampath="/gpfs/work/lie/run/test2/"  #"$HOME/setups/bale2002/test/"  
unmodpath="/gpfs/work/lie/oszmay_newphysics_minErrRange/" #"$HOME/setups/bale2002/" 
postpath="/gpfs/work/lie/bale2002/Plot/" #"$HOME/setups/bale2002/Plot/"
#rm -rf $parampath
#mkdir $parampath
cd $parampath

####recoded some part to make the last parameter constant but still rest working####
####simplified for parameter range definition#########

###needs input from slurm########3
ntask=$1 #how many cpu used  5
nid=$2  #cpu's number  2
echo "total cpu used are" $ntask
echo "current cpu is nr." $nid

p0="breakup_factor" 
#range0=($(seq 2000 1000 5000 )) #4000   10000     500 500 4000  2000 1000 8000
lk0=8 #5 #$( bc -l <<<"${#range0[@]}" ) #do not -1
delt0=300 #$( bc -l <<<"(${range0[-1]}-${range0[0]}) / ($lk0-1)" )  ##-1 need to be added now
min0=200 #500 #3000 #${range0[0]}

p1="coagulation_rate"
#range1=($(seq 50 10 80 )) #1 2 9   100    10 30 100  10 50 210  30 10 80
lk1=8 #10 #5 #$( bc -l <<<"${#range1[@]}" )
delt1=10 #$( bc -l <<<"(${range1[-1]}-${range1[0]}) / ($lk1-1)" )
min1=70 #50 #80  #120 #150#${range1[0]}

p2="ks" #"dens_lpm" #
#range2=($(seq 1.8 0.2 2.6 )) #0.01 0.01 0.5    0.01 0.1 2.01  0.01 0.5 2.01   1 0.1 2
lk2=8 #5 #$( bc -l <<<"${#range2[@]}" )
delt2=0.0002 #0.0002 #$( bc -l <<<"(${range2[-1]}-${range2[0]}) / ($lk2-1)" ) 
min2=0.0004 #0 #0.0003  #0.001 #${range2[0]}

p3="fractal_dimension1"
#range3=($(seq 1.5 0.2 2.1 ))  #1.5 0.5 2.5  1.8 0.1 2.2
lk3=5 #$( bc -l <<<"${#range3[@]}" )
delt3=0.2 #$( bc -l <<<"(${range3[-1]}-${range3[0]}) / ($lk3-1)" )
min3=1.8 #2.3 #${range3[0]}

p4="pws"
#range4=($(seq 1.45 0.02 1.55 )) #1.54   1.1 0.1 1.8   1.4 0.05 1.6
lk4=3 #4 # $( bc -l <<<"${#range4[@]}" )
delt4=0.1 #$( bc -l <<<"(${range4[-1]}-${range4[0]}) / ($lk4-1)" )
min4=1.6 #1.5 ${range4[0]}

p5="kws" #"kd" #
#range5=($(seq 220 20 280 ))  #347  100 200 1500   100 700 1500  200 50 400
lk5=5 #$( bc -l <<<"${#range5[@]}" )
delt5=500 #400 #$( bc -l <<<"(${range5[-1]}-${range5[0]}) / ($lk5-1)" )
min5=1700 #2000 #2300 #${range5[0]}

#p6="rho" #"Xsize" #"kd" #
#lk6=1 #5 #$( bc -l <<<"${#range6[@]}" )
#delt6=400 #0.1 #$( bc -l <<<"(${range6[-1]}-${range6[0]}) / ($lk6-1)" )
#min6=1600 #0.002 #0.5 #${range6[0]}


ntot=$( bc -l <<<"($lk0)*($lk1)*($lk2)*($lk3)*($lk4)*($lk5)" ) #*($lk6)" ) #total number of runs  
echo "total task number is" $ntot "."

nload=$(($ntot / $ntask))
echo "load for each cpu is" $nload "."
#########!!!!!!!!!!!!!!!!!!!
#exit  #!!!!!!for checking cpu number dividable before sbatch!!!!!!!!!!!!!!!
###!!!!!!!!####

if (( $nid > $(($ntask - 1)) )); then
    echo "cpu number larger than task number." 
    exit
fi

##when ntot is proportional to ntask###
#nid=$(( $nidgiven < $ntot ? $nidgiven : $ntot ))

range0=$(($nid * $nload))  #$(($(($nid * $nload))-$nload))  #current cpu's starting task number  
range1=$(($(($range0 + $nload)) -1)) #range0 + load -1 #current's ending task number
echo "task number is from" $range0 "to" $range1 "." #>> $(echo "$unmodpath/current.log")


num_variat=7 #4 #7 #6 #3  #varying 3 parameters
#args=()
variat_steps=("$lk0" "$lk1" "$lk2" "$lk3" "$lk4" "$lk5") # "$lk6") #step number -1 previously, no -1 now
variat_delt=("$delt0" "$delt1" "$delt2" "$delt3" "$delt4" "$delt5") # "$delt6") #interval
minv=("$min0" "$min1" "$min2" "$min3" "$min4" "$min5") # "$min6")
name=("$p0" "$p1" "$p2" "$p3" "$p4" "$p5") # "$p6")
dvl=1
echo ${variat_delt[@]} 


for i in $(seq $range0 1 $range1); do 
 dv=$i
 mkdir $parampath$dv
 cd $parampath$dv
 echo "task nr." $dv; 
 ln -s /gpfs/work/lie/bale2002/output.yaml
 ln -s /gpfs/work/lie/oszmay_newphysics_minErrRange/gotm.yaml
 cp $unmodpath"fabm.yaml" .
# cp $unmodpath"output.yaml" .
# cp $unmodpath"gotm.yaml" .
# cp $unmodpath"gotm_meteo.dat" .
 ln -s /gpfs/work/lie/bale2002/gotm_meteo.dat ##careful for the path!!!!!!!!

  dvl=1 #added
  for ((d=0; d<num_variat;d++)); do
    let vii="${variat_steps[$d]} * $dvl"
    #echo "vvi is:" $vvi "dvl is:" $dvl 
    vi=$(($dv % $vii))
    int=$(($vi / $dvl)) 
    vi=${int%.*}   #take the integer
    dv=$(($dv - $(($vi*$dvl))))  #from last task backwards
    dvl=$( bc -l <<<"${variat_steps[$d]}*$dvl" )   #added! 

#    dvl=$((${variat_steps[$d]} * $dvl)) #? why what is dvl multiplying? 
#    let f="${variat_delt[$d]} * ${variat_steps[$d]} + ${minv[$d]}"
#    let f="${variat_delt[$d]} * $vi + ${minv[$d]}"

    f=$( bc -l <<<"${variat_delt[$d]} * $vi + ${minv[$d]}" ) #changed to decimal friendly
#    echo $f ${name[$d]}
#    echo $d":"$vi":"$dv"\t"${name[$d]}":"$f 
    echo $i"-"${name[$d]}"-"$f>> $(echo "../current.log")  #$i "-" saving the run number out of the loop

    sed -i '/'${name[$d]}'/ c\''      '${name[$d]}': '$f fabm.yaml #'      '
  done

###here disabled running and postprocessing for just checking parameters##
  ../../gotm #run in current $d folder 
  python3 $postpath"postparameter2oszmay.py" $i  #new postparameterange 

  cd .. #back to test2 folder
  rm -rf $i  #testing removing folders


done
wait
#cd $unmodpath
#dv=12 #testing, dv as index from ntot 
exit

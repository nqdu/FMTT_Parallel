########################################################
#           Synthetic Traveltimes for A Given Model
# Before running this script, please make sure that:
# 	1. you have finished generating vtx model file, and modify
# 		the vtxfile variable below
#	2. you have generate reference time by Taup or  
#		generate it by using this script
########################################################

################################################
#	Variables defined
#################################################

# true and reference model defined
truemod=gridt.vtx
refmod=gridi.vtx

# whether to synthetic observed data(otimes.dat) (yes->true no->false)
syn=false

# no. of processes used
nproc=4

# command to synthetic data
cmd=../bin/fm3dt

#########################################################
# Variables end
###########################################################

# get command path
array=(${cmd//fm3dt/ })
path=''
for var in ${array[@]}
do
  path=${path}${var}
done

# if reference time required, synthetic it first
if [ $syn == true ];then
    echo "synthetic reference traveltime"
    cp $refmod gridc.vtx
    python parallel.py $cmd 1 $nproc 
    cp rtravel.out rtimes.dat
    cp rtimes.dat mkdata
    echo "synthetic true traveltime"
    cp $truemod gridc.vtx
    python parallel.py $cmd 2 $nproc
    cp rtravel.out mkdata
    cd mkdata
    ../${path}syntht
    cp otimes.dat ..
    cd ../../
else
    echo "synthetic traveltime"
    cp $refmod gridc.vtx
    python parallel.py $cmd 1 $nproc 
fi
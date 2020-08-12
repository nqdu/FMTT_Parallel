########################################################
#           Fmtt Inversion module
# Before running this script, please make sure that:
# 	1. you have finished generating initial vtx model file, and modify
# 		the vtxfile variable below
#	2. you have generate reference time by Taup or  
#		by using the synthetic module
########################################################

################################################
#	Variables defined
#################################################

# initial model defined
initmod=gridi.vtx

# no. of processes used
nproc=4

# command to synthetic theoretical traveltime
cmd=../bin/fm3dt

# maxiterations
maxiter=7

#########################################################
# Variables end
###########################################################

# copy initial model to gridc.vtx
cp $initmod gridc.vtx
mkdir -p models

# get command path
array=(${cmd//fm3dt/ })
path=''
for var in ${array[@]}
do
  path=${path}${var}
done

# remove residual file
if [  -f "residuals.dat" ];then
    rm residuals.dat
fi

# inversion begin
istep=1
for((i=1;i<=$maxiter;i++));
do
    echo "Iterations $i :"
    echo "----------------------------"
    echo "computing traveltimes by FMM ..."
    j=`echo $i+$istep |bc`
    python parallel.py $cmd $j $nproc

    echo 'solving linear systems by subspace inversion, rms = ...'
    echo $i+1|bc > subiter.in
    ${path}subinv
    #${path}fmttinv

    # save current results
    cp gridc.vtx models/grid_iter_$i.vtx
    
    # print rms
    ${path}residualst  >> residuals.dat
    tail -1 residuals.dat

    echo " "
done
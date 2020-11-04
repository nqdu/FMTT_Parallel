dampmin=0.01
dampmax=100.
ndamp=20

# no. of processes used, should be multiples of num_parts
nproc=16

# command to synthetic theoretical traveltime
cmd=../bin/fm3dt

# copy subinv.in
cp subinv.in subinv.in.temp

# synthetic data
python parallel.py $cmd 2 $nproc

# loop around damp
:>stats.dat
for((i=0;i<$ndamp;i++));do
    
done
# FMTT_Parallel
A parallel version of FMTT
Here I made several modifications:  
1. Parallelization of this module. See src/synthetic.sh and src/inversion.sh
2. If one ray pass through the sides of the 3-D model region, this ray will be discarded.
3. Add LSMR solver to solve least-square problems.

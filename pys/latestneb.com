# Inputfile for NEB calculation

# Default: bandpointnum=10 cyclelimit=32 theory=pm3 hightheory=none
#          theorychange="no",theorychange_cyclenum=20
#          algorythm=QuickMinVerlet 
#          option=(blank)
#          * option can be input as continuous phrases 
#            connected with  "___". 
#          e.g. option scrf=(pcm,solvent=ch2cl2)___int=grid=ultrafine
#          nproc=1,mem=1gb,qmsave=off,terminalsave=off,drag_abnormal=off
#          opt_turnfinenum=50,opt_timestep=turnfine,k_spring_max=0.05
#          k_spring=0.05,conver_coordinate=0.001,conver_force=0.001
#          logall=no
#          * logall can be "yes","no" and "ten". It ten, the log will 
#            be printerd per 10 calculation cycles. 
#          k_spring=0.05 (spring constant)
#          * variable spring constants can be used by "k_springs variable".
#            delta k can be specified using k_delta 0.02 etc.            


bandpointnum 10 
kyclelimit 128 
theory m062x/6-31G* 
nproc 2
hightheory m062x/6-31G*
algorythm QuickMinVerlet
qmsave off
logall ten 
k_spring variable
k_spring_max 0.1
k_delta 0.05

input terminal1
0 1
 C 0.0 0.0 0.0
 N 1.3 0.1 0.1
 H 2.3 0.0 0.0
end

input terminal2
0 1
 C 0.0 0.0 0.0
 N -1.3 0.1 0.13
 H 1.0 0.0 0.0
end


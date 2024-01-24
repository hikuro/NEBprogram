# Inputfile for NEB calculation

# Default Parameters: 
# bandpointnum = 10
# cyclelimit = 128
# opt_turnfinenum = 50
# opt_timestep = "turnfine"
#
# theory = "pm3"
# hightheory = "none"
# opttheory = "HF/3-21G"
# theorychange = "no"
# theorychange_cyclenum = 20
# algorythm = "quickmin_verlet"
# option = ""
#
# nproc = 1
# mem = "1gb"
#
# qmsave = "off"
# terminalsave = "off"
# drag_abnormal = "off"
# k_spring_max = 0.05
# k_spring = 0.05
# conver_coordinate = 0.0001
# conver_force = 0.0001
# logall = "no"
# k_spring = 0.05
# k_spring_max = 0.05
# k_delta = 0.02
# cineb_startnum = 5
# opt_timestep_max = 5E-15
# opt_maxstep = 0.1
#
#    * option can be input as continuous phrases 
#      connected with  "___". 
#      e.g. option scrf=(pcm,solvent=ch2cl2)___int=grid=ultrafine
#
#    * logall can be "yes","no" and "ten". It ten, the log will
#      be printerd per 10 calculation cycles.
#
#    * variable spring constants can be used by "k_spring variable".
#      delta k can be specified using k_delta 0.02 etc.

bandpointnum 2 
bandgen_method near_edge
bandgen_splitting on
bandgen_splitting_number 7 
bandgen_splitting_period 25

cyclelimit 300 
theory uff
hightheory pm7
theorychange_cyclenum 250
theorychange on
nproc 2
qmsave off
k_spring_max  0.05
terminalsave on
logall ten
k_spring variable
dihedral_force on
dihedral_force_constant 0.01
dihedral_force_endstep 300
terminalfix on

input terminal1
0 1
C    -1.089191   0.611765  -0.505681  
C    -2.405935   1.007557  -0.254690  
C    -0.692919  -0.679681  -0.154001  
C    -3.304081   0.122306   0.335378  
C    -1.592649  -1.573499   0.424231  
C    -2.901125  -1.169017   0.674425  
C    -0.124884   1.572113  -1.130663  
C     0.586890   2.405967  -0.051851  
C     1.853265   1.735569   0.491015  
C     1.568904   0.430797   1.264514  
N     1.695242  -0.775362   0.442187  
C     3.061434  -1.179760   0.085092  
C     0.687733  -1.158473  -0.434891  
O     0.932588  -1.936307  -1.336096  
H    -2.729077   2.010942  -0.522844  
H    -4.327783   0.438805   0.528410  
H    -1.267798  -2.585479   0.669727  
H    -3.607615  -1.860866   1.127425  
H     0.618614   1.037380  -1.759003  
H    -0.654273   2.245130  -1.836433  
H    -0.113523   2.624784   0.777609  
H     0.857507   3.392996  -0.477687  
H     2.560008   1.536669  -0.337613  
H     2.377929   2.444737   1.160412  
H     2.272423   0.337700   2.123560  
H     0.547720   0.458724   1.710793  
H     3.754952  -0.973237   0.913279  
H     3.078249  -2.270095  -0.114638  
H     3.420611  -0.674069  -0.826807  
end

input terminal2
0 1
C     1.089175   0.611757  -0.505689  
C     2.405917   1.007567  -0.254711  
C     0.692921  -0.679689  -0.153986  
C     3.304077   0.122337   0.335365  
C     1.592668  -1.573487   0.424252  
C     2.901140  -1.168986   0.674433  
C     0.124859   1.572088  -1.130682  
C    -0.586883   2.405985  -0.051883  
C    -1.853262   1.735624   0.491019  
C    -1.568910   0.430832   1.264489  
N    -1.695263  -0.775281   0.442108  
C    -3.061445  -1.179774   0.085083  
C    -0.687721  -1.158518  -0.434856  
O    -0.932549  -1.936467  -1.335973  
H     2.729045   2.010953  -0.522884  
H     4.327777   0.438852   0.528386  
H     1.267832  -2.585470   0.669760  
H     3.607641  -1.860820   1.127438  
H     0.654237   2.245076  -1.836489  
H    -0.618658   1.037338  -1.758985  
H    -0.857488   3.393008  -0.477741  
H     0.113546   2.624811   0.777561  
H    -2.377882   2.444799   1.160442  
H    -2.560036   1.536748  -0.337590  
H    -0.547724   0.458741   1.710766  
H    -2.272426   0.337713   2.123534  
H    -3.420248  -0.674932  -0.827445  
H    -3.078433  -2.270299  -0.113555  
H    -3.755174  -0.972238   0.912835  
end

connectivity
1   2   1.5 3   1.5 7   1           
2   4   1.5 15  1                   
3   5   1.5 13  1                   
4   6   1.5 16  1                   
5   6   1.5 17  1                   
6   18  1                           
7   8   1   19  1   20  1           
8   9   1   21  1   22  1           
9   10  1   23  1   24  1           
10  11  1   25  1   26  1           
11  12  1   13  1                   
12  27  1   28  1   29  1           
13  14  2                           
14                                      
15                                      
16                                      
17                                      
18                                      
19                                      
20                                      
21                                      
22                                      
23                                      
24                                      
25                                      
26                                      
27                                      
28                                      
29                                      
end




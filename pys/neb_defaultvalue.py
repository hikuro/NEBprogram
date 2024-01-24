def defaultvalue(keyword):
    ### Initialization
    if keyword == "bandpointnum":
        return int(10) 

    if keyword == "bandgen_method":
        return "none" #otheroptions, near_edge, dihedral

    # option for bandgen_dihedral
    if keyword == "bandgen_dihedral_rotnum":
        return int(1)

    if keyword == "bandgen_dihedral_threshold":
        return float(60)

    if keyword == "bandgen_dihedral_pathway":
        return "smallturn"

    if keyword == "bandgen_splitting":
        return "none"# This method can be used only when bandgen_method is near_edge

    if keyword == "bandgen_splitting_number":
        return int(3) # The band splits to 2*(3)+4 structures array

    if keyword == "bandgen_splitting_period":
        return int(10) # The top of the band will split every 5th iterations up to 2N+4

    if keyword == "dihedral_force_constant":
        return float(0.001)


    ### NEB options
    # Force calculation
    if keyword == "program":
        return "gaussian"

    if keyword == "theory":
        return "pm3"

    if keyword == "hightheory":
        return "none"

    if keyword == "opttheory":
        return "HF/3-21G"

    if keyword == "option":
        return ""

    if keyword == "nproc":
        return 1

    if keyword == "mem":
        return "4gb"

    if keyword == "multiprocessing":
        return "off"

    if keyword == "multiprocessing_totalcore":
        return 1

    # Theorychange
    if keyword == "theorychange":
        return "no"

    if keyword == "theorychange_cyclenum":
        return int(20) 

    # Optimization during NEB process
    if keyword == "cyclelimit":
        return int(128)

    if keyword == "opt_algorythm":
        return "quickmin_verlet"

    if keyword == "drag_abnormal":
        return "off"

    if keyword == "drag_start":
        return int(10)

    if keyword == "drag_threshold":
        return int(1000)

    if keyword == "opt_turnfinenum":
        return int(50)

    if keyword == "opt_timestep":
        return "turnfine"

    if keyword == "conver_coordinate":
        return float(0.0001)

    if keyword == "conver_force":
        return float(0.0001)

    if keyword == "k_spring":
        return float(0.05) 

    if keyword == "k_spring_max":
        return float(0.05)

    if keyword == "k_spring_max_init":
        return float(0.01)

    if keyword == "k_spring_max_shift":
        return 1 # If this is int value, the k_spring max starts from k_spring_max_init 
                     # and increases linearly up to k_spring_max at cyclenum "k_spring_max_shift."

    if keyword == "k_delta":
        return float(0.02)

    if keyword == "force_on_H":
        return "on"

    if keyword == "force_on_H_start":
        return int(9999) 

    if keyword == "cineb_startnum":
        return int(5)

    if keyword == "opt_timestep_max":
        return float(5E-15)

    if keyword == "opt_maxstep":
        return float(0.1)

    if keyword == "force_max":
        return float(0.01)

    if keyword == "dihedral_force":
        return "off"

    if keyword == "dihedral_force_endstep":
        return int(50)

    if keyword == "dihedral_force_constant":
        return float(0.001)

    if keyword == "dihedral_force_on_H":
        return "on"

    if keyword == "terminalfix":
        return "off"

    if keyword == "terminaloptimize":
        return "off"

    # Strict optimization of geometries
    if keyword == "optimize_theory":
        return "none"

    if keyword == "optimize_option":
        return ""

    if keyword == "optimize_nproc":
        return "1"

    if keyword == "optimize_mem":
        return "1gb"

    # Log save option
    if keyword == "qmsave":
        return "off"

    if keyword == "terminalsave":
        return "off"

    if keyword == "logall":
        return "no"

    if keyword == "forceprint":
        return "no"

    if keyword == "velocitysave":
        return "no"


# neb_mainprocess.py
# ver 0.5 2020/09/14 
# ver 0.6 2020/10/26
# written by Hiroaki Kurouchi

# Standard library
import os
import sys
import numpy as np
import shutil
import subprocess
import linecache
import copy
import datetime
import concurrent.futures

# Original Modules
import gausslog
import neb_reader
import neb_function
import neb_dihedral_forcemake
import neb_defaultvalue

### input paths and filenames ###
# The path must end with slash "/"
g16root      = "/local/apl/lx/g16a03/g16/"                  # CHECK HERE
#programdir   = "/home/users/to1/NEBprog/"                 # CHECK HERE
forceoption  = "force scf=(xqc,maxconv=35,fulllinear,nosym)"
g16env       = g16root + "bsd/g16.profile"        # CHECK HERE

try:
    subprocess.call("source " + g16env,shell=True)
except:
    print("g16env set failed")

mopacroot = "/Users/itsuu/MOPAC"
mopacprogram = "MOPAC2016.exe"
mopacforceoption = "XYZ GRADIENTS 1SCF PRTXYZ"


atom_to_weight = {'H':1.00783,'He':4.0026,'Li':6.941,'Be':9.012,'B':10.811,
            'C':12.,'N':14.007,'O':15.9994,'F':18.9984,'Ne':20.1797,
            'Na':22.989,'Mg':24.305,'Al':26.98154,'Si':28.0855,'P':30.9738,
            'S':32.066,'Cl':35.4527,'Ar':39.948,'K':39.0983,'Ca':40.078,
            'Sc':44.96,'Ti':47.867,'V':50.94,'Cr':51.9961,'Mn':54.938,
            'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
            'Ga':69.723,'Ge':72.64,'As':74.9216,'Se':78.96,'Br':79.904,
            'Pd':106.42,'I':126.90447}

origdir = os.getcwd()
if origdir[-1] != "/":
    origdir += "/"

# 基本的なループ
# Logfiles/bandarrayを読み込んで最適化する
def nebloop():
    inputfilename = neb_reader.get_inputfilename()

    # bandarrayには最新の構造情報が格納されているので読み込む
    atoms,chargespin,bandarray,energyarray = neb_reader.bandarray_read()
    newbandarray = 0
    cyclenum = 0
    cyclelimit = neb_reader.input_parameterread(inputfilename,"cyclelimit")
    qmsave = neb_reader.input_parameterread(inputfilename,"qmsave")
    bandgen_splitting = neb_reader.input_parameterread(inputfilename,"bandgen_splitting")
    bandgen_splitting_period = neb_reader.input_parameterread(inputfilename,"bandgen_splitting_period")
    terminalsave = neb_reader.input_parameterread(inputfilename,"terminalsave")
    drag_abnormal = neb_reader.input_parameterread(inputfilename,"drag_abnormal")
    drag_start = neb_reader.input_parameterread(inputfilename,"drag_start")
    velocityarray = [np.zeros([len(bandarray[0]),3]) for i in range(len(bandarray))]
    terminalfix = neb_reader.input_parameterread(inputfilename,"terminalfix")
    theorychange = neb_reader.input_parameterread(inputfilename,"theorychange")
    theorychange_cyclenum = neb_reader.input_parameterread(inputfilename,"theorychange_cyclenum")

    #計算途中で終わっていた場合にはrestartflagを読み込んでサイクル数を確認して再開
    if os.path.exists(origdir+"Logfiles/restartflag") == True:
        with open(origdir+"Logfiles/restartflag","r") as rr:
            if "nebdone" in rr.read():
                return 0
            else:
                cyclenum = int(rr.read())
        if os.path.exists(origdir+"Logfiles/velocityarray") == True:
            velocityarray = neb_reader.velocityarray_read()

    while True:
    # それぞれの構造を１つ前の隣接する構造にフィットさせ、速度に回転行列をかける。
        bandarray,velocityarray = neb_function.procrustes_superimposition(bandarray,velocityarray)

    # bandgen_splitting がONになっている場合、構造を分裂させていく
        if bandgen_splitting in ["on","yes"] and cyclenum % bandgen_splitting_period == 0:
            bandarray,velocityarray = splitting(bandarray,velocityarray,cyclenum)

    # それぞれの構造の力の計算、返ってくる値はNx3次元とN次元
        force_vector_array,energy_array = pes_analyze(cyclenum,atoms,chargespin,bandarray)

    # それぞれの構造にかかる力に対して接線方向に補正をかける(Method of Henkelman and Jonsson)
    # 関数内部では１次元化して処理しているが帰ってくる値はNx3次元
        neb_force_array = neb_forcemake(cyclenum,atoms,bandarray,force_vector_array,energy_array)

    # 構造の変更と新しいbandarrayの作成
        newbandarray,newvelocityarray = optimization(atoms,bandarray,velocityarray,neb_force_array,cyclenum)

    # 力の定数の評価、bandarray間での構造の変化量の評価
        terminated,evaluation = evaluate(newbandarray,bandarray,neb_force_array)


    # 基準を満たすか否か、cycleがどの程度回ったかで分岐
    # 最終的に得られたbandarrayとそのエネルギーを次に使うログファイルに保存
        if (theorychange in ["on","yes"] and  terminated == 1 
                            and cyclenum < theorychange_cyclenum):
                cyclenum = theorychange_cyclenum

        elif cyclenum >= cyclelimit or terminated == 1:
            logwrite(atoms,chargespin,bandarray,velocityarray,energy_array,
                     force_vector_array,neb_force_array,cyclenum,1)
            with open(origdir+"Logfiles/restartflag","w") as w:
                w.write("nebdone")
            
            shutil.copy(origdir+"Logfiles/bandarray",origdir+"Logfiles/bandarray_terminated")
            break

    # Cyclenum がdrag_startを超えてもなお構造に異常がある場合は強制的に構造を変化させる
        if 0 in energy_array and cyclenum > drag_start and drag_abnormal in ["on","yes"]:
            newbandarray = neb_function.nancheck_band2(bandarray,energy_array)

    # 基準を満たさなかった場合はcycleに１足して情報を更新して戻る
        if qmsave not in ["on", "yes"]:
           eracecomlog(bandarray,cyclenum) 
        cyclenum += 1
        bandarray = newbandarray
        velocityarray = newvelocityarray
        logwrite(atoms,chargespin,bandarray,velocityarray,energy_array,
                 force_vector_array,neb_force_array,cyclenum,0)            

    # 途中経過の保存
        bandarraywrite(atoms,chargespin,bandarray,energy_array)
        velocityarraywrite(atoms,chargespin,velocityarray)
        with open(origdir+"Logfiles/restartflag","w") as w:
            w.write(str(cyclenum))

        if terminalsave in ["yes", "on"]:
            terminalcoordsave(atoms,chargespin,bandarray) 

 
### 以下、関数

# origdirのinputfile,outputfileを消す関数
def eracecomlog(bandarray,cyclenum):
    inputfilename = neb_reader.get_inputfilename()
    terminalfix = neb_reader.input_parameterread(inputfilename,"terminalfix")
    program = neb_reader.input_parameterread(inputfilename,"program")
    for strnum in range(len(bandarray)):
        if terminalfix in ["on","yes"] and cyclenum == 0 and strnum in [0,len(bandarray)-1]:
            if strnum == 0:
                os.rename(origdir+"g16_0_0.out",origdir+"g16_terminal_1.out")
            else:
                os.rename(origdir+"g16_0_{}.out".format(len(bandarray)-1),origdir+"g16_terminal_2.out")

        elif program in ["none","g16","gaussian"]:
            try:
                os.remove(origdir+"g16_{}_{}.com".format(cyclenum,strnum)) 
                os.remove(origdir+"g16_{}_{}.out".format(cyclenum,strnum))
            except:
                pass

        elif program == "mopac":
            os.remove(origdir+"mopac_{}_{}.mop".format(cyclenum,strnum))
            os.remove(origdir+"mopac_{}_{}.out".format(cyclenum,strnum))


# 力の定数、エネルギーをそれぞれNx3次元のリストにまとめて返す関数
def pes_analyze(cyclenum,atoms,chargespin,bandarray):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")
    terminalfix = neb_reader.input_parameterread(inputfilename,"terminalfix")
#未実装    terminaloptimize = neb_reader.input_parameterread(inputfilename,"terminaloptimize")
    multiprocessing = neb_reader.input_parameterread(inputfilename,"multiprocessing")
    multiprocessing_totalcore = neb_reader.input_parameterread(inputfilename,"multiprocessing_totalcore")
    nproc = neb_reader.input_parameterread(inputfilename,"nproc")

    force_vector_array = [0 for i in range(len(bandarray))]
    energy_array = [0 for i in range(len(bandarray))]
    scratchdir = neb_reader.get_scratchdirpath()

    for strnum in range(len(bandarray)):
        commake(cyclenum,strnum,atoms,chargespin,bandarray[strnum])

#    sys.exit(1)
    startnumber = 0
    endnumber = len(bandarray)

    # terminalfix,terminaloptimizeの場合はterminal構造の計算をスキップする
    if terminalfix in ["on","yes"]  and cyclenum >= 1:
        g16_terminal1_output = "g16_terminal_1.out"
        g16_terminal2_output = "g16_terminal_2.out"
        force_vector_array[0],energy_array[0] = gausslog.forcereader(origdir+g16_terminal1_output)
        force_vector_array[-1],energy_array[-1] = gausslog.forcereader(origdir+g16_terminal2_output)
        force_vector_array[0] = np.array(force_vector_array[0])
        force_vector_array[-1] = np.array(force_vector_array[-1])
        startnumber = 1
        endnumber = len(bandarray) - 1

    #並列化計算と逐次計算で分岐
    if multiprocessing in ["on","yes"]:
#        executor = concurrent.futures.ProcessPoolExecutor()
        number_of_threads = int(multiprocessing_totalcore / nproc) 
        strnumarrays = [[] for i in range(number_of_threads)]
        for i in range(startnumber,endnumber):
            strnumarrays[i%number_of_threads].append(i)

        with concurrent.futures.ProcessPoolExecutor() as executor:
            for strnumarray in strnumarrays:
                executor.submit(multiprocess_forcecalc_run,  cyclenum,strnumarray)
    
    else: #逐次計算
        for strnum in range(startnumber,endnumber):
            g16_inputfilename = "g16_{}_{}.com".format(cyclenum,strnum) 
            g16_outputfilename = "g16_{}_{}.out".format(cyclenum,strnum)   
            forcecalc_run(g16_inputfilename,g16_outputfilename)

    for strnum in range(startnumber,endnumber):
        g16_outputfilename = "g16_{}_{}.out".format(cyclenum,strnum)
        force_vector_array[strnum],energy_array[strnum] = (
                gausslog.forcereader(scratchdir+g16_outputfilename))
        force_vector_array[strnum] = np.array(force_vector_array[strnum])
    
    energy_array = np.array(energy_array)
    return force_vector_array,energy_array

def multiprocess_forcecalc_run(cyclenum,strnumarray):
    for strnum in strnumarray:
        g16_inputfilename = "g16_{}_{}.com".format(cyclenum,strnum) 
        g16_outputfilename = "g16_{}_{}.out".format(cyclenum,strnum)
        forcecalc_run(g16_inputfilename,g16_outputfilename)
 

def neb_forcemake(cyclenum,atoms,bandarray,force_vector_array,energy_array): 
    inputfilename = neb_reader.get_inputfilename()
    ks_max = neb_reader.input_parameterread(inputfilename,"k_spring_max")
    ks_max_init = neb_reader.input_parameterread(inputfilename,"k_spring_max_init")
    ks_max_shift = neb_reader.input_parameterread(inputfilename,"k_spring_max_shift")
    ks = neb_reader.input_parameterread(inputfilename,"k_spring")
    dk = neb_reader.input_parameterread(inputfilename,"k_delta")
    force_max = neb_reader.input_parameterread(inputfilename,"force_max")
    force_on_H = neb_reader.input_parameterread(inputfilename,"force_on_H")
    force_on_H_start = neb_reader.input_parameterread(inputfilename,"force_on_H_start")
    dihedral_force = neb_reader.input_parameterread(inputfilename,"dihedral_force")
    dihedral_force_endstep = neb_reader.input_parameterread(inputfilename,"dihedral_force_endstep")
    dihedral_force_constant = neb_reader.input_parameterread(inputfilename,"dihedral_force_constant")
    cineb = neb_reader.input_parameterread(inputfilename,"cineb")
    N_cineb = neb_reader.input_parameterread(inputfilename,"cineb_startnum")
    terminalfix = neb_reader.input_parameterread(inputfilename,"terminalfix")

    if ks == "variable":
        ks = "variable"
    else:
        ks = float(ks)

    if float(dk) >=  ks_max:
        dk = ks_max/2 
    else:
        dk = float(dk)

    if ks_max_shift != "no":
        if cyclenum < ks_max_shift:
            ks_max = ks_max_init + (ks_max-ks_max_init)*(cyclenum/ks_max_shift)
        else:
            pass

    neb_force_array = []
    #計算で得られた力の定数にbandの張力を追加する
    for i in range(len(force_vector_array)):
        #ここで力の定数を１次元ベクトル化
        force_i = np.ravel(force_vector_array[i])

        #末端は力の定数に手を加えずに用いる
        if i == 0 or i >= len(force_vector_array)-1:
            force = copy.copy(force_i)
        else:
            #Definition of Vdiff: ndarray(2x1)
            # [energy_array[i]-energy_array[i-1],energy_array[i+1]-energy_array[i]]
            V_diff = energy_array[i:i+2] - energy_array[i-1:i+1]

            bandarray_diff = [bandarray[i] - bandarray[i-1],bandarray[i+1] - bandarray[i]]

            #以降、１次元化
            if np.sum(V_diff > 0) == 2:
                tangential_vector = np.ravel(bandarray_diff[1])
            
            elif np.sum(V_diff > 0) == 0: 
                tangential_vector = np.ravel(bandarray_diff[0])

            else:
                if energy_array[i+1] > energy_array[i-1]:
                    tangential_vector = ( np.min(np.abs(V_diff)) * np.ravel(bandarray_diff[0])
                                        + np.max(np.abs(V_diff)) * np.ravel(bandarray_diff[1]) )

                else:
                    tangential_vector = ( np.max(np.abs(V_diff)) * np.ravel(bandarray_diff[0])
                                        + np.min(np.abs(V_diff)) * np.ravel(bandarray_diff[1]) )

            #1次元接線ベクトルの正規化
            tangential_vector /= np.sum(tangential_vector ** 2) ** 0.5

            #1次元化された接線方向の力の計算
            # Variable spring法の場合
            if ks == "variable":
                E_max = np.max(energy_array)
                E_ref = np.max(np.array([energy_array[0],energy_array[-1]]))
                E_i = np.max(energy_array[i-1:i+1])
                if E_i > E_ref:
                    ks = ks_max - dk * ((E_max - E_i) / (E_max - E_ref))
                else:
                    ks = ks_max - dk

            #NEB長さが変わらない場合ここでチェックが必要
            force_tangential = ( ks * (np.linalg.norm(bandarray_diff[1])
                                    -  np.linalg.norm(bandarray_diff[0])) 
                                    *  tangential_vector) 

            #垂直方向の力の計算(符号に注意)
            force_vertical = force_i - np.dot(force_i,tangential_vector) * tangential_vector

            #Cyclenumが5を超えたらClimbing Image NEBを適用する
            if cyclenum >= N_cineb and i == np.argmax(energy_array) and cineb not in  ["no","off"]:
                force = force_i - 2.0 * np.dot(force_i,tangential_vector) * tangential_vector
            else:
                force = force_vertical + force_tangential

        # 力の定数が大きすぎる場合は閾値以下にする
        force = force - (force - force_max) * (force > force_max) - (force + force_max) * (force < -force_max) 

        # 1次元化したものをNx3次元に戻す
        force = force.reshape([int(len(force)/3),3])
        neb_force_array.append(force)

    # Dihedral angleに沿った力をかける
    if dihedral_force in ["on","yes"] and cyclenum < dihedral_force_endstep:
        dihedral_force = neb_dihedral_forcemake.forcemake(bandarray,atoms)
        for index,force in enumerate(neb_force_array):
            force += dihedral_force[index]

    #水素原子にbandからの力をかけない設定
    if force_on_H in ["off","no","halt"] and cyclenum < force_on_H_start:
        for i,atom in enumerate(atoms):
            if atom == "H":
                for j in range(len(neb_force_array)):
                    neb_force_array[j][i] *= 0

    elif force_on_H in ["half"] and cyclenum < force_on_H_start:
        for i,atom in enumerate(atoms):
            if atom == "H":
                for j in range(len(neb_force_array)):
                    neb_force_array[j][i] /= 2

    if terminalfix in ["on","yes"]:
        neb_force_array[0] *= 0
        neb_force_array[-1] *= 0

    return neb_force_array


# bandarrayに力を与えて動かす関数。いくつかのアルゴリズムに対応。
# 現在実装されているアルゴリズム: Quick-Min Verlet
def optimization(atoms,bandarray,velocityarray,neb_force_array,cyclenum):
    inputfilename = neb_reader.get_inputfilename()
    atomicweights = [ atom_to_weight[atoms[i]] for i in range(len(atoms)) ]
    atomicweight_array = np.ravel(np.array([atomicweights,atomicweights,atomicweights]).T)
    opt_algorythm = neb_reader.input_parameterread(inputfilename,"opt_algorythm")
    opt_timestep = neb_reader.input_parameterread(inputfilename,"opt_timestep")
    opt_maxstep = neb_reader.input_parameterread(inputfilename,"opt_maxstep")
    opt_turnfinenum = neb_reader.input_parameterread(inputfilename,"opt_turnfinenum")
    force_on_H = neb_reader.input_parameterread(inputfilename,"force_on_H")
    force_on_H_start = neb_reader.input_parameterread(inputfilename,"force_on_H_start")
    terminalfix = neb_reader.input_parameterread(inputfilename,"terminalfix")

    if opt_timestep in  ["none","turnfine"] and cyclenum > opt_turnfinenum:
        opt_timestep = 1E-15 * opt_turnfinenum  / float(cyclenum)
    elif opt_timestep in ["none","turnfine"]  and cyclenum <= opt_turnfinenum:
        opt_timestep = 1E-15
    elif float(opt_timestep) > 5E-15:
        opt_timestep = neb_defaultvalue.defaultvalue("opt_timestep_max")
    else:
        opt_timestep = float(opt_timestep)

    if (opt_algorythm == "quickmin_verlet" and force_on_H in ["off","no","halt","half"]
         and cyclenum < force_on_H_start):
        opt_algorythm = "verlet"

    # Gaussianのforce outputの単位はHartrees/Bohr, 原子量はg/molなので
    # 力を原子量で割ったときに単位がA/S^2 になるよう係数をかける。
    # Coeffの単位はBohr•A•g/Hartree•mol•s^2
    Coeff = 1E20 * 627.509 * 4184 * 1000 / 0.529177
    for i in range(len(neb_force_array)):
        neb_force_array[i] *= Coeff

    new_bandarray = []
    new_velocityarray = []
    if opt_algorythm in ["quickmin_verlet","verlet"]:
        for i in range(len(bandarray)):
            if terminalfix in ["on","yes"] and i in [0,len(bandarray)-1]:
                new_band_temp = bandarray[i]
                new_velocity_temp = np.zeros([len(bandarray[i]),3])
            else:
                band_temp = np.ravel(bandarray[i])
                velocity_temp = np.ravel(velocityarray[i])
                neb_force_temp = np.ravel(neb_force_array[i])
                acceleration = neb_force_temp / atomicweight_array
                new_velocity_temp = (acceleration * opt_timestep
                                  + (np.dot(velocity_temp,neb_force_temp) > 0 ) 
                                  *  np.dot(velocity_temp,neb_force_temp) * neb_force_temp
                                  /  np.dot(neb_force_temp,neb_force_temp))

                if np.dot(velocity_temp,neb_force_temp) <= 0 and opt_algorythm == "quickmin_verlet":
                    print("cyclenum:{} strnum:{}".format(cyclenum,i),end="") 
                    print("  oldvelovity was set to 0 because directions of force and velocity are not identical")
                delta_s = new_velocity_temp * opt_timestep + 0.5 * acceleration * (opt_timestep)**2
                delta_s = ( (delta_s > opt_maxstep) * opt_maxstep 
                           - (delta_s < (-1 * opt_maxstep)) * opt_maxstep 
                           +(np.abs(delta_s) <= opt_maxstep) * delta_s )
                new_band_temp = band_temp + delta_s


                new_band_temp = new_band_temp.reshape([int(len(new_band_temp)/3),3])
                new_velocity_temp = new_velocity_temp.reshape([int(len(new_velocity_temp)/3),3])

            new_bandarray.append(new_band_temp)
            new_velocityarray.append(new_velocity_temp)

    else:
        print("No algorythm here!")
        sys.exit()

    return new_bandarray, new_velocityarray

def theorychoose(cyclenum):
    inputfilename = neb_reader.get_inputfilename()
    theorychange = neb_reader.input_parameterread(inputfilename,"theorychange")
    theorychange_cyclenum = neb_reader.input_parameterread(inputfilename,"theorychange_cyclenum")
    theory = neb_reader.input_parameterread(inputfilename,"theory")
    hightheory = neb_reader.input_parameterread(inputfilename,"hightheory")

    if theorychange in ["on","yes"]:
        if cyclenum >= theorychange_cyclenum:
#            with open(origdir+"Logfiles/neb_log","a") as lw:
#                lw.write(" Theorychange:on,  theory: {}\n".format(hightheory))
            return hightheory

        else:
#            with open(origdir+"Logfiles/neb_log","a") as lw:
#                lw.write(" Theorychange: on,  theory: {}\n".format(theory))
            return theory
    else:
        return theory


def commake(cyclenum,strnum,atoms,chargespin,structure):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    nproc = neb_reader.input_parameterread(inputfilename,"nproc")
    mem = neb_reader.input_parameterread(inputfilename,"mem") 
    connectivity = neb_reader.connectivity_read()

    # theoryに関してはswitchingするので関数内で定義する
    theory = theorychoose(cyclenum)

#    theory = neb_reader.input_parameterread(inputfilename,"theory")
    option = neb_reader.input_parameterread(inputfilename,"option").split("___")
    isconnected = 0
    for contents in option:
        if "connectivity" in contents:
            isconnected = 1

    if len(connectivity) > 1 and isconnected == 0 and program in ["none","g16","gaussian"]:
        if option not in [ "", [""]]:
            option.append("geom=connectivity")
        else:
            option = ["geom=connectivity"]

    if program in ["none","g16","gaussian"]:
        g16_inputfilename = "g16_{}_{}.com".format(cyclenum,strnum)
        with open(origdir+g16_inputfilename,"w") as g16w:
            g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
            g16w.write(" #p {} {}\n".format(theory,forceoption))
            if option not in [ "", [""]]:
                for i in range(len(option)):
                    g16w.write(option[i]+" ")
                g16w.write("\n")
            g16w.write("\n NEBcycle:{} Structurenumber:{}\n\n".format(cyclenum,strnum))
            g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
            for i in range(len(atoms)):
                g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                            atoms[i],structure[i][0],structure[i][1],structure[i][2]))

            g16w.write("\n")
            if len(connectivity) > 1:
                for con in connectivity:
                    g16w.write(con)
            g16w.write("\n")
            g16w.write("\n")

    if program == "mopac":
        mopac_inputfilename = "mopac_{}_{}.mop".format(cyclenum,strnum)
        with open(origdir+mopac_inputfilename,"w") as comwrite:
            comwrite.write("{} {} threads={} charge={}".format(theory,mopacforceoption,nproc,chargespin[0]))
            if chargespin[1] == 2:
                comwrite.write(" DOUBLET")
            elif chargespin[1] == 3:
                comwrite.write(" TRIPLET")
            comwrite.write(" \n")
            # Title section
            comwrite.write("NEBcycle:{} Structurenumber:{}\n\n".format(cyclenum,strnum))
            for i in range(len(atoms)):
                comwrite.write(" {:2}  {:9> .7f} {:9> .7f} {:9> .7f} ".format(\
                     atoms[i],structure[i][0],structure[i][1],structure[i][2]))
                comwrite.write("\n")
        comwrite.write("\n")
        comwrite.flush()


# neb_initializeの関数とほぼ同じだがbandarrayの名前が異なる
def bandarraywrite(atoms,chargespin,bandarray,energy_array):
    with open(origdir+"Logfiles/bandarray","w") as bw:
        for num in range(len(bandarray)):
            bw.write("{}\n{} {} str:{} energy:{:9> .7f} a.u.\n".format(
                    len(bandarray[0]),chargespin[0],chargespin[1],num,energy_array[num]))
            for i in range(len(bandarray[0])):
                bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                atoms[i],bandarray[num][i][0],bandarray[num][i][1],bandarray[num][i][2]))


def velocityarraywrite(atoms,chargespin,velocityarray):
    with open(origdir+"Logfiles/velocityarray","w") as bw:
        for num in range(len(velocityarray)):
            bw.write("{}\n{} {} structurenum:{}\n".format(len(velocityarray[0]),chargespin[0],chargespin[1],num))
            for i in range(len(velocityarray[0])):
                bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                atoms[i],velocityarray[num][i][0],velocityarray[num][i][1],velocityarray[num][i][2]))

def terminalcoordsave(atoms,chargespin,bandarray):
    with open(origdir+"Logfiles/terminal1log","a") as bw:
        bw.write("{}\n{} {} structurenum:{}\n".format(len(bandarray[0]),chargespin[0],chargespin[1],0))
        for i in range(len(bandarray[0])):
            bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
            atoms[i],bandarray[0][i][0],bandarray[0][i][1],bandarray[0][i][2]))

    with open(origdir+"Logfiles/terminal2log","a") as bw:
        bw.write("{}\n{} {} structurenum:{}\n".format(len(bandarray[0]),chargespin[0],chargespin[1],"final"))
        for i in range(len(bandarray[0])):
            bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
            atoms[i],bandarray[-1][i][0],bandarray[-1][i][1],bandarray[-1][i][2]))


def logwrite(atoms,chargespin,bandarray,velocityarray,energy_array,
             force_vector_array,neb_force_array,cyclenum,terminated=0):
    inputfilename = neb_reader.get_inputfilename()
    bandarray_distance = [0]
    for i in range(1,len(bandarray)):
        bandarray_diff = bandarray[i] - bandarray[i-1]
        bandarray_distance.append(np.linalg.norm(bandarray_diff)+bandarray_distance[-1])

    energy_array_diff_kcal = [0.0]
    for i in range(1,len(bandarray)):
        temp = energy_array[i] - energy_array[0]
        energy_array_diff_kcal.append(627.51*temp)

    logall = neb_reader.input_parameterread(inputfilename,"logall")
    forceprint = neb_reader.input_parameterread(inputfilename,"forceprint")
    cyclelimit = neb_reader.input_parameterread(inputfilename,"cyclelimit")
    conver_force = neb_reader.input_parameterread(inputfilename,"conver_force")
    velocitysave = neb_reader.input_parameterread(inputfilename,"velocitysave")

    with open(origdir+"Logfiles/neb_log","a") as lw:
        if terminated == 1:
            if cyclenum >= cyclelimit:
                lw.write("\n\n -------- NEBcycle Not Converged --------\n")
            else:
                lw.write("\n\n ---------- NEBcycle Converged ----------\n")
            time_now = datetime.datetime.now()
            time_now_str = time_now.strftime("%Y/%m/%d/ %H:%M")
            lw.write(" NEB converged time: {}\n\n".format(time_now_str))

        time_now = datetime.datetime.now()
        time_now_str = time_now.strftime("%Y/%m/%d/ %H:%M:%S")
        lw.write(    "\n -------------- NEBcycle {} {}-------------\n".format(cyclenum,time_now_str))
        for i in range(len(energy_array)):
            lw.write(" str {: >3}: s= {:.2f} pes= {: >9.6f} a.u. pesdiff= {: >5.2f} kcal/mol\n".format(
                    i,bandarray_distance[i]/bandarray_distance[-1],energy_array[i],energy_array_diff_kcal[i]))

        if (terminated == 1) or logall in ["yes","on"] or (logall in ["10","ten","perten"] and int(cyclenum)%10 == 0):
            lw.write("\n --- Coordinates ---\n")
            for num in range(len(bandarray)):
                lw.write(" {}\n {} {} structurenum:{}\n".format(len(bandarray[0]),chargespin[0],chargespin[1],num))
                for i in range(len(bandarray[0])):
                    lw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                    atoms[i],bandarray[num][i][0],bandarray[num][i][1],bandarray[num][i][2]))

            lw.write("\n --- Raw Force vector \n")
            for num in range(len(force_vector_array)):
                lw.write(" \n {} {} structurenum:{}\n".format(
                         chargespin[0],chargespin[1],num))
                if forceprint in ["on","yes"]:
                    for i in range(len(bandarray[0])):
                        lw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                        atoms[i],force_vector_array[num][i][0],force_vector_array[num][i][1],force_vector_array[num][i][2]))
                else:
                    lw.write(" maximum force:{: >9.6f}, threshold:{: >9.6f} \n".format(
                             np.max(np.abs(force_vector_array[num][i])),conver_force))

            lw.write("\n --- NEB Force vector \n")
            Coeff = 1E20 * 627.509 * 4184 * 1000 / 0.529177
            for num in range(len(force_vector_array)):
                neb_force_array[num] /= Coeff
                lw.write(" \n {} {} structurenum:{}\n".format(
                         chargespin[0],chargespin[1],num))
                if forceprint in ["on","yes"]:
                    for i in range(len(bandarray[0])):
                        lw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                        atoms[i],neb_force_array[num][i][0],neb_force_array[num][i][1],neb_force_array[num][i][2]))
                else:   
                    lw.write(" maximum force {: >9.6f}, threshold:{: >9.6f} \n".format(
                        np.max(np.abs(neb_force_array[num][i])),conver_force))
            
    if velocitysave == "on":
        with open(origdir+"Logfiles/neb_velocity_log","a") as lw:
            lw.write("\n --- Velocity --- Cyclenum:{} ----\n".format(cyclenum))
            for num in range(len(velocityarray)):
                lw.write(" {}\n {} {} structurenum:{}\n".format(len(velocityarray[0]),chargespin[0],chargespin[1],num))
                for i in range(len(velocityarray[0])):
                    lw.write(" {:<2}  {:.6g}  {:.6g}  {:.6g} \n".format(
                    atoms[i],velocityarray[num][i][0],velocityarray[num][i][1],velocityarray[num][i][2]))

            
def evaluate(newbandarray,bandarray,neb_force_array):
    inputfilename = neb_reader.get_inputfilename()
    conver_coordinate = neb_reader.input_parameterread(inputfilename,"conver_coordinate")
    conver_force = neb_reader.input_parameterread(inputfilename,"conver_force")

    tempval = 1
    failnum = 0
    
    # 座標の評価
    for i in range(len(bandarray)): 
        if np.max(np.abs(newbandarray[i]-bandarray[i])) > conver_coordinate:
            tempval *= 0
            failnum += 1
        else:
            tempval *= 1

    # 力の評価
    Coeff = 1E20 * 627.509 * 4184 * 1000 / 0.529177
    for i in range(len(neb_force_array)):
        if np.max(np.abs(neb_force_array[i])) > conver_force * Coeff:
            tempval *= 0
            failnum += 1
        else:
            tempval *= 1

    return tempval,failnum

# 構造を分裂させる
def splitting(bandarray,velocityarray,cyclenum):
    inputfilename = neb_reader.get_inputfilename()
    bandgen_splitting_number = neb_reader.input_parameterread(inputfilename,"bandgen_splitting_number")
    bandgen_splitting_period = neb_reader.input_parameterread(inputfilename,"bandgen_splitting_period")
    newbandarray = []
    newvelocityarray = []
    if cyclenum < bandgen_splitting_period * bandgen_splitting_number:
        halfnum = int(len(bandarray) / 2)
        band_toadd1 = 0.99 * bandarray[halfnum-1] + 0.01 * bandarray[halfnum]
        band_toadd2 = 0.01 * bandarray[halfnum-1] + 0.99 * bandarray[halfnum]

        for i in range(len(bandarray)):
            newbandarray.append(bandarray[i])
            newvelocityarray.append(velocityarray[i])
            if i == halfnum-1:
                newbandarray.append(band_toadd1)
                newbandarray.append(band_toadd2)
                newvelocityarray.append(velocityarray[halfnum-1])
                newvelocityarray.append(velocityarray[halfnum])

        return newbandarray,newvelocityarray

    else:
        return bandarray,velocityarray


# g16計算をするだけの関数
def forcecalc_run(forcecalc_inputfilename,forcecalc_outputfilename):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    origdir = os.getcwd()
    if origdir[-1] != "/":
        origdir = origdir + "/"

    scratchdir = neb_reader.get_scratchdirpath()
    shutil.copy(origdir + forcecalc_inputfilename, scratchdir + forcecalc_inputfilename)

    # 計算を行うときのみスクラッチファイル生成ディレクトリに移動する
    os.chdir(scratchdir)

    if program in ["none","g16","gaussian"]:
        subprocess.call(g16root+"g16 < "+ scratchdir + forcecalc_inputfilename \
                   + " > " + scratchdir + forcecalc_outputfilename, shell=True)
    else:
        subprocess.call(mopacroot+mopacprogram+" "+scratchdir+forcecalc_inputfilename,shell=True)

    os.chdir(origdir)
    shutil.copy(scratchdir + forcecalc_outputfilename, origdir + forcecalc_outputfilename)



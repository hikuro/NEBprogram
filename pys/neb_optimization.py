# neb_optimization.py
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
import glob
import datetime

# Original Modules
import gausslog
import neb_reader
import neb_function
import neb_defaultvalue

### input paths and filenames ###
# The path must end with slash "/"
g16root      = "/local/apl/lx/g16a03/g16/"                  # CHECK HERE
#programdir   = "/home/users/to1/NEBprog/"                 # CHECK HERE
g16env       = g16root + "bsd/g16.profile"        # CHECK HERE

try:
    subprocess.call("source " + g16env,shell=True)
except:
    print("g16env set failed")

mopacroot = "/Users/itsuu/MOPAC"
mopacprogram = "MOPAC2016.exe"

origdir = os.getcwd()
if origdir[-1] != "/":
    origdir += "/"

# 基本的なループ
# Logfiles/bandarrayを読み込んで最適化する
def optimization_process():
    inputfilename = neb_reader.get_inputfilename()
    # bandarrayには最新の構造情報が格納されているので読み込む
    atoms,chargespin,bandarray,energyarray = neb_reader.bandarray_read("terminated")

    if os.path.exists(origdir+"Logfiles/restartflag") == True:
        with open(origdir+"Logfiles/restartflag","r") as rr:
            progress = rr.read()

    #極値を見つけて計算する
    maxima = []
    minima = []
    for i in range(2,len(energyarray)-2):
        if energyarray[i-1] < energyarray[i] and energyarray[i+1] < energyarray[i]:
            maxima.append(i)
        elif energyarray[i-1] > energyarray[i] and energyarray[i+1] > energyarray[i]:
            minima.append(i)
    
    print("maxima:",maxima,"minima:",minima)

    #極小値の構造を最適化
    if "optdone" not in progress:
        commake(minima,maxima,atoms,chargespin,bandarray) 
        optimization(keyword="geoopt")
        with open(origdir+"Logfiles/restartflag","a") as w:
            w.write("\noptdone")

    #IRC計算を試み、次いで最適化
    if "ircdone" not in progress:
        commake_irc(maxima) 
        optimization(keyword="irc")
        with open(origdir+"Logfiles/restartflag","a") as w:
            w.write("\nircdone")

    if "footdone" not in progress:
        commake_footopt(maxima)
        optimization(keyword="footopt")
        with open(origdir+"Logfiles/restartflag","a") as w:
            w.write("\nfootdone")


### 以下、関数
# 力の定数、エネルギーをそれぞれNx3次元のリストにまとめて返す関数
def optimization(keyword):
    inputfilename = neb_reader.get_inputfilename()
    scratchdir = neb_reader.get_scratchdirpath()
    inputfile_names = glob.glob("{}*com".format(keyword))

    for filename in inputfile_names:
        g16_inputfilename = copy.copy(filename)
        g16_outputfilename = filename.replace(".com",".out")
        optcalc_run(g16_inputfilename,g16_outputfilename)

                        

# Gaussianのみ使用可能
def commake(minima,maxima,atoms,chargespin,bandarray):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    nproc = neb_reader.input_parameterread(inputfilename,"optimize_nproc")
    mem = neb_reader.input_parameterread(inputfilename,"optimize_mem") 
    connectivity = neb_reader.connectivity_read()

    theory = neb_reader.input_parameterread(inputfilename,"optimize_theory")
    option = neb_reader.input_parameterread(inputfilename,"optimize_option").split("___")

    for i,structure in enumerate([bandarray[0],bandarray[-1]]):
        if program in ["none","g16","gaussian"]:
            g16_inputfilename = "geoopt_terminal_{}.com".format(i)
            with open(origdir+g16_inputfilename,"w") as g16w:
                g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
                g16w.write(" #p {} opt freq=hpmodes\n".format(theory))
                if option != [""] and option != "":
                    for i in range(len(option)):
                        g16w.write(option[i]+" ")
                    g16w.write("\n")
                g16w.write("\n terminal number:{} \n\n".format(i))
                g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
                for i in range(len(atoms)):
                    g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                                atoms[i],structure[i][0],structure[i][1],structure[i][2]))
                
                g16w.write("\n")
                g16w.write("\n")


    for i,strnum in enumerate(minima):
        if program in ["none","g16","gaussian"]:
            g16_inputfilename = "geoopt_minimum_{}.com".format(i)
            structure = bandarray[strnum]
            with open(origdir+g16_inputfilename,"w") as g16w:
                g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
                g16w.write(" #p {} opt freq=hpmodes\n".format(theory))
                if option != [""] and option != "":
                    for i in range(len(option)):
                        g16w.write(option[i]+" ")
                    g16w.write("\n")
                g16w.write("\n minimum number:{} bandarray number:{}\n\n".format(i,strnum))
                g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
                for i in range(len(atoms)):
                    g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                                atoms[i],structure[i][0],structure[i][1],structure[i][2]))

                g16w.write("\n")
                g16w.write("\n")


    for i,strnum in enumerate(maxima):
        if program in ["none","g16","gaussian"]:
            g16_inputfilename = "geoopt_maximum_{}.com".format(strnum)
            with open(origdir+g16_inputfilename,"w") as g16w:
                g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
                g16w.write(" #p {} opt=(qst3,calcall) freq=hpmodes  \n".format(theory))
                if option != [""] and option != "":
                    for i in range(len(option)):
                        g16w.write(option[i]+" ")
                    g16w.write("\n")

                for structure in [bandarray[0],bandarray[-1],bandarray[strnum]]:
                    g16w.write("\n maximun number:{} bandarray number:{}\n\n".format(i,strnum))
                    g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
                    for i in range(len(atoms)):
                        g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                                    atoms[i],structure[i][0],structure[i][1],structure[i][2]))
                    
                g16w.write("\n\n")

# Gaussianのみ使用可能
def commake_irc(maxima):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    nproc = neb_reader.input_parameterread(inputfilename,"optimize_nproc")
    mem = neb_reader.input_parameterread(inputfilename,"optimize_mem")

    theory = neb_reader.input_parameterread(inputfilename,"optimize_theory")
    option = neb_reader.input_parameterread(inputfilename,"optimize_option").split("___")
    for i,strnum in enumerate(maxima):
        if program in ["none","g16","gaussian"]:
            g16_opt_outputfilename = "geoopt_maximum_{}.out".format(strnum)
            chargespin,electronic_energy,atoms,structure,opt_status,Conversion = (
                    gausslog.outputreader(g16_opt_outputfilename,option=0))

            if opt_status == "Optimized":
                for direction in ["forward","reverse"]:
                    g16_irc_inputfilename = "irc{}_maximum_{}.com".format(direction,strnum)
                    with open(origdir+g16_irc_inputfilename,"w") as g16w:
                        g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
                        g16w.write(" #p {} irc=({},calcall,stepsize=3)  \n".format(theory,direction))
                        if option != [""] and option != "":
                            for i in range(len(option)):
                                g16w.write(option[i]+" ")
                            g16w.write("\n")

                        g16w.write("\n maximun number:{} bandarray number:{}\n\n".format(i,strnum))
                        g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
                        for i in range(len(atoms)):
                            g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                                        atoms[i],structure[i][0],structure[i][1],structure[i][2]))

                        g16w.write("\n\n")

            else:
                print(" Transition State Not converged")
                sys.exit(1)

# Gaussianのみ使用可能
def commake_footopt(maxima):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    nproc = neb_reader.input_parameterread(inputfilename,"optimize_nproc")
    mem = neb_reader.input_parameterread(inputfilename,"optimize_mem")
    theory = neb_reader.input_parameterread(inputfilename,"optimize_theory")
    option = neb_reader.input_parameterread(inputfilename,"optimize_option").split("___")


    for i,strnum in enumerate(maxima):
        if program in ["none","g16","gaussian"]:
            for direction in ["forward","reverse"]:
                g16_irc_outputfilename = "irc{}_maximum_{}.out".format(direction,strnum)
                chargespin,electronic_energy,atoms,structure,Pt = (
                    gausslog.IRCreader(g16_irc_outputfilename))
                g16_footopt_inputfilename = "footopt{}_maximum_{}.com".format(direction,strnum)
                with open(origdir+g16_footopt_inputfilename,"w") as g16w:
                    g16w.write(" %nproc={}\n %mem={}\n".format(nproc,mem))
                    g16w.write(" #p {} opt=calcall freq=hpmodes  \n".format(theory))
                    if option != [""] and option != "":
                        for i in range(len(option)):
                            g16w.write(option[i]+" ")
                        g16w.write("\n")

                    g16w.write("\n maximun number:{} bandarray number:{}\n\n".format(i,strnum))
                    g16w.write(" {} {}\n".format(chargespin[0],chargespin[1]))
                    for i in range(len(atoms)):
                        g16w.write(" {} {:>9.6f} {:>9.6f} {:>9.6f}\n".format(
                                    atoms[i],structure[i][0],structure[i][1],structure[i][2]))

                    g16w.write("\n\n")



# g16計算をするだけの関数
def optcalc_run(optcalc_inputfilename,optcalc_outputfilename):
    inputfilename = neb_reader.get_inputfilename()
    program = neb_reader.input_parameterread(inputfilename,"program")

    origdir = os.getcwd()
    if origdir[-1] != "/":
        origdir = origdir + "/"

    scratchdir = neb_reader.get_scratchdirpath()
    shutil.copy(origdir + optcalc_inputfilename, scratchdir + optcalc_inputfilename)

    # 計算を行うときのみスクラッチファイル生成ディレクトリに移動する
    os.chdir(scratchdir)

    if program in ["none","g16","gaussian"]:
        subprocess.call(g16root+"g16 < "+ scratchdir + optcalc_inputfilename \
                   + " > " + scratchdir + optcalc_outputfilename, shell=True)
    else:
        subprocess.call(mopacroot+mopacprogram+" "+scratchdir+optcalc_inputfilename,shell=True)

    os.chdir(origdir)
    shutil.copy(scratchdir + optcalc_outputfilename, origdir + optcalc_outputfilename)



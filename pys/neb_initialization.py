# neb_initialization.py
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
import datetime

# Original modules
import neb_reader
import neb_function

### input paths and filenames ###
# The path must end with slash "/"
origdir = os.getcwd()
if origdir[-1] != "/":
    origdir += "/"

#基本的な操作
def initialize():
    inputfilename = neb_reader.get_inputfilename()
    bandgen_method = neb_reader.input_parameterread(inputfilename,"bandgen_method")
    bandgen_splitting = neb_reader.input_parameterread(inputfilename,"bandgen_splitting")
    # 初期構造の読み込み
    terminal,transition,atoms,chargespin = neb_reader.input_structureread()

    # 初期構造を回転させterminal1に合わせ、
    # terminal0,transition0–n,terminal1の順番に並べる(Procrustes Superimposition)
    originalbandarray = neb_function.rotstructures_reorder(atoms,terminal,transition)

    # 初期構造生成の方法により分岐する
    # none: 単純に座標の重み付けをして線形結合
    # dihedral: 目立って変化している二面角をScanで変化させて作成
    # Bandの生成、書き込み、Logfileへの書き込み
    if bandgen_splitting in ["on","yes"]:
        bandarray = [originalbandarray[0],originalbandarray[-1]]
        bandarraywrite(atoms,chargespin,bandarray)

    elif bandgen_method in ["none","near_edge","nearedge"]:
        bandarray = neb_function.bandgen(originalbandarray)
        bandarraywrite(atoms,chargespin,bandarray)

    elif bandgen_method in ["dihedral"]:
        import neb_bandgen_dihedral
        bandarray = neb_bandgen_dihedral.bandgen(originalbandarray,atoms)
        bandarraywrite_multiple(atoms,chargespin,bandarray)
        initialize_resultwrite_multiple(atoms,chargespin,bandarray)

    elif bandgen_method in ["read"]:
        if os.path.exists(origdir+"Logfiles/bandarray") == False:
            print("Initial Bandarray structures are required!")
            sys.exit()

### 以下、関数

def initialize_resultwrite(atoms,chargespin,bandarray):
    time_now = datetime.datetime.now()
    time_now_str = time_now.strftime("%Y/%m/%d/ %H:%M")
    
    with open(origdir+"Logfiles/neb_log","a") as nl:
        nl.write(" NEBprog Logfile starttime: {}\n\n".format(time_now_str))
        nl.write(" Process:1 Initialization Done \n")
        nl.write(" -- Band coordinates --\n")
        for i in range(len(bandarray)):
            nl.write(" {}\n {} {}\n".format(len(bandarray[0]),chargespin[0],chargespin[1]))
            for j in range(len(bandarray[0])):
                nl.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f}\n".format(
                         atoms[j],bandarray[i][j][0],bandarray[i][j][1],bandarray[i][j][2]))
 
def bandarraywrite(atoms,chargespin,bandarray):
    with open(origdir+"Logfiles/bandarray","w") as bw:
        for num in range(len(bandarray)):
            bw.write("{}\n{} {} structurenum:{}\n".format(len(bandarray[0]),chargespin[0],chargespin[1],num))
            for i in range(len(bandarray[0])): 
                bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                atoms[i],bandarray[num][i][0],bandarray[num][i][1],bandarray[num][i][2]))

    with open(origdir+"Logfiles/bandarray_init","w") as bw:
        for num in range(len(bandarray)):
            bw.write("{}\n{} {} structurenum:{}\n".format(len(bandarray[0]),chargespin[0],chargespin[1],num))
            for i in range(len(bandarray[0])):
                bw.write(" {:<2} {: >9.6f} {: >9.6f} {: >9.6f} \n".format(
                atoms[i],bandarray[num][i][0],bandarray[num][i][1],bandarray[num][i][2]))


def writeparameter():
    cineb = neb_reader.input_parameterread(inputfilename,"cineb")
    with open(origdir+"Logfiles/neb_log","a") as nl:
        pass

















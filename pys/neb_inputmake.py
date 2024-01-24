# インプットファイル群があるディレクトリ内で本プログラムを走らせる
# settingsにインプットファイルの情報を入力しておく
# ディレクトリ内にはgjf,com,out,logのいずれかの拡張子のみをもつ構造ファイルを用意する。
import sys
import os
import glob
import numpy as np

import gausslog
import neb_function
import neb_reader

rundir = os.getcwd()
if rundir[-1] != "/":
    rundir += "/"

# str2のメチル基を検出してstr1との誤差が最小になるよう回転させる
def rotatemethyl(str1,str2,atoms):
    Adjacency_matrix = neb_function.definebonds(atoms,str1,read=0)
    methylhydrogen = get_methylhydrogen(Adjacency_matrix,atoms)

    if len(methylhydrogen) != 0:
        str2_modified = neb_function.rotation_optimization(str2,str1,velocity=0) 
        for Hset in methylhydrogen:
            exchange_H = [[i for i in range(len(str2))],[i for i in range(len(str2))]]
            exchange_H[0][Hset[0]] = Hset[1]
            exchange_H[0][Hset[1]] = Hset[2]
            exchange_H[0][Hset[2]] = Hset[0]
            exchange_H[1][Hset[0]] = Hset[2]
            exchange_H[1][Hset[1]] = Hset[0]
            exchange_H[1][Hset[2]] = Hset[1]
            #比較の前に更にrotation_optimizationでフィッティングしたほうが本来は正確
            dev0 = np.std(str2_modified-str1)
            dev1 = np.std(str2_modified[exchange_H[0],:]-str1)
            dev2 = np.std(str2_modified[exchange_H[1],:]-str1)
            if dev0 >= dev1 >= dev2 or dev1 >= dev0 >= dev2:
                str2 = str2_modified[exchange_H[1],:]
                print("Rotated")
            elif dev0 >= dev2 >= dev1 or dev2 >= dev0 >= dev1:
                str2 = str2_modified[exchange_H[0],:]
                print("Rotated2")

    return str2

def get_methylhydrogen(Adjacency_matrix,atoms):
    # First, get methyl carbon
    methyl_H = []
    for i in range(len(Adjacency_matrix)):
        num_H = 0
        if atoms[i] == "C" and np.sum(Adjacency_matrix[i]) == 4:
            temp_H = []
            for j in range(len(Adjacency_matrix)):
                if Adjacency_matrix[i][j] == 1 and atoms[j] == "H":
                    num_H += 1
                    temp_H.append(j)
            if num_H == 3:
                methyl_H.append(temp_H)

    return methyl_H


def inputwrite(inputfilename,str1,str2,atoms,chargespin):
    with open("settings","r") as r:
        settinginfo = r.read()

    if "fitmethyl" in settinginfo:
        print("trying fitmethyl:",inputfilename)
        str2 = rotatemethyl(str1,str2,atoms)

    with open(inputfilename,"w") as w:
        w.write(settinginfo)
        w.write("\n\ninput terminal1\n{} {}\n".format(chargespin[0],chargespin[1]))
        for i in range(len(atoms)):
            w.write(" {:<2} ".format(atoms[i]))
            for j in range(3):
                w.write(" {: >9.6f} ".format(str1[i][j]))
            w.write("\n")

        w.write("end\n\n")
        w.write("\n\ninput terminal2\n{} {}\n".format(chargespin[0],chargespin[1]))
        for i in range(len(atoms)):
            w.write(" {:<2} ".format(atoms[i]))
            for j in range(3):
                w.write(" {: >9.6f} ".format(str2[i][j]))
            w.write("\n")

        w.write("end\n\n")

def structurewrite(structure,atoms,chargespin):
    with open(rundir+"inputstructures","w") as s:
        for coordinates in structure:
            s.write("{}\n{} {}\n".format(len(atoms),chargespin[0],chargespin[1]))
            for i,atom_coordinate in enumerate(coordinates):
                s.write(" {:<2} ".format(atoms[i]))
                for j in range(3):
                    s.write(" {: >9.6f} ".format(atom_coordinate[j]))
                s.write("\n")

def make_submissionfile(directoryname):
    #multiprocessing_totalcoreの情報を取得し記入する
    multiprocessing_totalcore = neb_reader.input_parameterread(rundir+"settings","multiprocessing_totalcore")
    with open(directoryname+"/runneb.csh","w") as s:
        s.write("#!/bin/csh -f\n")
        s.write("#PBS -l select=1:ncpus={}".format(int(multiprocessing_totalcore)))
        s.write(":mpiprocs=1:ompthreads=1:jobtype=core\n")
        s.write("#PBS -l walltime=72:00:00\n\n")
        s.write("module purge\n\n")
        s.write("if ($?PBS_O_WORKDIR) then\n")
        s.write("cd ${PBS_O_WORKDIR}\n")
        s.write("endif\n\n")
        s.write("source /local/apl/lx/g16a03/g16/bsd/g16.login\n")
        s.write("mkdir /work/users/to1/tmp.$$\n")
        s.write("set WORK=/work/users/to1/tmp.$$\n\n")
        s.write("python3 /home/users/to1/NEBprog/neb_starter.py neb.com ${WORK}\n")
        s.write("exit 0\n")

if __name__=="__main__":
    files = (glob.glob("*"))
    files.sort()
    filetype = "none"
    print("rindir:",rundir)
    print(files)

    try:
        files.remove("settings")
    except:
        print("No setting file")
        sys.exit()

    try:
        files.remove("inputstructures")
    except:
        pass

    if "com" in files[0] or "gjf" in files[0]:
        filetype = "com"
    elif "log" in files[0] or "out" in files[0]:
        filetype = "out"

    structure = []
    atoms = []
    chargespin = [0,0]
    for filename in files:
        if filetype == "com":
            temp_structure,atoms,chargespin = gausslog.inputfilereader(filename)
        elif filetype == "out":
            temp_structure,atoms,dummy = gausslog.structurereader(filename)
            chargespin = gausslog.chargespinreader(filename)

        structure.append(temp_structure)
        
#    print(structure)
    if chargespin[0] == 0 and chargespin[1] == 0:
        print("Bad chargespin")
        sys.exit()

    # 読み込んだ後はfile内の構造を一応書き出しておく
    structurewrite(structure,atoms,chargespin)    
    
    # 構造同士で組を作り、ディレクトリ作成。NC2個のディレクトリが作成される。
    for i in range(len(structure)-1):
        for j in range(i+1,len(structure)):
            directoryname = rundir + "transition{}-{}".format(str(i),str(j))
            os.mkdir(directoryname)
            inputwrite(directoryname+"/neb.com",structure[i],structure[j],atoms,chargespin)
            make_submissionfile(directoryname)







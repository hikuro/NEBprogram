# neb_function.py
# ver 0.5 2020/09/14 
# ver 0.6 2020/10/26

# Standard library
import os
import sys
import copy
import numpy as np
import shutil
import subprocess
import linecache

# Original modules
import neb_reader

### input paths and filenames ###
# The path must end with slash "/"
atom_to_weight = {'H':1.00783,'He':4.0026,'Li':6.941,'Be':9.012,'B':10.811,
            'C':12.,'N':14.007,'O':15.9994,'F':18.9984,'Ne':20.1797,
            'Na':22.989,'Mg':24.305,'Al':26.98154,'Si':28.0855,'P':30.9738,
            'S':32.066,'Cl':35.4527,'Ar':39.948,'K':39.0983,'Ca':40.078,
            'Sc':44.96,'Ti':47.867,'V':50.94,'Cr':51.9961,'Mn':54.938,
            'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
            'Ga':69.723,'Ge':72.64,'As':74.9216,'Se':78.96,'Br':79.904,
            'Pd':106.42,'I':126.90447}
### input paths and filenames ###
# The path must end with slash "/"

# Set atom sets 
period1 = {'H','He'}
period2 = {'Li','B','C','N','O','F','Ne'}
period3 = {'Na','Mg','Al','Si','P','S','Cl'}
period4 = {'K','Ca','Ti','Br'}
period5 = {'Zr','Sn','Sb','I'}

origdir = os.getcwd()
if origdir[-1] != "/":
    origdir += "/"

def rotstructures_reorder(atoms,terminal,transition):
    newterminal,newtransition = rotstructures(atoms,terminal,transition)
    bandarray = []
    bandarray.append(newterminal[0])
    for i in range(len(transition)):
        bandarray.append(transition[i])
    bandarray.append(newterminal[1])

    return bandarray

def rotstructures(atoms,terminal,transition):
    newterminal = []
    newtransition = []
    for structure in terminal:
        newterminal.append(rotation_optimization(structure,terminal[0]))
    if len(transition) > 0:
        for structure in transition:
            newtransition.append(rotation_optimization(structure,terminal[0]))

    return newterminal,newtransition

# str_to_rotate(Ri)、str_standard(Rs)、回転行列(A)があったとき、
# |ARi - Rs| を最小化する回転行列Aを求めてARiを返す関数
def rotation_optimization(str_to_rotate_raw,str_standard_raw,velocity=0):
    #重心を原点に合わせる
    str_to_rotate = copy.deepcopy(str_to_rotate_raw)
    str_standard = copy.deepcopy(str_standard_raw)

    str_standard -= np.sum((str_standard),axis=0)/ len(str_standard)
    str_to_rotate -=  np.sum((str_to_rotate),axis=0)/ len(str_standard)

    # 焼き鈍し法に準ずるが不安定状態への遷移確率は０とする（経験的に問題ないため）
    # 特異値分解を使っても同じ結果が得られる
    rotmat = np.identity(3)
    cycle = 1
    cool = 0.9998
    T = 10
    cyclelimit = 1000
    Dev = np.std(str_to_rotate - str_standard)
    while cycle < cyclelimit:
        Rot = func_Rotation((np.random.rand(3) -0.5) * T / cycle)
        newRot = np.dot(Rot,rotmat)
        newDev = np.std(np.dot(newRot,str_to_rotate.T).T - str_standard)
        if newDev < Dev:
            Dev = newDev
            rotmat = newRot
        T *= cool
        cycle += 1

    structure = np.dot(rotmat,str_to_rotate.T).T

    if type(velocity) is np.ndarray:
        newvelocity = np.dot(rotmat,velocity.T).T
        return structure, newvelocity

    elif type(velocity) is int:
        return structure

def func_Rotation(Var):
    RotX = np.array([[1,0,0],[0,np.cos(Var[0]),-np.sin(Var[0])],[0,np.sin(Var[0]),np.cos(Var[0])]])
    RotY = np.array([[np.cos(Var[1]),0,np.sin(Var[1])],[0,1,0],[-np.sin(Var[1]),0,np.cos(Var[1])]])
    RotZ = np.array([[np.cos(Var[2]),-np.sin(Var[2]),0],[np.sin(Var[2]),np.cos(Var[2]),0],[0,0,1]])
    Rot = np.dot(np.dot(RotX,RotY),RotZ)
    return Rot

# bandの端点と中間構造を指定することでbandの要素を生成する
def bandgen(originalbandarray):
    inputfilename = neb_reader.get_inputfilename()
    bandpointnum = neb_reader.input_parameterread(inputfilename,"bandpointnum")
    bandgen_method = neb_reader.input_parameterread(inputfilename,"bandgen_method")

    devidenum = 1 + int((bandpointnum - len(originalbandarray)) / (len(originalbandarray)-1))
    #print("bandgen, originalbandarray",originalbandarray)
    bandarray_temp = []
    if bandgen_method == "none":
        for i in range(len(originalbandarray)-1):
            for j in range(devidenum):
                bandarray_temp.append((originalbandarray[i] 
                                  + j/devidenum *(originalbandarray[i+1] - originalbandarray[i])))
    elif bandgen_method in ["near_edge","nearedge"]:
        for j in range(int(devidenum/2)):
            bandarray_temp.append((originalbandarray[0]
                              + j/(10*devidenum) *(originalbandarray[-1] - originalbandarray[0])))
        for j in range(int(devidenum/2)):
            bandarray_temp.append((originalbandarray[-1]
                              + (int(devidenum/2)- j)/(10*devidenum) *(originalbandarray[0] - originalbandarray[-1])))

    bandarray_temp.append(originalbandarray[-1])
    bandarray = procrustes_superimposition(bandarray_temp)

    return bandarray

# i >= 1 のとき、i=0~N-1のN個の構造で構成されるbandの要素RiをRi-1の構造に対して重ねる
def procrustes_superimposition(bandarray,velocityarray=0):
    new_bandarray = []
    new_bandarray.append(bandarray[0])
    if type(velocityarray) is int:
        for i in range(1,len(bandarray)):
            temp_array = rotation_optimization(bandarray[i],new_bandarray[-1],0)
            new_bandarray.append(temp_array)

        new_bandarray = nancheck_band(new_bandarray)
        return new_bandarray

    elif type(velocityarray) is list:
        new_velocityarray = []
        new_velocityarray.append(velocityarray[0])
        for i in range(1,len(bandarray)):
            temp_bandarray,temp_velocityarray = rotation_optimization(bandarray[i],
                                                 new_bandarray[i-1],velocityarray[i])
            new_bandarray.append(temp_bandarray)
            new_velocityarray.append(temp_velocityarray)

        new_bandarray = nancheck_band(new_bandarray)
        new_velocityarray = nancheck_velocity(new_velocityarray)
        return new_bandarray,new_velocityarray

#エラーが起きた場合にリカバリーするための関数
def nancheck_band(bandarray):
    for i in range(len(bandarray)):
        if np.isnan(bandarray[i][0][0]) == True:
            small_terminal = 0
            large_terminal = 0 
            for j in range(0,i):
                if np.isnan(bandarray[i][0][0]) == True:
                    small_terminal = j - 1
                    break
            for k in range(i+1,len(bandarray)):
                if np.isnan(bandarray[i][0][0]) == False:
                    large_terminal = k
                    break
            bandarray[i] = (( (large_terminal - i) * bandarray[small_terminal] 
                            + (i - small_terminal) * bandarray[large_terminal] ) 
                            / (large_terminal - small_terminal))

    return bandarray

def nancheck_band2(bandarray,energy_array):
    drag_threshold = neb_reader.input_parameterread(inputfilename,"drag_threshold")

    for i in range(1,len(bandarray)-1):
        if energy_array[i] == 0:
            if energy_array[i-1] != 0:
                if i > int(len(bandarray)) and energy_array[i+1] != 0:
                    bandarray[i] = 0.5 * bandarray[i+1] + 0.5 * bandarray[i]
                else:
                    bandarray[i] = 0.5 * bandarray[i-1] + 0.5 * bandarray[i]
            elif energy_array[i+1] != 0:
                bandarray[i] = 0.5 * bandarray[i+1] + 0.5 * bandarray[i]

    return bandarray

def nancheck_velocity(velocityarray):
    for i in range(len(velocityarray)):
        velocityarray[i] = np.nan_to_num(velocityarray[i])

    return velocityarray

# 水素を省いた隣接行列の作成
def definebonds(atomarray,coordinates,read=1):
    if read == 1:
        inputfilename = neb_reader.get_inputfilename()
        dihedral_force_on_H = neb_reader.input_parameterread(inputfilename,"dihedral_force_on_H")
    else:
        dihedral_force_on_H = "on"

    Adjacency_matrix = np.zeros([len(coordinates),
                                 len(coordinates)])

    BondLength = 0.0
    for i in range(len(Adjacency_matrix)):
        for j in range(i+1,len(Adjacency_matrix)):
            BondLength = distance(coordinates[i],coordinates[j])
            if (atomarray[i] in period2 and atomarray[j] in period2):
                if  BondLength - 1.8 < 0: # So as to calculate TS
                    Adjacency_matrix[i][j] = Adjacency_matrix[j][i] = 1

            elif ((atomarray[i] in period2 and atomarray[j] in period3)
              or (atomarray[j] in period2 and atomarray[i] in period3)):
                if BondLength - 2.4 < 0:
                    Adjacency_matrix[i][j] = Adjacency_matrix[j][i] = 1

            elif ((atomarray[i] in period2 and atomarray[j] in period4)
              or (atomarray[j] in period2 and atomarray[i] in period4)):
                if BondLength - 2.4 < 0:
                    Adjacency_matrix[i][j] = Adjacency_matrix[j][i] = 1

            elif ((atomarray[i] in period2 and atomarray[j] in period5)
              or (atomarray[j] in period2 and atomarray[i] in period5)):
                if BondLength - 2.8 < 0:
                    Adjacency_matrix[i][j] = Adjacency_matrix[j][i] = 1

            elif ((atomarray[i] not in period1 and atomarray[j] == "H")
                and  dihedral_force_on_H in ["yes","on"]):
                if  BondLength - 1.6 < 0:
                    Adjacency_matrix[i][j] = Adjacency_matrix[j][i] = 1

    return Adjacency_matrix

# A function to measure a distance of atoms
def distance(atom1,atom2):
    distance_sq = 0.0
    for i in range(3):
        distance_sq += (atom1[i]-atom2[i])**2
    distance = distance_sq ** 0.5

    return distance






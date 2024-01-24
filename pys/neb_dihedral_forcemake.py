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

# Original Modules
import gausslog
import neb_reader
import neb_function

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

# Set atom sets 
period1 = {'H','He'}
period2 = {'Li','B','C','N','O','F','Ne'}
period3 = {'Na','Mg','Al','Si','P','S','Cl'}
period4 = {'K','Ca','Ti','Br'}
period5 = {'Zr','Sn','Sb','I'}

origdir = os.getcwd()
if origdir[-1] != "/":
    origdir += "/"

def forcemake(bandarray,atoms):
    inputfilename = neb_reader.get_inputfilename()
    dihedral_force_constant = neb_reader.input_parameterread(inputfilename,"dihedral_force_constant")

    # ignoreHの場合、骨格のみ考えるために水素原子を無視して
    # 構造の二面角を計算する
    # bandarrayの最初の構造で隣接行列を作成するので注意
    dihedralangle_array = dihedral_angle_get(bandarray,atoms)

    # 隣接するarray(要素N)の全ての二面角の差を計算する(N個目の要素はゼロ行列)
    dihedralangle_array_difference = dihedral_angle_compare(dihedralangle_array)
#    return dihedralangle_array,dihedralangle_array_difference

    # ねじれの力を計算する
    force_on_angle = dihedral_force_calculation(dihedralangle_array_difference,bandarray)
    force_on_angle *= dihedral_force_constant

    return force_on_angle

### 以下、関数
def dihedral_force_calculation(dihedralangle_array_difference,bandarray):
    # Dihedral Angleに基づいた力の方向の計算
    force_direction_array = np.zeros([len(bandarray),len(bandarray[0]),len(bandarray[0][0])])
    for angle_index,angles in enumerate(dihedralangle_array_difference[0]):
        #k-i-j-lの順で結合指定
        atom_k = int(angles[0][0])
        atom_i = int(angles[0][1])
        atom_j = int(angles[0][2])
        atom_l = int(angles[0][3])
        for str_index,structure in enumerate(bandarray):
            #まずatom_kへの力の計算
            vector_i_to_k = structure[atom_k] - structure[atom_i]
            vector_i_to_j = structure[atom_j] - structure[atom_i]
            angle_diff = dihedralangle_array_difference[str_index][angle_index][1]
            force_k = (np.sin(angle_diff) * np.cross(vector_i_to_k,vector_i_to_j) 
                       / (np.linalg.norm(vector_i_to_j) * np.linalg.norm(vector_i_to_k)))
#            force_k = (np.sin(angle_diff) * np.cross(vector_i_to_j,vector_i_to_k) 
#                       / (np.linalg.norm(vector_i_to_j) * np.linalg.norm(vector_i_to_k)))
            #次にatom_lへの力の計算（正規化した後に角度に応じた重み付けしている）
            vector_j_to_l = structure[atom_k] - structure[atom_i]
            force_l = (np.sin(angle_diff) * np.cross(vector_i_to_j,vector_j_to_l)
                       / (np.linalg.norm(vector_i_to_j) * np.linalg.norm(vector_j_to_l)))
#            force_l = (np.sin(angle_diff) * np.cross(vector_i_to_j,vector_i_to_k)
#                       / (np.linalg.norm(vector_i_to_j) * np.linalg.norm(vector_j_to_l)))

            #力の定数を行列に加算
            force_direction_array[str_index][atom_k] += force_k
            force_direction_array[str_index][atom_l] += force_l

    return force_direction_array


# 二面角(k-i-j-l)の解析
# リスト内に[[k,i,j,l],degree]が格納される。i<jなので一意的に二面角が指定される。
def dihedralanalyze(atomarray,coordinates,Adjacency_matrix):
    dihedrals = []
    for i in range(len(Adjacency_matrix)):
        for j in range(i+1,len(Adjacency_matrix)):
            if Adjacency_matrix[i][j] == 1 and i != j:
                for k in range(len(Adjacency_matrix)):
                    if Adjacency_matrix[i][k] == 1 and j != k and i != j:
                        for l in range(len(Adjacency_matrix)):
                            if Adjacency_matrix[j][l] == 1 and i != l and k != l:
                                Dihedral_angle = dihedral(coordinates,k,i,j,l)
                                dihedrals.append([[k,i,j,l],Dihedral_angle])
           
    return dihedrals


def dihedral_angle_get(bandarray,atoms):
    inputfilename = neb_reader.get_inputfilename()
    dihedral_force_on_H = neb_reader.input_parameterread(inputfilename,"dihedral_force_on_H")
    Adjacency_matrix = neb_function.definebonds(atoms,bandarray[0]) 

    dihedralangle_array = []
    for coordinates in bandarray:
        dihedrals = dihedralanalyze(atoms,coordinates,Adjacency_matrix)
        dihedralangle_array.append(dihedrals)

    return dihedralangle_array


# 二面角計算用関数(単位はラジアン)
def dihedral(structure,atom1,atom2,atom3,atom4):
    # set atom2 to original point
    coord1 = structure[atom1]-structure[atom2]
    coord3 = structure[atom3]-structure[atom2]
    coord4 = structure[atom4]-structure[atom2]
    # set atom 3 to X axis
    r_xyz = np.linalg.norm(coord3)
    r_xy  = np.linalg.norm(coord3[0:2])
    sint  = coord3[1] / r_xy
    cost  = coord3[0] / r_xy
    sinp  = r_xy    / r_xyz
    cosp  = coord3[2] / r_xyz
    rot_mat_yx = np.array([[ cost, sint, 0],
                           [-sint, cost, 0],
                           [    0,    0, 1]])
    rot_mat_zx = np.array([[ sinp, 0, cosp],
                           [    0, 1,    0],
                           [-cosp, 0, sinp]])
    coord1 = np.dot(rot_mat_zx,np.dot(rot_mat_yx,coord1.T)).T
    coord4 = np.dot(rot_mat_zx,np.dot(rot_mat_yx,coord4.T)).T
    coord1[0] = coord4[0] = 0.0
    # dihedral angleの定義
    # np.sign(np.cross(coord1,coord4)[0])をかけるまでは0<dihedral<piの間
    # coord2->coord3(正)の方向に関してcoord4がcoord1より右回り側にあれば正値、逆なら負値
    dihedral = (np.arccos(np.dot(coord1[1:],coord4[1:])
                /(np.linalg.norm(coord1[1:])*np.linalg.norm(coord4[1:])))
                * np.sign(np.cross(coord1,coord4)[0]) )

    return dihedral

# i-1, i, i+1番目の3つの構造のdihedral angle array からi+1, i-1とiとの差の差を計算して出力する
def dihedral_angle_compare(dihedralangle_array):
    nullinfo = []
    for i in range(len(dihedralangle_array[0])):
        nullinfo.append([dihedralangle_array[0][i][0],0])

    dihedralangle_array_difference = [nullinfo]
    for i in range(1,len(dihedralangle_array)-1):
        temp = []
        for anglenum in range(len(dihedralangle_array[0])):
            diff = (dihedralangle_array[i+1][anglenum][1] + dihedralangle_array[i-1][anglenum][1]
                    - 2 * dihedralangle_array[i][anglenum][1])
            dihedralangle_array_difference_temp = [dihedralangle_array[i][anglenum][0],diff]
            temp.append(dihedralangle_array_difference_temp)
        dihedralangle_array_difference.append(temp)

    dihedralangle_array_difference.append(nullinfo)

    return  dihedralangle_array_difference



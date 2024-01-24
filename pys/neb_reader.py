# neb_reader.py
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
import neb_defaultvalue


origdir = os.getcwd()
if origdir[-1] != "/":
    origdir = origdir + "/"

def get_scratchdirpath():
    with open(origdir+"Logfiles/scratchdir","r") as w:
        val = w.read()
    return val

def get_inputfilename():
    with open(origdir+"Logfiles/inputfilename","r") as w:
        val = w.read()
    return val

# This function reads inputfile
def input_structureread():
    # inputfileは、bandpointnum, cyclenum, theory1and2, 
    # terminalcoordinate1and2i,transition1~ の情報を持つ
    inputfilename = get_inputfilename()
    terminal = []
    transition = []
    atoms = []
    chargespin = []
    bandpointnum = input_parameterread(inputfilename,"bandpointnum")
    transitionnum = input_parameterread(inputfilename,"transition")

    for keyword in ["terminal1","terminal2"]:
        chargespin_temp,atoms_temp,coordinates = input_structureread_module(inputfilename,keyword)
        terminal.append(coordinates)
        atoms.append(atoms_temp)
        chargespin.append(chargespin_temp)

    for i in range(transitionnum):
        keyword = "transition"+str(i+1)
        chargespin_temp,atoms_temp,coordinates = input_structureread_module(inputfilename,keyword)
        transition.append(coordinates)
        atoms.append(atoms_temp)
        chargespin.append(chargespin_temp)

    # Check whether all atom orders are unified
    for atomorder in atoms:
        if atomorder != atoms[0]:
            print("Atom order should be the same")
            sys.exit()
    
    for chargespinvals in chargespin:
        if chargespinvals != chargespin[0]:
            print("Chargespins should be the same")
            sys.exit()

    return terminal,transition,atoms[0],chargespin[0]


def input_parameterread(inputfilename,val):
    linenum  = sum(1 for line in open(inputfilename))
    counted_number = 0
    returnval = "none"

    for i in range(1,linenum+1):
        Lineinfo = linecache.getline(inputfilename,i).split()
        try:
            if Lineinfo[0] == val:
                if Lineinfo[1] != "=":
                    returnval = Lineinfo[1]
                elif Lineinfo[1] == "=":
                    returnval = Lineinfo[2]

            elif Lineinfo[0] == "input" and val in Lineinfo[1]:
                counted_number += 1

        except:
            pass

    if val in ["transition"]:
        return counted_number

    else:
        if returnval in ["variable"]:
            return returnval

        defaultval = neb_defaultvalue.defaultvalue(val)
        if returnval == "none":
            returnval = defaultval
        else:
            if type(defaultval) is int:
                returnval = int(returnval)
            elif type(defaultval) is float:
                returnval = float(returnval)

        return returnval

def input_structureread_module(inputfilename,keyword):
    linenum  = sum(1 for line in open(inputfilename))
    flag = 0
    chargespin = []
    atoms = []
    coordinates = []
    for i in range(1,linenum+1):
        Lineinfo = linecache.getline(inputfilename,i).split()
        if len(Lineinfo) <= 1:
            Lineinfo += ["",""]

        if flag == 1:
            if Lineinfo[0] == "end":
                break
            elif Lineinfo[0].isnumeric():
                chargespin = [int(Lineinfo[0]),int(Lineinfo[1])]
            else:
                atoms.append(Lineinfo[0])
                coordinates.append([float(Lineinfo[num]) for num in range(1,4)])

        if Lineinfo[0] == "input" and Lineinfo[1] == keyword:
            flag = 1

    coordinates = np.array(coordinates)

    return chargespin,atoms,coordinates

def connectivity_read():
    inputfilename = get_inputfilename()
    linenum  = sum(1 for line in open(inputfilename))
    connectivity = []
    flag = 0
    for i in range(1,linenum+1):
        line = linecache.getline(inputfilename,i)
        Lineinfo = line.split()
        if len(Lineinfo) <= 1:
            Lineinfo += ["",""]

        if flag == 1:
            if Lineinfo[0] == "end":
                break
            else:
                connectivity.append(line)

        if Lineinfo[0] in  ["connectivity","conectivity","connect"]:
            flag = 1

    return connectivity


def bandarray_read(val=0):
    if val == 0:
        filename = origdir+"Logfiles/bandarray"
    elif val == "terminated":
        filename = origdir+"Logfiles/bandarray_terminated"

    atoms,chargespin,bandarray,energyarray = general_read(filename)
    return atoms,chargespin,bandarray,energyarray

def velocityarray_read():
    filename = origdir+"Logfiles/velocityarray"
    dummy,dummy2,velocityarray,dummy3 = general_read(filename)
    return velocityarray    

def general_read(filename):
    linenum  = sum(1 for line in open(filename))
    atomnum = int(linecache.getline(filename,1).split()[0])
    strnum = int(linenum / (atomnum + 2))
    atoms = []
    temp = linecache.getline(filename,2).split()
    chargespin = [int(temp[0]),int(temp[1])]
    array = []
    energyarray = []
    for i in range(strnum):
        tempstr = []
        for j in range(atomnum):
            Lineinfo = linecache.getline(filename,i*(atomnum+2)+j+3).split()
            if i == 0:
                atoms.append(Lineinfo[0])
            tempstr.append([float(Lineinfo[1]),float(Lineinfo[2]),float(Lineinfo[3])])
        array.append(np.array(tempstr))

    try:
        for i in range(strnum):
            Lineinfo = linecache.getline(filename,i*(atomnum+2)+2).split("energy:")
            Lineinfo = Lineinfo[1].split()
            energyarray.append(float(Lineinfo[0]))

    except:
        pass

    return atoms,chargespin,array,energyarray











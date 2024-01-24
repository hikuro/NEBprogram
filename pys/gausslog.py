#################################### gausslog.py #######################################
#
# ver 07/15/2017
# ver 06/11/2018 revised for QuasiPROGDYN
# ver 06/21/2018 oniomlogmake is implemented
# ver 09/07/2018 revised for MOLPRO2016 
# ver 10/20/2020 revised for NEBprog
#
# written for python 3.6
# This program reads gaussian09 output file and returns calculated parameters
# MOLPRO output files are also read and edited.
#
# written by Hiroaki Kurouchi
#
#########################################################################################

import sys
import numpy as np
import linecache

# define parameters
number_atom  = {
   '1':'H','2':'He','3':'Li','4':'Be','5':'B',
   '6':'C','7':'N','8':'O','9':'F','10':'Ne',
   '11':'Na','12':'Mg','13':'Al','14':'Si','15':'P',
   '16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca',
   '21':'Sc','22':'Ti','23':'V','24':'Cr','25':'Mn',
   '26':'Fe','27':'Co','28':'Ni','29':'Cu','30':'Zn',
   '31':'Ga','32':'Ge','33':'As','34':'Se','35':'Br',
   '36':'Kr','37':'Rb','38':'Sr','39':'Y','40':'Zr',
   '41':'Nb','42':'Mo','43':'Tc','44':'Ru','45':'Rh',
   '46':'Pd','47':'Ag','48':'Cd','49':'In','50':'Sn',
   '51':'Sb','52':'Te','53':'I','54':'Xe','55':'Cs','56':'Ba'}

atom_weight = {'H':1.00784,'He':4.0026,'Li':6.941,'Be':9.012,'B':10.811,
               'C':12.0,'N':14.007,'O':15.9994,'F':18.9984,'Ne':20.1797,
               'Na':22.989,'Mg':24.305,'Al':26.98154,'Si':28.0855,'P':30.9738,
               'S':32.066,'Cl':35.4527,'Ar':39.948,'K':39.0983,'Ca':40.078,
               'Sc':44.96,'Ti':47.867,'V':50.94,'Cr':51.9961,'Mn':54.938,
               'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
               'Ga':69.723,'Ge':72.64,'As':74.9216,'Se':78.96,'Br':79.904,
               'Pd':106.42,'I':126.90447}


def titlereader(filename,title1):
    linenumtot  = sum(1 for line in open(filename))
    linenum = 1
    title = "Error for some reason :-<"
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) == 0:
            lineinfo = "empty"
        if lineinfo[0] == title1:
            title = lineinfo
            break
        linenum += 1
    linecache.clearcache()
    return title

def structurereader(filename):
    # this function reads standard orientation of the molecule.
    # the last structure in the logfile is read
    Natoms      = atomnum(filename)
    linenumtot  = sum(1 for line in open(filename))
    linenum = 1
    flagnum = 0
    linenum = 1
    thermoflag = 0
    structure  = np.zeros([Natoms,3])
    AtomNumber = [6 for i in range(Natoms)]
    Atoms      = [i for i in range(Natoms)]
    Atomweight = np.zeros([Natoms])
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 1:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "Standard" and lineinfo[1] == "orientation:":
            flagnum = 1
        if lineinfo[0] == "Rotational" :
            flagnum = 0
        try:
            if int(lineinfo[0]) > 0 and flagnum == 1:
                AtomNumber[int(lineinfo[0])-1]       = int(lineinfo[1])
                for i in range(3):
                    structure[int(lineinfo[0])-1][i] = float(lineinfo[i+3])
        except:
            pass
        if lineinfo[1] == "Thermochemistry":
            thermoflag = 1
        if lineinfo[0] == "Molecular":
            thermoflag = 0
        try:
            if lineinfo[0] == "Atom" and thermoflag == 1:
                Atomweight[int(lineinfo[1])-1] = float(lineinfo[8])
        except:
            pass
        linenum += 1
    for atm in range(Natoms):
        Atoms[atm] = number_atom[str(AtomNumber[atm])]
    linecache.clearcache()
    return structure,Atoms,Atomweight

def chargespinreader(logfilename):
    linenum     = sum(1 for line in open(logfilename))
    charge_spin = [0,0]
    # Here we get charge and spin values
    for line in range(linenum):
        Linedata = linecache.getline(logfilename,line).split()
        try:
            if Linedata[0] == "Charge" and Linedata[3] == "Multiplicity":
                charge_spin = [int(Linedata[2]),int(Linedata[5])]
                break
        except:
            pass

    return charge_spin

def inputstructurereader(filename):
    # this function reads standard orientation of the molecule.
    # the last structure in the logfile is read
    linenumtot  = sum(1 for line in open(filename))
    Natoms      = atomnum(filename)
    flagnum = 0
    linenum = 1
    structure  = np.zeros([Natoms,3])
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 1:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "Input" and lineinfo[1] == "orientation:":
            flagnum = 1
        if lineinfo[0] == "Distance" or lineinfo[0] == "Standard":
            flagnum = 0
        try:
            if int(lineinfo[0]) > 0 and flagnum == 1:
                for i in range(3):
                    structure[int(lineinfo[0])-1][i] = float(lineinfo[i+3])
        except:
            pass
        linenum += 1

    linecache.clearcache()
    return structure

def inputfilereader(filename):
    # This function reads gaussian input file
    linenumtot  = sum(1 for line in open(filename))
    structure  = []
    atoms = []
    flag = 0
    linenum = 1
    chargespin = [0,0]
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) != 0:
            if flag == 2 and len(lineinfo) == 4:
                atoms.append(lineinfo[0])
                structure.append([float(lineinfo[1]),float(lineinfo[2]),float(lineinfo[3])])
            elif flag == 2 and len(lineinfo) == 2:
                chargespin = [int(lineinfo[0]),int(lineinfo[1])]
                
        else:
            flag += 1

        linenum += 1

    structure = np.array(structure)
    return structure,atoms,chargespin

def trajstructurereader(filename,filename2,number):
    Natoms      = atomnum(filename2)
    structure  = np.zeros([Natoms,3])
    linenumtot  = sum(1 for line in open(filename))
    for i in range(Natoms):
        lineinfo = linecache.getline(filename,linenumtot-number*(Natoms+2)+3+i).split()
        for j in range(3):
            structure[i][j] = float(lineinfo[j+1])
    linecache.clearcache()
    return structure

def getinputstructures(filename,strnum):
    Natoms      = atomnum(filename)
    inputstructures = []
    linenum = 0
    for i in range(strnum):
        inputstr,linenum = inputstructurereader_multiple(filename,linenum)
        inputstructures.append(inputstr)
    return inputstructures

def inputstructurereader_multiple(filename,startline):
    # this function reads standard orientation of the molecule.
    # the last structure in the logfile is read
    linenumtot  = sum(1 for line in open(filename))
    Natoms      = atomnum(filename)
    flagnum = 0
    linenum = startline
    structure  = np.zeros([Natoms,3])
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if flagnum > 0 and len(lineinfo) == 1 and "----" in lineinfo[0]:
           flagnum += 1 
        if len(lineinfo) <= 1:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "Input" and lineinfo[1] == "orientation:":
#            print("input orientation found, linenum = ",linenum)
            flagnum = 1
#        if flagnum == 1:
#            print(lineinfo)
#        if flagnum > 0 and ("------" in lineinfo[0]):
#            print("----- found, flagnum = ",flagnum)
#            flagnum += 1
        try:
            if int(lineinfo[0]) > 0 and flagnum > 1:
                for i in range(3):
                    structure[int(lineinfo[0])-1][i] = float(lineinfo[i+3])
        except:
            pass
        if flagnum == 4:
            break
        linenum += 1

    linecache.clearcache()
    return structure,linenum


def mopacinputstructurereader(filename):
    # this function reads standard orientation of the molecule.
    # the last structure in the logfile is read
    linenumtot  = sum(1 for line in open(filename))
    Natoms      = atomnum_mopac(filename)
    flagnum = 0
    linenum = 1
    structure  = np.zeros([Natoms,3])
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 1:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "ATOM" and lineinfo[1] == "CHEMICAL":
            flagnum = 1
        if lineinfo[0] == "CARTESIAN" or lineinfo[0] == "COORDINATES":
            flagnum = 0
            break
        try:
            if int(lineinfo[0]) > 0 and flagnum == 1:
                for i in range(3):
                    structure[int(lineinfo[0])-1][i] = float(lineinfo[2*i+2])
        except:
            pass
        linenum += 1
    linecache.clearcache()

    return structure

def grrmstructurereader(filename):
    # this function reads standard orientation of the molecule.
    # the last structure in the logfile is read
    linenumtot  = sum(1 for line in open(filename))
    Natoms      = atomnum_grrm(filename)
    flagnum = 0
    linenum = 1
    atoms = ["null" for i in range(Natoms)]
    structure  = np.zeros([Natoms,3])
    atom_number = 0
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 1:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "empty" and flagnum == 1:
            break
        if len(lineinfo) == 4 and flagnum == 1:
            atoms[atom_number] = lineinfo[0]
            for i in range(3):
                structure[atom_number][i] = float(lineinfo[i+1])
            atom_number += 1
        if lineinfo[0] == "Geometry" and lineinfo[1] == "(Origin":
            flagnum = 1
        linenum += 1
    linecache.clearcache()

    return structure,atoms


def forcereader(filename):
    Natoms      = atomnum(filename)
    linenumtot  = sum(1 for line in open(filename))
    forceflag = 0
    potentialE = 0
    linenum = 1
    force  = np.zeros([Natoms,3])
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 3:
            lineinfo = ["empty" for i in range(10)]
        if lineinfo[0] == "Energy=" and lineinfo[2] == "NIter=":
            potentialE = float(lineinfo[1])
        if lineinfo[0] == "SCF"  and lineinfo[1] == "Done:":
            potentialE = float(lineinfo[4])
        if lineinfo[0] == "ONIOM:"  and lineinfo[1] == "extrapolated":
            potentialE = float(lineinfo[4])
        if lineinfo[0] == "E2"  and lineinfo[3] == "EUMP2":
            potentialE = float(lineinfo[5][:-4])*1000
        # tempstring is not implemented
        if lineinfo[0] == "Total":
            try:
                if float(lineinfo[4]) < 0:
                    potentialE = float(lineinfo[4])
            except:
                pass
                
        if lineinfo[0] == "Center" and lineinfo[1] == "Atomic" and lineinfo[2] == "Forces":
            forceflag = 1
        try:
            if int(lineinfo[0]) > 0 and forceflag == 1:
                for i in range(3):
                    force[int(lineinfo[0])-1][i] = float(lineinfo[i+2])
        except:
            pass

        if lineinfo[1] == "Forces:":
            forceflag = 0

        linenum += 1
    linecache.clearcache()
    return force,potentialE 

def mopacforcereader(filename):
    Natoms = atomnum_mopac(filename)
    linenumtot  = sum(1 for line in open(filename))
    forceflag = 0
    potentialE = 0
    force  = np.zeros([Natoms,3])
    linenum = 1
    cart = {"X":0,"Y":1,"Z":2}
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) <= 3:
            lineinfo = ["empty" for i in range(10)]
        try:
            if lineinfo[7] == "KCAL/ANGSTROM":
                force[int(lineinfo[1])-1][cart[lineinfo[4]]] = - float(lineinfo[6]) * 0.52917725 / 627.509
        except:
            pass
        if lineinfo[0] == "FINAL" and lineinfo[1] == "HEAT":
            potentialE = float(lineinfo[5]) / 627.509
        linenum += 1

    return force,potentialE

def atomnum(filename):
    Natoms = 0
    linenum = 1 
    linenumtot  = sum(1 for line in open(filename))
    while linenum < linenumtot:
        lineinfo = linecache.getline(filename,linenum).split()
        if len(lineinfo) == 0:
            lineinfo = "empty"
        if lineinfo[0] == "NAtoms=":
            Natoms = int(lineinfo[1])
            break
        linenum += 1

    linenum = 1
    if Natoms == 0:
        flag = 0
        while linenum < linenumtot:
            lineinfo = linecache.getline(filename,linenum).split()
            if len(lineinfo) == 0:
                lineinfo = "empty"
            if lineinfo[0] == "Center":
                flag = 1
            try:
                if int(lineinfo[0]) > Natoms and flag == 1:
                    Natoms = int(lineinfo[0])
            except:
                pass
            try:
                if int(lineinfo[0]) < Natoms and flag == 1:
                    break
            except:
                pass
            linenum += 1
    #print(Natoms)        
    return Natoms


def atomnum_mopac(filename):
    Natoms = 0
    linenum = 1
    linenumtot  = sum(1 for line in open(filename))
    if Natoms == 0:
        flag = 0
        while linenum < linenumtot:
            lineinfo = linecache.getline(filename,linenum).split()
            if len(lineinfo) == 0:
                lineinfo = "empty"
            if lineinfo[0] == "ATOM":
                flag = 1
            try:
                if int(lineinfo[0]) > Natoms and flag == 1:
                    Natoms = int(lineinfo[0])
            except:
                pass
            try:
                if int(lineinfo[0]) < Natoms and flag == 1:
                    break
            except:
                pass
            linenum += 1
    #print(Natoms)        
    return Natoms


def logreader(logfilename):
    charge_spin,electronic_energy,atoms,coordinates,\
    opt_status,Conversion = outputreader(logfilename)
    if np.sum(coordinates) == 0:
        charge_spin,electronic_energy,atoms,coordinates,\
        opt_status,Conversion = outputreader(logfilename,1)
    return charge_spin,electronic_energy,atoms,coordinates,opt_status,Conversion

def outputreader(logfilename,option=0):
    linenum     = sum(1 for line in open(logfilename))
    charge_spin = [0,0]
    electronic_energy = 0.0
    opt_status = 0
    # Here we get electronic energy
    for line in range(linenum):
        Linedata = linecache.getline(logfilename,line).split()
        try:
            if Linedata[0] == "SCF" and Linedata[1] == "Done:":
                electronic_energy = float(Linedata[4])
        except:
            pass

    # Here we get charge and spin values
    for line in range(linenum):
        Linedata = linecache.getline(logfilename,line).split()
        try:
            if Linedata[0] == "Charge" and Linedata[3] == "Multiplicity":
                charge_spin = [int(Linedata[2]),int(Linedata[5])]
                break
        except:
            pass

    # Get Conversion status
    # Conversion = [Maximum_Force, RMS_Force,Maximum_Displacement, RMS_Displacement]
    # in each value of list, [Value,Threshold,Converged?]
    Conversion = [[0.0,0.0,"Null"],[0.0,0.0,"Null"],[0.0,0.0,"Null"],[0.0,0.0,"Null"]]
    for line in range(linenum):
        Linedata = linecache.getline(logfilename,line).split()
        [Linedata.append("Null") for i in range(5)]
        if Linedata[0] == "Maximum" and Linedata[1] == "Force":
            for i in range(3):
                Conversion[0][i] = Linedata[i+2]
        if Linedata[0] == "RMS" and Linedata[1] == "Force":
            for i in range(3):
                Conversion[1][i] = Linedata[i+2]
        if Linedata[0] == "Maximum" and Linedata[1] == "Displacement":
            for i in range(3):
                Conversion[2][i] = Linedata[i+2]
        if Linedata[0] == "RMS" and Linedata[1] == "Displacement":
            for i in range(3):
                Conversion[3][i] = Linedata[i+2]
        if Linedata[1] == "Optimized" and Linedata[2] == "Parameters":
            break

    # Here we check optimization status and get atoms and coordinates
    optflag = 0
    stdflag = 0
    atoms   = []
    coordinates = []
    for line in range(linenum):
        Linedata = linecache.getline(logfilename,line).split()
        if len(Linedata) <= 3:
            [Linedata.append("null") for i  in range(10)]
        if Linedata[1] == "Optimized" and Linedata[2] == "Parameters":  #and optflag == 0:
            opt_status = "Optimized"
            coordinates = []
            optflag = 1
        if Linedata[1] == "Non-Optimized" and Linedata[2] == "Parameters": #and optflag == 0:
            opt_status = "Non-Optimized"
            optflag = 1
            coordinates = []
        if option == 0: 
            if optflag == 1 and Linedata[0] == "Standard" and Linedata[1] ==  "orientation:":
                stdflag = 1
        elif option == 1:
            if optflag == 1 and Linedata[0] == "Input" and Linedata[1] ==  "orientation:":
                stdflag = 1
        try:
            if optflag == 1 and stdflag == 1 and int(Linedata[0]) > 0:
                atoms.append(number_atom[Linedata[1]])
                coordinates.append([float(Linedata[3]),float(Linedata[4]),float(Linedata[5])])
        except:
            pass
        if (Linedata[0] == "Rotational" and stdflag == 1) or \
            (Linedata[0] == "Distance" and stdflag == 1):
            optflag = 0
            stdflag = 0
            break

    coordinates = np.array(coordinates)
    linecache.clearcache()

    return charge_spin,electronic_energy,atoms,coordinates,opt_status,Conversion

def IRCreader(filename):
    linenum  = sum(1 for line in open(filename))
    strflag = 0
    Pt = 0
    electronic_energy = 0
    atomarray = []
    coordinates = []
    charge_spin = [0,0]
    for line in range(linenum):
        Linedata = linecache.getline(filename,line).split()
        [Linedata.append("Null") for i in range(10)]
        if Linedata[0] == "Charge" and Linedata[3] == "Multiplicity":
            charge_spin[0] = int(Linedata[2])
            charge_spin[1] = int(Linedata[5])
        if Linedata[0] == "Input" and Linedata[1] ==  "orientation:":
            strflag = 1
            atomarray = []
            coordinates = []
        if Linedata[0] == "Pt":
            Pt = Linedata[1]
        if "---" in Linedata[0] and strflag == 2:
            strflag = 0
        if strflag >= 1 and Linedata[0].isdigit() == True:
            atomarray.append(number_atom[Linedata[1]])
            coordinates.append([float(Linedata[3]),float(Linedata[4]),float(Linedata[5])])
            strflag = 2
        if line == linenum-1:
            if Pt == 0:
                print("\n\n\n           IRC   FAILED             \n\n\n")
                sys.exit()

    return charge_spin,electronic_energy,atomarray,coordinates,Pt

def Num_of_imag_freq(filename):
    linenum  = sum(1 for line in open(filename))
    num_of_imag_freq = 100
    for line in range(linenum):
        Linedata = linecache.getline(filename,line).split()
        [Linedata.append("null") for i  in range(3)]
        if Linedata[0] == "Frequencies" and Linedata[1] == "--":
            Firstfreq  = Linedata[2]
            Secondfreq = Linedata[3]
            print("")
            if float(Firstfreq) > 0:
                num_of_imag_freq = 0 
            elif float(Firstfreq) < 0 and float(Secondfreq) > 0:
                num_of_imag_freq = 1
            else:
                num_of_imag_freq = 2
            break
    return num_of_imag_freq

def thermalinfo(filename):
    linenum  = sum(1 for line in open(filename))
    ZPE = 0
    UCor = 0
    HCor = 0
    GCor = 0
    EZPE = 0
    UE = 0
    HE = 0
    GE = 0
    for line in range(linenum):
        RawLinedata = linecache.getline(filename,line)
        Linedata = RawLinedata.split()
        [Linedata.append("null") for i  in range(8)]
        if "Zero-point correction=" in RawLinedata:
            ZPE = float(Linedata[2])
        elif "Thermal correction to Energy=" in RawLinedata:
            UCor = float(Linedata[4])
        elif "Thermal correction to Enthalpy=" in RawLinedata:
            HCor = float(Linedata[4])
        elif "Thermal correction to Gibbs Free Energy=" in RawLinedata:
            GCor = float(Linedata[6])
        elif "zero-point Energies=" in RawLinedata:
            EZPE = float(Linedata[6])
        elif "thermal Energies=" in RawLinedata:
            UE = float(Linedata[6])
        elif "thermal Enthalpies=" in RawLinedata:
            HE = float(Linedata[6])
        elif "thermal Free Energies=" in RawLinedata:
            GE = float(Linedata[7])
    
    return ZPE,UCor,HCor,GCor,EZPE,UE,HE,GE

def rootsectionread(filename):
    linenum  = sum(1 for line in open(filename))
    rootinfo = ""
    inputreadflag = 0
    for line in range(linenum):
        RawLinedata = linecache.getline(filename,line).replace('\n', '')
        Linedata = RawLinedata.split()
        [Linedata.append("null") for i  in range(8)]
        if Linedata[0] == "#" or Linedata[0] == "#p" and inputreadflag == 0:
            rootinfo = RawLinedata
            addnum = 1
            while True:
                RawLinedata = linecache.getline(filename,line+addnum).replace('\n', '')
                if "---" in RawLinedata:
                    inputreadflag = 1
                    break
                rootinfo += RawLinedata
                addnum += 1
        if inputreadflag == 1:
            break
    rootinfo = rootinfo.lower()
    return rootinfo

def printgeo(atomarray,coordinates):
    for N in range(len(atomarray)):
        print(" {:<2} ".format(atomarray[N])+" {: >9.6f}".format(coordinates[N][0])
              +" {: >9.6f}".format(coordinates[N][1])+" {: >9.6f}".format(coordinates[N][2]))
    print("")


def RootProcessing(option_interpreted,rootinfo_separated,option="opt"):
    rootlength = 0
    print(" ",end="")
    for i in range(len(rootinfo_separated)):
        if option in rootinfo_separated[i]:
            print(option_interpreted,end=" ")
            rootlength += len(option_interpreted)
        elif "freq" in rootinfo_separated[i]:
            pass 
        elif ")" in rootinfo_separated[i] and "(" in rootinfo_separated[i]:
            print(rootinfo_separated[i],end=" ")
            rootlength += len(rootinfo_separated[i])
        elif ")" not in rootinfo_separated[i] and "(" in rootinfo_separated[i]:
            rootinfo_separated[i] += rootinfo_separated[i+1]
            print(rootinfo_separated[i],end=" ")
            rootlength += len(rootinfo_separated[i])
        elif ")" in rootinfo_separated[i] and "(" not in rootinfo_separated[i]:
            pass
        else:
            print(rootinfo_separated[i],end=" ")
            rootlength += len(rootinfo_separated[i])
        if rootlength > 50 and i < len(rootinfo_separated)-1:
            print("")
            rootlength = 0
    

def IRCoptimizemake(filename,option):
    charge_spin,electronic_energy,atomarray,coordinates,Pt = IRCreader(filename)
    rootinfo = rootsectionread(filename)
    if "irc" not in rootinfo.lower():
        print("its not IRC calculation")
        sys.exit()
    rootinfo_separated = rootinfo.split()
    RootProcessing("opt=(tight,calcfc) freq=hpmodes",rootinfo_separated,"irc")
    print("\n\n "+filename+" "+option+" Pt= "+str(Pt)+"\n")
    print(" "+str(charge_spin[0])+" "+str(charge_spin[1]))
    printgeo(atomarray,coordinates)


if __name__ == "__main__":
    filename = sys.argv[1]
    inputstructures = getinputstructures(filename,15)
    print(inputstructures)


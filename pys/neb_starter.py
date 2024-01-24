# neb_starter.py
# ver 0.6 2020/10/26
#written by Hiroaki Kurouchi

# Standard library
import os
import sys
import hashlib
import datetime

# Original modules
import neb_initialization
import neb_mainprocess
import neb_optimization

#ここで重要なパスをファイルとして書き込む
def prepare_for_calculation(origdir,programdir,scratchdir,inputfilename):
    # インプットファイルの中身のハッシュ値が同じであればリスタートする。
    # 書き換わっていた場合はLogfilesをLogfiles_oldにリネームして新規スタート。
    with open(origdir+inputfilename,"r") as f:
        inputdata = f.read()

    hash_input = hashlib.sha224(inputdata.encode()).hexdigest()
    time_now = datetime.datetime.now()

    if os.path.exists(origdir+"Logfiles") == False:
        print("No Logfile directory. Start new calculation.")
        os.mkdir(origdir+"Logfiles")
        with open(origdir+"Logfiles/hash_inputfile","w") as w:
            w.write(hash_input)
        with open(origdir+"Logfiles/scratchdir","w") as w:
            w.write(scratchdir)
        with open(origdir+"Logfiles/programdir","w") as w:
            w.write(programdir)
        with open(origdir+"Logfiles/inputfilename","w") as w:
            w.write(inputfilename)
        


    elif os.path.exists(origdir+"Logfiles/hash_inputfile") == True:
        with open(origdir+"Logfiles/hash_inputfile","r") as r:
            if hash_input == r.read():
                print("The same input file. Continue calculation.")

            else:
                print("Different input file. Start new calculation.")
                os.rename(origdir+"Logfiles",origdir+"Logfiles_"+time_now.strftime("%Y%m%d%H%M"))
                os.mkdir(origdir+"Logfiles")
                with open(origdir+"Logfiles/hash_inputfile","w") as w:
                    w.write(hash_input)
                with open(origdir+"Logfiles/scratchdir","w") as w:
                    w.write(scratchdir)
                with open(origdir+"Logfiles/programdir","w") as w:
                    w.write(programdir)
                with open(origdir+"Logfiles/inputfilename","w") as w:
                    w.write(inputfilename)


def nebprocess():
    ### 初期化プロセス 
    neb_initialization.initialize()

    ###  Force計算と構造最適化のループ 
    neb_mainprocess.nebloop()

    ### 極大値、初期構造の最適化を試みる。QST3を利用 
    ### その後、極大値からIRC計算・最適化計算を走らせる 
    neb_optimization.optimization_process()
    
#    neb_functions.nebsummarize(bandarray_stationary_opt, bandarray_stationary_en)

    sys.exit()

if __name__ == '__main__':
    # This program must be run in the origdir.
    # Run this program as "python3 (path)/neb_starter.py  (inputfilename) (scratchdirpath)"

    origdir = os.getcwd()
    if origdir[-1] != "/":
        origdir += "/"

    args = sys.argv
#    programdir = str(args[1])
    inputfilename = str(args[1])
    scratchdir = str(args[2])

    # Get programdir from args[0]
    programdir = args[0].replace("neb_starter.py","")
    if programdir == "":
        programdir = origdir
    else:
        os.chdir(programdir)
        programdir = os.getcwd()
        os.chdir(origdir)
        if programdir[-1] != "/":
            programdir += "/"

    if ".." in scratchdir:
        os.chdir(scratchdir)
        scratchdir = os.getcwd()
        os.chdir(origdir)
        if scratchdir[-1] != "/":
            scratchdir += "/"

    if os.path.exists(origdir+inputfilename) == False:
        print("specify correct inputfilename")
        sys.exit()

    if os.path.exists(scratchdir) == False:
        print("specify scratchdir using absolute path")
        sys.exit()

    if scratchdir[-1] != "/":
        scratchdir = scratchdir + "/"

    print("programdir:{} origdir=origdir:{} scratchdir:{}".format(programdir,origdir,scratchdir))

    prepare_for_calculation(origdir,programdir,scratchdir,inputfilename)
    nebprocess()
#    os.removedirs(scratchdir)
    sys.exit()



import numpy as np
import os
from multiprocessing import Pool
import sys

class SourceReceiverPair:
    def __init__(self, elat = 0.0, elon = 0.0, edep = 0.0, phase="P"):
        self.elat = elat
        self.elon = elon
        self.edep = edep
        self.phase=phase

    def get_stations(self,stainfo,with_elevation=True):
        n = len(stainfo)
        self.num_station = n
        self.cords= np.zeros((n,3),dtype=float)

        for i in range(n):
            elev,lat,lon = list(map(lambda x:float(x),stainfo[i].split()[:3]))
            self.cords[i,0] = lat
            self.cords[i,1] = lon
            if not with_elevation:
                elev = 0.0
            self.cords[i,2] = elev 

    def __str__(self):
        out = "The event is at {},{},{} ".format(str(self.elat),str(self.elon),str(self.edep))  + "\n"
        out += "There are {} station in this event \n ".format(self.num_station)
        for i in range(self.num_station):
            out += " {},{} \n".format(self.cords[i,0],self.cords[i,1])
        return out

    
def initialize(srcfile,rcvfile,with_elevation=True):
    # initialize source-receivers pairs
    Pairs = []

    # read receiver files
    with open(rcvfile,"r") as f:
        stainfo = f.readlines()
    
    # compute receiver index for each array
    na = int(stainfo[0].split()[0]) # no. of arrays
    stas_each_arr = []
    idx = 1
    for _ in range(na):
        nsta = int(stainfo[idx].split()[0])
        #print(nsta)
        idx = idx + 1
        stas_each_arr.append(stainfo[idx:idx+nsta])
        idx = idx + nsta 
    
    # read sources 
    f = open(srcfile,"r")
    f.readline()
    for i in range(na):
        line = f.readline()
        nsrc = int(line)
        for _ in range(nsrc):
            line = f.readline()
            elat,elon,edep = list(map(lambda x:float(x),line.split()))
            line = f.readline()
            phase = line.split()[0]
            pair = SourceReceiverPair(elat,elon,edep,phase=phase)
            pair.get_stations(stas_each_arr[i],with_elevation=with_elevation)
            Pairs.append(pair)
    f.close()

    return Pairs
    
def prepare_files(nproc,source,receiv):
    # read station and source information
    Pairs = initialize(source,receiv,True)
    nsrc = len(Pairs)

    for pid in range(nproc):
        # compute index
        start = nsrc // nproc * pid
        end = nsrc // nproc * (pid+1)
        if pid == nproc - 1 :
            end= nsrc

        # mkdir
        os.system("mkdir -p fm3dt" +str(pid) )

        # copy files to new directorys
        os.system('cp -r ak135.hed ak135.tbl ak135.vel fm3dt'+str(pid))
        os.system('cp *.vtx itimes.in aktsurf.in fm3dt.in fm3dt' + str(pid))

        # write sources and receivers to this new directory
        os.chdir("fm3dt"+str(pid))
        na = end - start # no. of arry
        f = open(source,"w")
        g = open(receiv,"w")
        f.write("%d\n"%(na))
        g.write("%d\n"%(na))
        for i in range(start,end):
            p = Pairs[i]
            f.write("1\n")
            f.write("%f %f %f\n"%(p.elat,p.elon,p.edep))
            f.write(p.phase+'\n')
            g.write("%d\n"%(p.num_station))
            for j in range(p.num_station):
                g.write("%f %f %f\n"%(p.cords[j,2],p.cords[j,0],p.cords[j,1]))

        f.close()
        g.close()

        os.chdir('..')

def run_fm3dt(args):
    istep = args[-1]
    cmd = args[-2]
    pid = args[0]

    # cd to corresponding directory
    filename = 'fm3dt'+str(pid)
    os.system("cp gridc.vtx "+filename)
    os.chdir(filename)

    # if for the first time, prepare all files
    cmd_path = '../'+cmd.split('fm3dt')[0]
    if istep == 1:
        os.system(cmd_path+'aktsurf')
        os.system(cmd_path+'itimes')
    
    # run fm3dt
    print("computing traveltime and frechet kernel, process: ",pid+1)
    os.system('../'+cmd)
    os.chdir('..')

def merge(nproc,filename):
    cmd = 'cat '
    for i in range(nproc):
        cmd += 'fm3dt'+str(i)+'/'+filename + ' '
    cmd += '> ' + filename
    os.system(cmd)

def main(argv):
    if len(argv) == 1:
        cmd  = "../bin/fm3dt"
        istep = 1
        nproc = 4
    else:
        cmd = argv[1]
        istep = int(argv[2])
        nproc= int(argv[3])

    # read source and receiver information
    srcfile = "sources.dat"
    rcvfile = "receivers.dat"
    sfile = open(srcfile,"r")
    line = sfile.readline()
    sfile.close()
    nsrc = int(line)

    # prepare works for each process
    args = []
    for i in range(nproc):
        start = nsrc // nproc * i 
        end = nsrc // nproc * (i+1)
        if i== nproc -1 :
            end = nsrc
        arg = (i,start,end,srcfile,rcvfile,cmd,istep)
        args.append(arg)
    
    # if iteration for the first time, synthetic times to bottom
    if istep == 1:
        print("computing traveltimes to bottom of research area ...")
        prepare_files(nproc,srcfile,rcvfile)

    # prepare parallel computation
    pool = Pool(processes=nproc)
    pool.map(run_fm3dt,args)
    pool.close()
    pool.join()

    # merge traveltime and frechet kernel
    merge(nproc,"frechet.out")
    merge(nproc,"rtravel.out")

if __name__ == "__main__":
    main(sys.argv)


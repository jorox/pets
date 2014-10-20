import numpy as np
import sys

def main():
    """
    Takes in a series of dump files and produces the average layer pressure
    input parameters: <scriptname> filename coln(1-based) natoms minx:dx:maxx outname
    """
    if len(sys.argv) < 5:
        print "!! error insufficient arguments !!"
        print "Needed: filename col nsteps natoms"
        return 1
    
    nargs = len(sys.argv)
    nfiles = nargs-4-1          #-1 = scriptname, -3 = extra arguments
    
    fls = sys.argv[1:-4]     # file names  
    col = int(sys.argv[-4])      # column marking pressure data (or any other for that matter!)
    natoms = int(sys.argv[-3])   # number of atoms
    rng = sys.argv[-2]           # range of x
    outname = sys.argv[-1]       # output file name
    
    rng = rng.split(":")
    nsteps = nfiles              #number of steps
    minx = float(rng[0])
    dx = float(rng[1])     #layer thickness
    maxx = float(rng[2])
    
    print "\n... files = " + str(nfiles)
    print "... column   = " + str(col)
    print "... # of steps = " + str(nsteps)
    print "... " + str(natoms) + " atoms \n"
    print "... working :"
    nbins = int(np.ceil((maxx-minx)/dx)) #number of bins
    freq = []; avpr = []
    for i in range(nsteps):
        freq.append([0]*nbins)        #frequency for each bin
        avpr.append([0]*nbins)        #average pressure for each bin
    binmins = np.arange(0,nbins)*dx+minx   #lower edge of each bin. Note that: len(binmins) = nbins-1
    freqout = open(outname+".freq", "w")
    prssout = open(outname+".prss", "w")
    
##### loop over all files
    itstep = 0 #step index
    for flname in fls:
        print "    " + flname + "... ",
        sys.stdout.flush()
        fl = open(flname)
        
        cntr = 0 #atom counter
        atmsec = False # signals beginning of atoms section
        
        atoms_x = [0.]*natoms
        atoms_P = [0.]*natoms
        per = ""

    ##### read atom data from file
        print " reading ...",
        while 1:
            if cntr == natoms:
                break
            line = fl.readline()
            #print line
            if atmsec:
                line = line.strip().split()
                id = int(line[0])-1 #first item is always the id
                x = float(line[col-2]) # coordinate always preceeds the pressure column
                P = float(line[col-1]) # col is added as base-1 numbering
                
                atoms_x[id] = x
                atoms_P[id] = P
                
                cntr += 1
            
            if line == "ITEM: TIMESTEP\n":
                step = fl.readline()
                fl.readline()  #ITEM: Number of Atoms
                if natoms != int(fl.readline()):
                    print "error number of atoms not correct"
                    return 1
            
            if line[:11] == "ITEM: ATOMS":
                atmsec = True
                continue
        
    ##### bin the atoms for this timestep/file
        print "binning... " + str(itstep),
        sys.stdout.flush()
        for i in range(natoms):
            binned = False
            
            if float(i)%50 == 0:
                statusdone = "%06.2f"%(float(i)/natoms*100) + "%"
                print statusdone,
            for j in range(1,nbins):
                wall = binmins[j]
                if atoms_x[i]<= wall:
                    binned = True
                    freq[itstep][j-1] +=1
                    avpr[itstep][j-1] +=atoms_P[i]
                    break
            if float(i)%50 == 0:
                print "\b"*9,
            if binned: continue
            freq[itstep][-1] += 1
            avpr[itstep][-1] += atoms_P[i]
        #--- end of atom bin loop
        itstep +=1
        print ""
        
##### write results to file
    print "writing ..."
    ### t1 t2 t3 ..... tn
    ### b1
    ### b2
    ### b3
    ### ...
    ### bn
    for ibin in range(nbins):
        
        cntr = (binmins[ibin]+dx)
        freqout.write("\n%1.4f" % (cntr))
        prssout.write("\n%1.4f" % (cntr))
        
        for it in range(nsteps):
            itbfreq = freq[it][ibin]
            itbprss = avpr[it][ibin]
             
            if itbfreq < 1: # found an empty bin ... 1
                if itbprss != 0.0: # check that there is no pressure values ...1a
                    print "error incompatible freq and pressure zero values"
                    return 1
                    # set to 1 to avoid division by zero
                #print " warning/" + str(ibin) + " " + str(it),
                itbfreq=1
                
            freqout.write(" %1.0f"%(freq[it][ibin])) # we use this value to print because itbfreq is adjusted to avoid zero division
            prssout.write(" %1.5e"%(itbprss/float(itbfreq)))
            
    freqout.close()
    prssout.close()
    print "... exited successfully :)"

if __name__ == "__main__":
    main()
    

import numpy as np
import sys
import matplotlib.pyplot as plt

def main():
    """
    input parameters: <scriptname> filename coln(1-based) nsteps natoms
    """
    if len(sys.argv) < 5:
        print "!! error insufficient arguments !!"
        print "Needed: filename col nsteps natoms"
        return 1
    
    flname = sys.argv[1]
    col = int(sys.argv[2])
    nsteps = int(sys.argv[3])
    natoms = int(sys.argv[4])
    
    print "... filename = " + flname
    print "... column   = " + str(col)
    print "... # of steps = " + str(nsteps)
    print "... " + str(natoms) + " atoms"
    
    fl = open(flname)
    
    cntr1 = 0  #step counter
    cntr2 = 0 #atom counter
    atmsec = False # signals beginning of atoms section
    
    atoms_x = [0]*natoms
    atoms_P = [0]*natoms
    per = ""
    print "... working >>>>"
    
    while 1:
        
        if cntr1 == nsteps:
            # gathered all steps
            temp = atoms_x 
            atoms_x = [x/cntr1 for x in temp]
            temp = atoms_P
            atoms_P = [x/cntr1 for x in temp]
            print "... done"
            break
        
        line = fl.readline()
        #print line
        if atmsec:
            line = line.strip().split()
            id = int(line[0])-1 #first item is always the id
            x = float(line[1]) #second item is always the coordinate
            P = float(line[col-1])
            
            atoms_x[id] += x
            atoms_P[id] += P
            
            cntr2 += 1
            if cntr2 == natoms:
                # Done adding atoms
                atmsec = False
        
        if line == "ITEM: TIMESTEP\n":
            # found a new timestep
            cntr1 += 1
            step = fl.readline()
            print "\b"*len(per),            
            per = "%3.2f"%(float(cntr1)/nsteps*100)+"%"            
            print per,
            fl.readline()  #ITEM: Number of Atoms
            if natoms != int(fl.readline()):
                print "error"
                return 1
        
        if line[:11] == "ITEM: ATOMS":
            atmsec = True
            cntr2 = 0
            continue
        
    plt.plot(atoms_x,atoms_P,'r.')
    plt.imshow
    
if __name__ == "__main__":
    main()
    
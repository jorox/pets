import sys
import numpy as np
def main():
    # Parse command line arguments
    #   arg[1] = input file name
    #   arg[2] = output file name
    #   -d D  = diameter of the wire
    #   -l L  = length of nanowire
    #   -q T  = atoms in file have a charge
    #   <-edge> in,out = True, False
    #   <-ortho> y,n    = Change to orthogonal struct?
    #   <-lammps> y,n   = outpuf file format
    #   <-map> y,n      = map file?
    #   <-rough> n l r   = create n modulations l Angstroms wide of ratio r 
    
    n = [1]*3
    flname = sys.argv[1]
    print "... number of args = " + str(len(sys.argv))
    outname = sys.argv[2]
    nw_diam = 0.0 # diameter of the nanowire
    nw_length = 0.0
    doedge = True
    dorough = False
    ## extra options imported from vec2lammpsbox needed to make the map file
    domap = False
    doalloy = False
    lammps = True
    ischarge = False  #******making a map file is not yet supported***********#
    
### Parse other input parameters
    args = sys.argv[3:]
    nargs = len(args)
    if nargs>0:
        for iarg in range(nargs):
            if args[iarg] == "-d":
                nw_diam = float(args[iarg+1])
                continue
            if args[iarg] == "-edge":
                if args[iarg+1] == "no" or args[iarg+1] == "n" or args[iarg+1] == "N" or args[iarg+1] == "0":
                    doedge = False
                    continue
            if args[iarg] == "-l":
                nw_length = float(args[iarg+1])
                continue
            if args[iarg] == "-lammps":
                if args[iarg+1] == "n" or args[iarg+1] == "N":
                    lammps = False
            if args[iarg] == "-rough":
                dorough = True
                mod_n = int(args[iarg+1])
                mod_l = float(args[iarg+2])
                mod_ratio = int(args[iarg+3])          
            if args[iarg] == "-q":
                ischarge = True
                continue
                
            #if args[iarg] == "-map":
            #    if args[iarg+1] == "y" or args[iarg+1] == "Y":
            #        domap = True
    
#### load text from input file, ignore lines with 
    ABC = np.loadtxt(flname)
    # unit cell vectors will be the first three lines 
    #   take only the first three entries 
    #   (trailing zeros for book-keeping)
    A = ABC[0][:3]
    B = ABC[1][:3]
    C = ABC[2][:3]
    # build the lammps basis a,b,c using the relations 
    #   specified on the website
    a = [norm(A),0.0,0.0]
    Ahat = unit(A)
    b = [dot(B,Ahat), norm(np.cross(Ahat,B)), 0]
    ABhat = unit(np.cross(A,B))
    c = [dot(C,Ahat), dot(C,np.cross(ABhat,Ahat)), 0.]
    c[2] = np.sqrt(norm(C)**2-c[0]**2-c[1]**2)
    
#### Calculate the required number of unit cells
    
    # Calculate the nanowire cross-sectional area
    nw_area = np.pi*nw_diam**2/4.0
    # Calculate the area of the sectional face (for now along the c-axis)
    sq_area = norm(a)*norm(b)
    # Calculate the ratio        a*b*ratio = area_nw
    #   round to nearest integer ratio     = round(ratio)
    #   find the square root     sqrt_ratio= a_ratio*b_ratio
    ratio = np.ceil(nw_area/sq_area)
    ratio = int(np.ceil(np.sqrt(ratio)))
    #   the number of unit cells needed is ceiling(sqrt(ratio))
    #   the number of unit cells for the length is also given easily
    z_ratio = int(np.ceil(nw_length/norm(c)))
    n = [ratio,ratio,z_ratio]
    
#### If there are more lines these will be the atoms in the unit cell ****which should be the case****

    if len(ABC) > 3:
        nuc = len(ABC[3:])
        ntypes = 1  # at least one atom type
        ## find the number of atom types by taking the maximum
        for elem in ABC[3:]:
            if elem[3] > ntypes: ntypes = elem[3] # number of atom types
        print "... found atoms (" + str(nuc) + ") in text file, changing basis"
        douc = True
    # if found count number of atoms in unit cell 
    #    and transform. Build more if needed
    # ucatms = [ [x1,y1,z1,t1], [x2,y2,z2,t2], .... , [xn,yn,zn,tn]]
    if douc:
        ucatms = ABC[3:]
        atoms = transUC(ucatms,[A,B,C], [a,b,c])
        if n[0]>1 or n[1]>1 or n[2]>1:
            print "... need more atoms adding " + str(n) 
            atoms = buildmore(n,[a,b,c],atoms)
    print "     >>> " + str(len(atoms)) + " atoms"
    sysm = [[n[0]*a[0], 0 ,0], [n[1]*b[0], n[1]*b[1], 0], [n[2]*c[0], n[2]*c[1], n[2]*c[2]]]
    print "     >>> %1.4f, %1.4f, %1.4f nm"% (norm(sysm[0])/10, norm(sysm[1])/10, norm(sysm[2])/10)
    
### build nanowire by cropping atoms outside the circle diameter
    kill_list = []
    # First find midpoint of the square region
    mdpt_x = norm(a)*float(n[0])/2.0
    mdpt_y = norm(b)*float(n[1])/2.0
    nw_radius = nw_diam/2.0
    # loop over atoms
    for ia in range(len(atoms)):
        atm = atoms[ia]
        # calculate radius of atom
        atm_R = ((mdpt_x-atm[0])**2 + (mdpt_y-atm[1])**2)**0.5
        # if outside, remove it
        if atm_R > nw_radius:
            kill_list.append(ia)
    print "######## DEBUG ########"
    #print kill_list
    print len(atoms)
    print len(kill_list)
    print "#######################"
    for ia in range(len(kill_list)-1,-1,-1):
        try:
            atoms.pop(kill_list[ia])
        except:
            print "error",
           # print "problem with index %1.0f, atom %1.0f"%(ia, kill_list[ia])
    print "... removed %1.0f atoms for nanowire, %1.0f remaining"%(len(kill_list), len(atoms))
    
#### modulate the nanowire
   #       d
   # <---------->
   # {   }      {   }
   # <-->     
   #   l
   #                             _________L___________
   #                             {__}___{__}___{  }___ 
   #                             <-----><-----><----->
   #                                 d      d     d       ==> n = 3
   
    if dorough:
        mod_r2 = (nw_diam/2.0/mod_ratio)**2 #cutoff radius
        mod_d = norm(c)*n[2]/mod_n          #number of modulations
        print ".... roughening "
        print "     >>> period = %1.5f nm"%(mod_d/10.0)
        print "     >>> ratio = %1.0f (%1.4fnm)"%(mod_ratio, np.sqrt(mod_r2)/10.0)
        print "     >>> width = %1.5f nm"%(mod_l/10)
        if mod_d < mod_l:
            print "!!!error too many modulations!!!!"
            print "...... modulation length is %1.4f"%mod_d
            return 1
        kill_list=[]
        for nmod in range(mod_n):
            for ia in range(len(atoms)):
                atm = atoms[ia]
                atm_r2 = (mdpt_x-atm[0])**2 + (mdpt_y-atm[1])**2  # radius of atom
                if nmod*mod_d <= atm[2] < nmod*mod_d+mod_l:
                    if  atm_r2 > mod_r2: 
                        kill_list.append(ia)
        
        for ia in range(len(kill_list)-1,-1,-1):
            atoms.pop(kill_list[ia])
        
        print "... moded %1.0f atoms for nanowire, %1.0f remaining"%(len(kill_list), len(atoms))
    
#### write to file in either lammps format or in same format as input text
    fout = open(outname,'w')
    
    if lammps:
        fout.write("Cell with dimensions " + str(n))
        if doalloy: fout.write(" changed " + str(base) + " to " + str(alloy) +  
                               " with conc. " + str(float(len(changeList))/len(atoms)) + " at.")
        fout.write("\n\n")
        fout.write(str(len(atoms)) + " atoms\n%1.0f atom types\n\n"%(ntypes))
    
        fout.write("%10.8f %10.8f xlo xhi" % (0.0, n[0]*a[0]))
        fout.write("\n%10.8f %10.8f ylo yhi" % (0.0, n[1]*b[1]))
        fout.write("\n%10.8f %10.8f zlo zhi" % (0.0, n[2]*c[2]))

        if b[0] != 0.0 or c[0] != 0.0 or c[1] != 0.0:
            print "!! non-orthogonal box"
            fout.write("\n%10.8f %10.8f %10.8f xy xz yz\n\n" % (n[1]*b[0], n[2]*c[0], n[2]*c[1]))
    else:
        fout.write(str(len(atoms)) + "\n")
    if douc:
        ia = 0
        if lammps: fout.write("\n\nAtoms\n")
        else: fout.write("Atoms")
        for elem in atoms:
            ia +=1
            if lammps:
                if ischarge:
                    fout.write("\n%1.0f %1.0f %1.3f %10.8f %10.8f %10.8f" % (ia, elem[3], 1.011, elem[0], elem[1], elem[2]))
                else:
                    fout.write("\n%1.0f %1.0f %10.8f %10.8f %10.8f" % (ia, elem[3], elem[0], elem[1], elem[2]))
            else:
                fout.write("\n%1.0f %10.8f %10.8f %10.8f" % (elem[3], elem[0], elem[1], elem[2]))
                
    fout.close()
    print "... Done writing to file " + outname
#### Done writing to file
#### Write to map file
    if douc and domap:
        outname = "map." + outname
        fout = open(outname, 'w')
        
        fout.write("%1.0f %1.0f %1.0f %1.0f" % (n[0], n[1], n[2], nuc)) #header
        fout.write("\n#l1 l2 l3 k type atom_id")
        
        inv_abc = np.linalg.inv(np.transpose(np.array([a,b,c])))
        ia = 0
        for elem in atoms:
            ia +=1
            type = elem[3]
            elem = elem[:3]
            elem = np.dot(inv_abc, np.transpose(np.array(elem))) # transform to fractional coordinates
            elem = [int(x) for x in elem]                        # take the integer part of the coordinate
            
            fout.write("\n%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f" % 
                       (elem[0], elem[1], elem[2], (ia-1)%nuc, type, ia))
        print "... Done writing map file: " + outname
        fout.close()
        
def dot(a,b):
    res = 0
    for x,y in zip(a,b):
        res += x*y
    return res

def norm(a):
    res = 0;
    for x in a:
        res += x**2
    res = res**0.5
    return res

def unit(a):
    return a/norm(a)
  
def transUC(atms, oldbasis, newbasis):
    """
    change from old basis which is the input file basis to the new basis i.e. the lammps basis
    Does not check whether new atoms are needed
    """
    
    newbasis = np.transpose(np.array(newbasis)) #form the V matix

    from numpy.linalg import inv
    inv_newbasis = inv(newbasis) # invert
    
    CBM = []
    for u in oldbasis:
        u = np.dot(inv_newbasis,np.transpose(u))
        CBM.append(u) # collect the vectors
    CBM = np.transpose(np.array(CBM)) # turn into matrix, transpose so that you get column vectors
    # CBM is the representation of the old basis vectors in the new vector space
    
    newatms = []
    for elem in atms:
        type = elem[3]
        elem = elem[:3]
        elem = np.transpose(np.array(elem))
        newatms.append(np.dot(CBM,elem).tolist())
        newatms[-1].append(type)  # restore type
    return newatms

def buildmore(n,basis,uc):
    numofatoms = n[0]*n[1]*n[2]*len(uc)  #total number of atoms
    atoms = []

    for ix in range(n[0]):
        for seed in uc:
            atoms.append([seed[0]+ix, seed[1], seed[2], seed[3]])
    
    seedmax = len(atoms)
    for iy in range(1,n[1]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1]+iy, atoms[ia][2], atoms[ia][3]])
            
    seedmax = len(atoms)
    for iz in range(1,n[2]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1], atoms[ia][2]+iz, atoms[ia][3]])
    
    basis = np.array(basis)
    for ia in range(len(atoms)):
        temp = (atoms[ia][0]*basis[0] + 
                atoms[ia][1]*basis[1] + 
                atoms[ia][2]*basis[2]).tolist()
        temp.append(atoms[ia][3])    # add the type
        atoms[ia] = temp 
        
    return atoms
     
    
if __name__=="__main__":
    main()

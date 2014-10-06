#!/usr/bin/python

"""example.py

Compute autocorrelation function from several files for a specific E-field value


"""
import numpy as np
#import matplotlib.pyplot as plt
#import pyfftw
import sys
#import Gnuplot, Gnuplot.funcutils
def wait():
	m.getch()

def main():
	"""
	Argument list:
	1-to-x : file names
				The filename is ctritical. It should have the form xxxx_id.yyy
				The prefactor file will use each files' id to find a corresponding match in the prefactor file
	x+1    : pre-factor file
				The prefactor file should have a specific format. Namely, the last entry should be the prefactor
				and the one just before should be it the timestep
	x+2    : which column in the files to use
				!!!0-based!!!
	x+3    : number of rows to consider
				if nrows<1 then all rows are considered
	"""
    numOfFiles = len(sys.argv)-4
    listOfFiles = sys.argv[1:-3]     #names of files
    prfcFile = sys.argv[-3]
    coln = int(sys.argv[-2])
    nrows = int(sys.argv[-1]) 

    print "... Number of files: ", numOfFiles, nrows
    print "... Argument list: "
    for elem in listOfFiles:
       print "...    ", elem
    print "... Prefactor file at: ", prfcFile

    corrOfFiles = []      # a list of arrays
    avcondc = []
    prefacs = np.loadtxt(prfcFile, delimiter=" ")
    
    for k in range(numOfFiles):
        print "K is ================== ", str(k)
        print listOfFiles[k]
        if nrows > 0: correlation_data = autocorr(listOfFiles[k],coln,nrows)
        else: correlation_data = autocorr(listOfFiles[k],coln)

        tmp = listOfFiles[k][::-1]
        rndkey = listOfFiles[k][-tmp.find("_"):-tmp.find(".")] #find the last occurance of "_" and "."
        print ".... key = ", rndkey
        for elem in prefacs:    
            if elem[0] == float(rndkey):
                prfc = float(elem[-1])
                tstep = float(elem[-2])
                break

        conduc = prfc*cumtrapz(correlation_data, tstep*1e-12)
        avcondc.append(conduc)

        corrOfFiles.append(correlation_data)
        fout_corr = open(listOfFiles[k]+".corr",'w')
        print "writing to file... "
        fout_corr.write("#corr k \n")
        fout_corr.write("#(Jm/s)^2 W/mK \n")
        for p in range(len(correlation_data)-1): # N=len(corr_data) => len(conduc)=N-1
            fout_corr.write(str(correlation_data[p])+" "+str(conduc[p])+"\n")
        fout_corr.close()
        print "Done writing to file... ", listOfFiles[k]+".corr"

# average data and output to file
    sigmcond = np.std(avcondc,0)/numOfFiles**0.5
    avcorr = np.sum(corrOfFiles,0)/numOfFiles
    avcondc = np.sum(avcondc,0)/numOfFiles
    time = np.arange(0,len(avcorr))*tstep*1e-12
    n = len(listOfFiles[0])
    fout = open(listOfFiles[0][0:n-tmp.find("_")-2]+".avg",'w') # find last occurance of _ and replace with .avg
    print "saving data to file... " + fout.name
    fout.write("#time correlation conductivity error")
    for p in range(len(avcondc)):
        fout.write("\n" + str(time[p]) + " " + str(avcorr[p]) + " " + str(avcondc[p]) + " " + str(sigmcond[p]))
    fout.close()
    
def autocorr(flname,coln,nrows=-1):

    #### Data extraction ####
    # Load data from files
    # Select only one column
    print 'Processing file... '+ flname + "for " + str(nrows) + " rows"
         
    data = np.loadtxt(flname, skiprows=2 ,delimiter=" ", 
                      usecols=(coln,) )
    if nrows > 0:
    print "++++++croppping++++++"
    print str(len(data))
        data = data[:nrows]
    # get the length of the data
    Nsteps = len(data)
    # pad to a 2^n length with zeros
    Npw2 = 2**np.ceil(np.log2(Nsteps)); Nzrs = Npw2-Nsteps
    print("Old size " + str(len(data)))
    print("Nzrs = " + str(Nzrs))        
    data = np.append(data, [0]*Nzrs, 0)     
    
    print("Padded to new size "+str(len(data))) 
    
    #### FFT ####
    # rfft produces Npw2/2+1 points
    fft_v = np.fft.rfft(data)
    fft_v = np.abs(fft_v)**2 
    # we need to pad with Npw2/2-1 points fft_v[-1] contains terms for n/2 and -n/2
    fft_v = np.append(fft_v, fft_v[-1:1:-1])  # PSD will be the same H(-f)=H*(f)
    data =  np.fft.ifft(fft_v)
    data = np.real(data[:Nsteps-1])
    data = data/np.arange(Nsteps-1,0,-1)
    return data
    
    #######################################################################
    #### IGNORE THESE #####################################################
    #######################################################################
    N = len(R_xx)                       # Reset the number of samples #####
    cumInt = (R_xx[:-1] + R_xx[1:])/2.  # int(f,a,b)=h[f(a)+f(b)]/2   #####
    #10-15 fs / 1e-20 A^2
    cumInt = 1.0e-15*1e-20*np.cumsum(cumInt)/kB*np.mean(prefac[:,4])
    efield = flname.split("_")[2]
    
    # Save data to binary file    
    flname = 'k2_'+ str(efield)    
    flR_xx = 'Rxx2_'+ str(efield)
    fout = open(flR_xx, 'w')
    for item in R_xx:
        fout.write(str(item)+"\n")
    print "saved to " + flR_xx
    fout.close()
    fout = open(flname,'w')
    for item in cumInt:
        fout.write(str(item)+"\n")
    print "saved to " + flname
    
    if ppO!=None:    
        plt.plot(R_xx)
        plt.plot(cumInt)
        ppO.savefig()
    
    return cumInt
    ########################################################################
    
def plot_cumInt(flList):
#    flname = 'k_'+str(flList[0])+'.npy'   
#    d = np.load(flname)    
    for indx in flList:
        flname = 'k_'+str(indx)+'.npy'        
        d = np.load(flname)
        plt.plot(d, label=str(indx)+' V/A')
    plt.legend()
    plt.ylabel("R (K/W)")    
    plt.show()
        
    
def cumtrapz(y,dx):
# returns the cumulative integral of 1-D array using trapezoidal integration
    y = np.array(y)
    inty = y[:-1]+y[1:]
    inty = dx/2*inty
    inty = np.cumsum(inty)
    return inty

if __name__ == "__main__":
    main()

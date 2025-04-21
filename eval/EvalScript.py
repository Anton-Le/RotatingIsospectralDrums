#!/bin/python
import numpy as np
import scipy as sp
import os, shutil, sys, re
from matplotlib import pyplot as plt
from scipy.constants import codata
from mpl_toolkits.mplot3d import Axes3D
''' The script runs a set of numerical experiments
in an automated fashion '''


#default parameters
refDir="refData";
workDir=".";
Nref = range(2,7);
omega = np.concatenate([np.array([0]),2*np.logspace(1,8,50)]);
IOR = [1, 2.42];
Domain = ["ID1", "ID2"]
shift = [ np.array([-0.47619,0.66667]), np.array([-1.0/3,-0.047619]), np.array([-0.99048,0.438095]), np.array([-1.00952,2.70476]) ];
tol = 1e-11;
Neig = 10;

c = sp.constants.value('speed of light in vacuum');
#the maximum distance from the centre of rotation a point can have
R0 = 1;

########################################################################
#EVALUATION OF THE EXPERIMENTS
########################################################################
#1. Load the list of experiments
f = open("ExperimentCases.npy")
ExperimentList = np.load(f);
if workDir == '.':
    workDir = os.path.abspath(workDir);
#2. Go into each data directory
for i in xrange(len(ExperimentList)):
        #create experiment directory
        expDir = workDir + "/Exp_"+str(i);
        #check if directory exists
        if os.path.isdir(expDir) == False:
            print "ERROR! The experiment directory ", expDir, " does not exist"
            sys.exit(1);
        #change into directory
        os.chdir(expDir);
        #dict for storing accumulated differences
        results = {}
        #dict for storing accumulated eigenvalues
        eigvals1 = {}
        eigvals2 = {}
        #iterate over the refinement levels
        for N in ExperimentList[i]['nref']:        
                #2.1 Load the data for both subdomains
                datafile1 = open("output_N%d_%s_1.dat"%(N,ExperimentList[i]['domain']))
                datafile2 = open("output_N%d_%s_2.dat"%(N,ExperimentList[i]['domain']))
                omega = []
                EigVal1 = []
                EigVal2 = []
                lineIdx = 0;
                for line in datafile1:
                    data = line.split();
                    #2.1.1 Load the first row into omega
                    if lineIdx == 0:
                        for el in data:
                            omega.append(float(el));
                    #2.1.2 Load the rest of the rows into the eigenvalue array
                    else:
                        EigValRow = []
                        for el in data:
                            EigValRow.append(float(el.split('+')[0]));
                        EigVal1.append(EigValRow);
                    lineIdx +=1;
                #same procedure for the second domain
                lineIdx = 0;
                for line in datafile2:
                    data = line.split();
                    #skip the omegas
                    #2.1.2 Load the rest of the rows into the eigenvalue array
                    if lineIdx == 0:
                        lineIdx+=1;
                        continue;
                    else:
                        EigValRow = []
                        for el in data:
                            EigValRow.append(float(el.split('+')[0]));
                        EigVal2.append(EigValRow);
                    lineIdx+=1;
                #2.2 For all eigenvalues (using stride 2 due to multiplicity) do
                #2.2.1 Plot |lambda_1 -lambda_2|(omega) vs R0*omega/c
                EVArr1 = np.array(EigVal1);
                EVArr2 = np.array(EigVal2);
                EVDiff = abs(EVArr1[0:-1:2,:]-EVArr2[0:-1:2,:])
                
                omega = np.array(omega);
                #R0 = np.sqrt(3.**2+1**2);
                R0 = 1;
                xVal = abs(R0*omega/c);

                Neig = np.shape(EVDiff)[0];
                #save data
                results["%d"%(N)] = EVDiff;
                eigvals1["%d"%(N)] = EVArr1[0:-1:2,:];
                eigvals2["%d"%(N)] = EVArr2[0:-1:2,:];
                #plot
                fig, host = plt.subplots(figsize=(11,8))
                #skip the first line of data - stationary eigenvalues
                for k in range(Neig/2):
                    host.loglog(xVal[1:], EVDiff[Neig-1-k, 1:], label="$\lambda_%d$"%(k));
                host.loglog(xVal[1:], EVDiff[1,1:], 'k+--',label="$\lambda_9$")
                plt.xlabel("$\\frac{\omega R}{c}$")
                plt.ylabel("$|(\lambda^{(1)})^2 - (\lambda^{(2)})^2|$")
                host.grid();
                plt.title("Deviation of the eigenvalues. meshsize=%g, ior=%g"%(2./2**N, ExperimentList[i]['ior']))
                plt.legend(loc='best')
                plt.show();
                plt.savefig("EVDiff_vs_RotFreq_RefLevel_%d_ior_%g.png"%(N,ExperimentList[i]['ior']),format='png',papertype='a4', orientation='landscape', dpi=600);
                plt.savefig("EVDiff_vs_RotFreq_RefLevel_%d_ior_%g.eps"%(N,ExperimentList[i]['ior']),format='eps',papertype='a4', orientation='landscape', dpi=1000);
                plt.close();
                

                #plot the evolution of the 9th eigenvalue
                fig, host = plt.subplots(figsize=(11,8))
                host.semilogx(xVal[1:], -EVArr1[2,1:], 'r:^', xVal[1:], -EVArr2[2,1:], 'b:^')
                plt.xlabel("$\\frac{\omega R}{c}$")
                plt.ylabel("$|\lambda^2|(\omega)$")
                host.grid();
                plt.title("Behaviour of the analytic eigenvalues. meshsize=%g, ior=%g"%(2./2**N, ExperimentList[i]['ior']))
                plt.legend(("Domain 1", "Domain 2"), loc='best')                
                plt.show();
                plt.savefig("AnalyticEigenvalueEvolution_RefLevel_%d_ior_%g.png"%(N,ExperimentList[i]['ior']),format='png',papertype='a4', orientation='landscape', dpi=600);
                #save as eps
                plt.savefig("AnalyticEigenvalueEvolution_RefLevel_%d_ior_%g.eps"%(N,ExperimentList[i]['ior']),format='eps',papertype='a4', orientation='landscape', dpi=1000);
                plt.close();
       
        #append the angular frequencies
        results['omega'] = omega;
        eigvals1['omega'] = omega;
        eigvals2['omega'] = omega;
        #store results
        savefile = file("Results_Exp_%d.npz"%(i), 'w');
        
        #collect eigenvalue differences into blocks
        Datablob = np.zeros((len(ExperimentList[i]['nref']), np.shape(EVDiff)[0],np.shape(EVDiff)[1]) )
        for k in range(len(ExperimentList[i]['nref'])):
            Datablob[k,:,:] = results["%d"%( ExperimentList[i]['nref'][k] )];
            
        np.savez(savefile, omega=results['omega'], N=ExperimentList[i]['nref'], Data=Datablob );
        savefile.close();
        #collect eigenvalues into blocks
        savefile = file("Eigenvalues_Subdomain1_Exp_%d.npz"%(i), 'w');
        Datablob = np.zeros((len(ExperimentList[i]['nref']), np.shape(EVDiff)[0],np.shape(EVDiff)[1]) )
        for k in range(len(ExperimentList[i]['nref'])):
            Datablob[k,:,:] = eigvals1["%d"%( ExperimentList[i]['nref'][k] )];
        np.savez(savefile, omega=eigvals1['omega'], N=ExperimentList[i]['nref'], Data=Datablob )    ;
        savefile.close();
            
        savefile = file("Eigenvalues_Subdomain2_Exp_%d.npz"%(i), 'w');
        Datablob = np.zeros((len(ExperimentList[i]['nref']), np.shape(EVDiff)[0],np.shape(EVDiff)[1]) )
        for k in range(len(ExperimentList[i]['nref'])):
            Datablob[k,:,:] = eigvals2["%d"%( ExperimentList[i]['nref'][k] )];
        np.savez(savefile, omega=eigvals2['omega'], N=ExperimentList[i]['nref'], Data=Datablob )    ;
        savefile.close();     
       
        #plot the behaviour of the ground state and the 9th eigenstate as functions of N
        GSEVDiff = np.zeros((len(ExperimentList[i]['nref']), len(results['omega'])-1) );
        NinthEVDiff = np.zeros((len(ExperimentList[i]['nref']), len(results['omega'])-1) );
        
        for k in range(len(ExperimentList[i]['nref'])):
            GSEVDiff[k] = results["%d"%(ExperimentList[i]['nref'][k])][Neig-1,1:];
            NinthEVDiff[k] = results["%d"%(ExperimentList[i]['nref'][k])][1,1:];
            
        #prepare for 3D plot
        meshSize = 2.0/np.logspace(min(ExperimentList[i]['nref']),max(ExperimentList[i]['nref']),base=2,num=len(ExperimentList[i]['nref']) );
        xVal = abs(results['omega'][1:]*R0/c );
        X, Y = np.meshgrid(np.log10(xVal), np.log2(meshSize) );
        fig = plt.figure(figsize=(11,8));
        ax = fig.gca(projection='3d');
        ax.plot_wireframe(X,Y, np.log10(GSEVDiff), cstride=0, color='red',label="$\lambda_1$");
        ax.plot_wireframe(X,Y, np.log10(GSEVDiff), rstride=0, color='red', linestyle='dotted');
                
        ax.plot_wireframe(X,Y, np.log10(NinthEVDiff), cstride=0, color='blue',label="$\lambda_9$");
        ax.plot_wireframe(X,Y, np.log10(NinthEVDiff), rstride=0, color='blue', linestyle='dotted');
        xticks = ax.get_xticks()
        xticklabels = ["$10^{%d}$"%(x) for x in xticks]
        yticks = [2.0/(2**x) for x in ExperimentList[i]['nref'] ]
        yticklabels = ["$\\frac{2}{2^{%d}}$"%(x) for x in ExperimentList[i]['nref'] ]
        zticks = ax.get_zticks()
        zticklabels = ["$10^{%d}$"%(x) for x in zticks];
        ax.set_yticks(np.log2(yticks));
        ax.set_xticklabels(xticklabels, size='12');
        ax.set_yticklabels(yticklabels, size='13',  va='center');
        ax.set_zticklabels(zticklabels);
        ax.set_xlabel("$\\frac{R_0\omega}{c}$", size='large', labelpad=8)
        ax.set_ylabel("mesh width", labelpad=6)
        ax.set_zlabel("$|(\lambda^{(1)})^2 - (\lambda^{(2)})^2|$", labelpad=8)
        ax.autoscale_view(tight=True, scalez=False)
        ax.view_init(25,-125)
        ax.legend()
        plt.show()
        plt.savefig("%s_IOR%.3f_Anisospectrality_vs_meshWidth_vs_Parameters.png"%(ExperimentList[i]['domain'], ExperimentList[i]['ior']), format='png');
            
        #go back to the working directory
        os.chdir(workDir);
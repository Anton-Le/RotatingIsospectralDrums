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
samples = 50
omega = np.concatenate([np.array([0]),3*np.logspace(1,8,samples)]);
IOR = [1, 2.42];
Domain = ["ID1", "ID2"]
shift = [ np.array([-0.47619,0.66667]), np.array([-1.0/3,-0.047619]), np.array([-0.99048,0.438095]), np.array([-1.00952,2.70476]) ];

tol = 1e-11;
Neig = 10;
executable='fem_eigenvaluedecomposition_P3_Full_MW04Simplified'


if __name__=='__main__':
    #1. parse command line parameters
    while len(sys.argv)>1:
        option = sys.argv[1]
        del sys.argv[1];
        if option=='-refdir':
            refDir = str(sys.argv[1]);  del sys.argv[1];
        elif option=='-workdir':
            workDir = str(sys.argv[1]);  del sys.argv[1];
        elif option=='-nref':
            rawData = str(sys.argv[1]); del sys.argv[1];
            data = re.split(',', rawData);
            Nref = [int(x) for x in data];
        elif option=='-ior':
            rawData = str(sys.argv[1]); del sys.argv[1];
            data = re.split(',', rawData);
            IOR = [float(x) for x in data];
        elif option == '-domain':
            rawData = str(sys.argv[1]); del sys.argv[1];
            Domain = [str(x) for x in rawData.split(',')];
        elif option == '-runfile':
            rawData = str(sys.argv[1]); del sys.argv[1];
            executable = rawData;
        elif option == '-samples':
            rawData = int(sys.argv[1]); del sys.argv[1];
            samples = rawData;
            
    omega = np.concatenate([np.array([0]),3*np.logspace(1,8,samples)]);
    #Check whether the reference data for the given domains exists
    if refDir == 'refData':
        refDir = os.path.abspath(".")+"/"+refDir;     
    if os.path.isdir(refDir) == False:
        print "Error! The reference directory", refDir, " does not exist";
        sys.exit(1);
        #absolute path to wroking directory
    if workDir == '.':
        workDir = os.path.abspath(workDir);
    #check whether the geometry files are there
    files = os.listdir(refDir);
    for D in Domain:
        Datalist = [D+"_1_Coord.dat", D+"_2_Coord.dat", D+"_1_Elements.dat", D+"_2_Elements.dat",
        D+"_1_Boundary.dat",D+"_2_Boundary.dat"];
        for data in Datalist:
            filepath = refDir +'/'+data;
            if os.path.exists(filepath) == False:
                print "Reference data ", data, "does not exist!"
                sys.exit(1);
    #write the file of angular frequencies
    os.chdir(refDir);
    omegaFile = file("AngularFrequencies.dat", 'w');
    for i in omega:
        omegaFile.write("%.6g\n"%(i));
    omegaFile.close();
    os.chdir(workDir);
    #2. Create a list of all cases to be treated.
    ExperimentList = [];
    for D in Domain:
        for n in IOR:
            #set up experimental data
            experiment = {};
            experiment['domain']=D;
            experiment['ior']=n;
            experiment['omega']=refDir+"/AngularFrequencies.dat";
            experiment['nref']=Nref;
            
            ExperimentList.append(experiment);
       
    #3. Store the list of all cases in a file in the
    #    current working directory
    experimentFile = np.save("ExperimentCases", ExperimentList);
    #4. Run the numerical experiments
    for i in xrange(len(ExperimentList)):
        #create experiment directory
        expDir = workDir + "/Exp_"+str(i);
        if os.path.isdir(expDir):
            print "ERROR! The experiment directory ", expDir, " exists"
            sys.exit(1);
        os.mkdir(expDir);
        #change into directory
        os.chdir(expDir);
        for subID in [1,2]:
            ElementFile = ExperimentList[i]['domain']+"_%d_Elements.dat"%(subID);
            CoordFile = ExperimentList[i]['domain']+"_%d_Coord.dat"%(subID);
            BoundaryFile = ExperimentList[i]['domain']+"_%d_Boundary.dat"%(subID);
            OmegaFile = ExperimentList[i]['omega']
            #link the mesh files from the refernce directory
            #link the coordinate file
            src = refDir+"/"+CoordFile;
            dst = expDir+"/"+CoordFile;
            os.symlink(src,dst);
            #link the boundary file
            src = refDir+"/"+BoundaryFile;
            dst = expDir+"/"+BoundaryFile;
            os.symlink(src,dst);
            #link the element file
            src = refDir+"/"+ElementFile;
            dst = expDir+"/"+ElementFile;
            os.symlink(src,dst);
            #link the angular frequencies file
            dst = expDir+"/AngularFrequencies.dat";
	    if os.path.isfile(dst) == False:
	            os.symlink(OmegaFile,dst);
            #create parameter file
            params = '''CF=%s
EF=%s
BFD=%s
n=%f
omega=%s
tol=%g
FP=%s
Neig=%d
xShft=%f
yShft=%f\n'''%(CoordFile, ElementFile, BoundaryFile, ExperimentList[i]['ior'],\
            "AngularFrequencies.dat", tol, D+"_"+str(subID), Neig, shift[subID-1][0],\
            shift[subID-1][1]);
            #write parameter file
            paramFile = file("parameters.txt", 'w')
            paramFile.write(params)
            paramFile.close()           
            #run experiment
            for N in ExperimentList[i]['nref']:
                cmd = refDir+"/"+executable + " --N "+str(N);
                failure = os.system(cmd);
                if failure:
                    print "Simulation with N =",N, " failed. Exiting..."
                    sys.exit(1);

                #rename output
                os.rename("output.dat", "output_N%d_%s_%d.dat"%(N,ExperimentList[i]['domain'], subID))
        #change back to working directory
        os.chdir(workDir);

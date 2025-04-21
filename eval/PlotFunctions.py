from matplotlib import pyplot as plt
import numpy as np
import matplotlib.tri as tri
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import AxesGrid

def readEigenmodes(filename):
    f = open(filename,'r');
    eigval = [];
    for line in f:
        data = line.split();
        ev = [];
        for x in data:
            x = x.replace( '(', ' ' );
            x = x.replace( ",-", '-' );
            x = x.replace( ',', '+' );
            x = x.replace( ')', 'j' );
            ev.append(complex(x));
        eigval.append(ev);
    f.close()
    return np.copy( np.array(eigval) );

def readTriangulation(vertexFilename, elementFilename):
    vertFile = open(vertexFilename, 'r')
    vertices = []
    for line in vertFile:
        data = line.split();
        vertex = [float(data[0]), float(data[1]) ];
        vertices.append(vertex);
    
    vertFile.close()

    triFile = open(elementFilename, 'r')
    triangles = []
    for line in triFile:
        data = line.split();
        triangle = [int(data[0]), int(data[1]), int(data[2])]
        triangles.append(triangle)

    triFile.close()
    
    vertexArray = np.asarray(vertices)
    triangleArray = np.asarray(triangles)
    triangulationOwn = tri.Triangulation(vertexArray[:,0], vertexArray[:,1], triangles=triangleArray)
    return triangulationOwn;

def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    return at

def PlotEigenmodes(Modes, ModeDiff, triangulation, mode_cbar_normalizer, modeDiff_cbar_normalizer,
    xticklabels=[], yticklabels=[], show=False, save=True, title = None, filename="Eigenmodes.png",
    drawOrigin=False):
    #fetch the number of eigenmodes and check whether the numbers match for modes and 
    #their deviations from the static case
    Nmodes = np.shape(Modes)[1]
    assert( Nmodes == np.shape(ModeDiff)[1] );
    
    fig = plt.figure();
    grid = AxesGrid(fig, 111, nrows_ncols=(2,3), axes_pad=(0.0,0.0),\
                share_all=True, label_mode="1", cbar_location="right", cbar_mode="edge")
    #set aspect ratio
    for k in range(2*Nmodes):
        grid[k].set_aspect('equal');
    #draw
    for k in range(Nmodes):
        im1 = grid[k].tripcolor(triangulation, np.real(Modes[:,k]), shading='gouraud', cmap=plt.cm.viridis, norm=mode_cbar_normalizer)
        im2 = grid[k+Nmodes].tripcolor(triangulation, np.real(ModeDiff[:,k]), shading='gouraud', cmap=plt.cm.viridis, norm=modeDiff_cbar_normalizer)
    
    #set color bars
    grid.cbar_axes[0].colorbar(im1);
    grid.cbar_axes[1].colorbar(im2);
    #set the ticks, if provided
    if len(xticklabels) != 0:
        grid[0].set_xticklabels(xticklabels);
    if len(yticklabels) != 0:
        grid[0].set_yticklabels(yticklabels);
        
    #create and add labels to the plot to be used later on in referring
    #the modes to partiular values of the angular frequency
    imLabels = [];
    for k in range(Nmodes):
        imLabels.append("$\omega_%d$"%(k+1))
    
    for ax, im_title in zip(grid.axes_row[0], imLabels):
        t = add_inner_title(ax, im_title, loc=1)
        t.patch.set_alpha(0.5)
        
    if title != None:
        grid[1].set_title(title);
    if drawOrigin == True:
            circ = mpatches.Circle((0,0), 0.1, ec=None, fc='red', alpha=0.6);
            grid[0].add_artist(circ);
            
    if show == True:
        plt.show();
    if save == True:
        fig.savefig(filename, format='png', dpi=1000);

def readEigenvalueData(filenameDifferences, filenameEigenvalues=None):
    data = np.load(filenameDifferences);
    EVDiff = data['Data'];
    omega = data['omega'];
    EVDiffVal = EVDiff[-1];
    if filenameEigenvalues != None:
        data = np.load(filenameEigenvalues);
        EV = data['Data'];
        Eigenvalues = EV[-1];
        return [omega, EVDiffVal, Eigenvalues];        
    else:
        return [omega, EVDiffVal];
        

def Plot4(x, EVDiff, Eigval1, Eigval2, GSIdx, NSIdx, ls1, ls2, mark1, mark2, color1, color2):
    fig, axesArray = plt.subplots(2, sharex=True)
    
    axesArray[0].loglog(x, EVDiff[GSIdx], color=color1, linestyle=ls1, marker=mark1)
    axesArray[0].loglog(x, EVDiff[NSIdx], color=color1, linestyle=ls1, marker=mark2)
    
    axesArray[0].legend(("$\lambda_1$", "$\lambda_9$"), loc='best', framealpha=0.7)
    
    axesArray[1].semilogx(x, Eigval1[GSIdx], linestyle=ls1, marker=mark1, color=color1, label="Domain 1" );
    axesArray[1].semilogx(x, Eigval2[GSIdx], linestyle=ls2, marker=mark1, color=color2, label="Domain 2" );
    
    axesArray[1].semilogx(x, Eigval1[NSIdx], linestyle=ls1, marker=mark2, color=color1)
    axesArray[1].semilogx(x, Eigval2[NSIdx], linestyle=ls2, marker=mark2, color=color2)
    
    axesArray[1].legend(loc='center left', framealpha=0.7);
    #grid lines
    axesArray[0].grid()
    axesArray[1].grid()
    axesArray[0].set_xlim(min(x), max(x));
    
    fig.subplots_adjust(hspace=0.3);
    
    #add shading
    axesArray[0].axvspan(5*1e-1, max(x), facecolor='none', hatch='x', alpha=0.1)
    axesArray[1].axvspan(5*1e-1, max(x), facecolor='none', hatch='x', alpha=0.1)
    
    #labels
    axesArray[0].set_ylabel("$|(\lambda^{(1)})^2 - (\lambda^{(2)})^2|$", size='large')
    axesArray[1].set_ylabel("$\lambda^2$", size='large')
    axesArray[1].set_xlabel("$\\frac{R_0\omega}{c}$", size='xx-large')
    
    fig.savefig("ID1_Anisospectrality_And_Eigenvalues_Full_Vacuum.eps", format='eps', dpi=resolution);
    #plt.show()

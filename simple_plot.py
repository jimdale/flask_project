import os
import glob
import numpy as np
import scipy
import matplotlib
import matplotlib.pylab as plt
import datetime
import time
import random
import astropy
from astropy.io import fits
from astropy import units as u
from astropy import table
import pdb
plt.rcParams['figure.figsize'] = (12,8)


def plotData(NQuery,FigureStrBase,variables,xMin,xMax,yMin,yMax,zMin,zMax) :
#def plotData(NQuery,FigureStrBase,xMin,xMax,yMin,yMax,zMin,zMax) :
    """
    This is where documentation needs to be added

    Parameters
    ----------
    NQuery
    FigureStrBase : str
        The start of the output filename, e.g. for "my_file.png" it would be
        my_file
    xMin
    xMax
    yMin
    yMax
    zMin
    zMax
    """
    
    # deletes plotting screen
    plt.clf()

    # reads entire table
    d = table.Table.read("merged_table.ipac", format='ascii.ipac')
    # Assigns columns to separate variables
    # TODO: make it more general, to have arbitrary plots
    Author = d['Names']
    Run = d['IDs']
    print variables[0]
    x_ax = d[variables[0]]
    y_ax = d[variables[1]]
    z_ax = d[variables[2]]
    IsSim = (d['IsSimulated'] == 'True')
    
    # selects surface density points wthin the limits
    Use_x_ax = (x_ax > xMin) & (x_ax < xMax)
    Use_y_ax = (y_ax > yMin) & (y_ax < yMax)
    Use_z_ax = (z_ax > zMin) & (z_ax < zMax)
    # intersects the three subsets defined above
    Use = Use_x_ax & Use_y_ax & Use_z_ax
    Obs = (~IsSim) & Use
    Sim = IsSim & Use
    
    # defines a set of authors of the selected datapoints
    UniqueAuthor = set(Author[Use])
    NUniqueAuthor = len(UniqueAuthor)
    
    # TODO: add controls on the colors 
    #colors = random.sample(matplotlib.colors.cnames, NUniqueAuthor)
    colors = list(plt.cm.jet(np.linspace(0,1,NUniqueAuthor)))
    random.seed(12)
    random.shuffle(colors)
    
    # sets log axes
    plt.loglog()
    # TODO: add controls on the symbols
    # sets round markers for obs's and square m for sim's
    markers = ['o','s']
    for iAu,color in zip(UniqueAuthor,colors) :
        UsePlot = (Author == iAu) & Use
        ObsPlot = ((Author == iAu) & (~IsSim)) & Use 
        SimPlot = ((Author == iAu) & (IsSim)) & Use
        # plots colored markers for obs's and sim's
        if any(ObsPlot):
            plt.scatter(x_ax[ObsPlot], y_ax[ObsPlot], marker=markers[0],
                        s=(np.log(np.array(z_ax[ObsPlot]))-np.log(np.array(zMin))+0.5)**3.,
                        color=color, alpha=0.5)
        if any(SimPlot):
            plt.scatter(x_ax[SimPlot], y_ax[SimPlot], marker=markers[1],
                        s=(np.log(np.array(z_ax[SimPlot]))-np.log(np.array(zMin))+0.5)**3.,
                        color=color, alpha=0.5)
    # draws black contour
    if any(Obs):
        plt.scatter(x_ax[Obs], y_ax[Obs], marker=markers[0],
                    s=(np.log(np.array(z_ax[Obs]))-np.log(np.array(zMin))+0.5)**3.,
                    facecolors='none', edgecolors='black',
                    alpha=0.5)
    if any(Sim):
        plt.scatter(x_ax[Sim], y_ax[Sim], marker=markers[1],
                    s=(np.log(np.array(z_ax[Sim]))-np.log(np.array(zMin))+0.5)**3.,
                    facecolors='none', edgecolors='black',
                    alpha=0.5)
	
    plt.xlabel('$\Sigma$ [M$_{\odot}$ pc$^{-2}$]', fontsize=16)
    plt.ylabel('$\sigma$ [km s$^{-1}$]', fontsize=16)

    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # adding fake points to show the size
    axes_limits = ax.axis()
    xax_limits = axes_limits[:2]
    yax_limits = axes_limits[2:]

    # TODO: write a function with this section
    # TODO: change position based on user input
    xfake = [0.1,0.1,0.1] #[xax_limits[0] + xax_limits[0]*2.,xax_limits[0] + xax_limits[0]*2.,xax_limits[0] + xax_limits[0]*2.]
    yfake = [0.85,0.9,0.95,] #[yax_limits[1] - yax_limits[1]*0.01,yax_limits[1] - yax_limits[1]*0.3,yax_limits[1] - yax_limits[1]*0.6]
    radius = np.array([1e-1,1e0,1e1]) #*u.pc #(zMin + zMax)*0.5
    print xfake, yfake
    print radius, zMin, zMax

    # Put a legend to the right of the current axis
    ax.legend(UniqueAuthor, loc='center left', bbox_to_anchor=(1.0, 0.5), prop={'size':12}, markerscale = .7, scatterpoints = 1)

    plt.scatter(np.array(xfake), np.array(yfake), marker='s',
	      s=(np.log(np.array(radius))-np.log(np.array(zMin.value))+0.5)**3., transform=ax.transAxes,
	      facecolors='g')
    for xf,yf,rad in zip(xfake,yfake,radius):
      plt.text(xf + 0.05,yf-0.01, str(rad) + ' ' + str(zMin.unit), transform=ax.transAxes)

    #plt.xlim((xMin.to(u.M_sun/u.pc**2).value,xMax.to(u.M_sun/u.pc**2).value))
    #plt.ylim((yMin.to(u.km/u.s).value,yMax.to(u.km/u.s).value))
    #print xMin.value,xMax.value
    plt.xlim(xMin.value,xMax.value)
    plt.ylim(yMin.value,yMax.value)
    plt.savefig(FigureStrBase+'.png',bbox_inches='tight',dpi=150)
    plt.savefig(FigureStrBase+'.pdf',bbox_inches='tight',dpi=150)
    plt.show()
    
def clearPlotOutput(FigureStrBase,TooOld) :
    
    for fl in glob.glob(FigureStrBase+"*.png") + glob.glob(FigureStrBase+"*.pdf"):
        now = time.time()
        if os.stat(fl).st_mtime < now - TooOld :
            os.remove(fl)
    
def timeString() :
    
    TimeString=datetime.datetime.now().strftime("%Y%m%d%H%M%S%f")
    return TimeString
  
if __name__ == "__main__":
  # TODO: change units according to the axes
  xMin = 1e-1*u.M_sun/u.pc**2
  xMax = 1e5*u.M_sun/u.pc**2
  yMin = 1e-1*u.km/u.s
  yMax = 3e2*u.km/u.s
  zMin = 1e-2*u.pc
  zMax = 1e3*u.pc

  variables = ['SurfaceDensity','VelocityDispersion','Radius']
  print variables
  FigureStrBase = ''
  for var in variables:
    FigureStrBase += var + '_'
  FigureStrBase = FigureStrBase[0:-1]
  NQuery=timeString()
  TooOld=300

  clearPlotOutput(FigureStrBase,TooOld)

  #FigureStrBase='Output_Sigma_sigma_r_'
  plotData(NQuery,FigureStrBase,variables,xMin,xMax,yMin,yMax,zMin,zMax)

  #d.show_in_browser(jsviewer=True)



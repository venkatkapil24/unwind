import numpy as np
import sys, ast
import xmltodict as x2d
from scipy.interpolate import CubicSpline

def input_A( path2inputfile ):
    """Imports the drift matrix from the i-pi input file."""
    try:
        #Parses the inputfile as dictionary.
        with open(path2inputfile) as x:
            xmldict = x2d.parse(x.read())
        #Obtains the elements and shape of the array.
        #TODO: Implement units conversion
        arrel =  np.asarray(ast.literal_eval(xmldict["simulation"]["system"]["motion"]["dynamics"]["thermostat"]["A"]["#text"]))
        arrsh =  np.asarray(ast.literal_eval(xmldict["simulation"]["system"]["motion"]["dynamics"]["thermostat"]["A"]["@shape"]))
        #Returns the reshaped array.
        return  arrel.reshape(arrsh)*41.341373
    except:
        print("Drift matrix not found.", sys.exc_info()[0])
        raise

def input_C( path2inputfile ):
    """Imports the drift matrix from the i-pi input file."""
    try:
        #Parses the inputfile as dictionary.
        with open(path2inputfile) as x:
            xmldict = x2d.parse(x.read())
        #Obtains the elements and shape of the array.
        #TODO: Implement units conversion
        arrel =  np.asarray(ast.literal_eval(xmldict["simulation"]["system"]["motion"]["dynamics"]["thermostat"]["C"]["#text"]))
        arrsh =  np.asarray(ast.literal_eval(xmldict["simulation"]["system"]["motion"]["dynamics"]["thermostat"]["C"]["@shape"]))
        #Returns the reshaped array.
        return  arrel.reshape(arrel.reshape(arraysh))
    except:
        try:
        #Parses the inputfile as dictionary.
            with open(path2inputfile) as x:
                xmldict = x2d.parse(x.read())
        #Assumes canonaical sampling and computes the covariance matrix.
            arrsh =  np.asarray(ast.literal_eval(xmldict["simulation"]["system"]["motion"]["dynamics"]["thermostat"]["A"]["@shape"]))
            temperature = float(ast.literal_eval(xmldict["simulation"]["system"]["ensemble"]["temperature"]["#text"]))
            temperature = temperature/3.1577464e5
        #TODO: Implement units conversion
            return np.eye(arrsh[-1])
        except:
            print("Temperature not found.", sys.exc_info()[0])
            raise

def input_vvac(path2inputfile, mrows, stride):
    """Imports the vvac file and extracts the ."""
    try:
        dvvac=np.genfromtxt(path2inputfile, usecols=((2,3)))
        if( mrows == -1 ):
            mrows = len(dvvac)
        return dvvac[:mrows][::stride]
    except:
        print("error in inporting the vvac file",  sys.exc_info()[0])
        raise

def output_vvac(xy,path2outfile, refvvac, action):
    """Imports the vvac file and extracts the ."""
    try:
        xorg=refvvac[:,0]
        xred=xy[0]
        yred=xy[1]
        if(action == "org"):
            x=xorg
            y=CubicSpline(xred, yred, extrapolate=True)(xorg)
            np.savetxt(path2outfile,np.vstack((x, y)).T)
        elif(action == "red"):
            x=xred
            y=yred
            np.savetxt(path2outfile,np.vstack((x, y)).T)
    except:
        print("error in printing the vvac",  sys.exc_info()[0])
        raise

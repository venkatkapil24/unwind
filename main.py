#  __________________________________________
# |  __         __             _           _ |
# | / /         \ \           (_)         | ||
# || |_   _ _ __ | | __      ___ _ __   __| ||
# || | | | | '_ \| | \ \ /\ / / | '_ \ / _` ||
# || | |_| | | | | |  \ V  V /| | | | | (_| ||
# || |\__,_|_| |_| |   \_/\_/ |_|_| |_|\__,_||
# | \_\         /_/                          |
# |__________________________________________|

# Written by Venkat Kapil <venkat.kapil[at]gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                          
"""Given the parameters of a Generalized Langevin Equation and the vibrational density of states predicts the velocity-velocity autcorrelation obtained by the dynamics. Conversely, given the velocity-velocity autocorrelation function removes the disturbance affected by the thermostat and returns the underlying vibrational density of states."""

import argparse,sys
from glelibs import *
from gleio import *

# adds description of the program.
parser=argparse.ArgumentParser(description="Given the parameters of a Generalized Langevin Equation and the vibrational density of states predicts the velocity-velocity autcorrelation obtained by the dynamics. Conversely, given the velocity-velocity autocorrelation function removes the disturbance affected by the thermostat and returns the underlying vibrational density of states. ")

# adds arguments.
parser.add_argument("-a","--action", nargs=1, choices=["conv","deconv"], default=None, help="choose conv if you want to obtain the response of the thermostat on the vibrational density of states; choose deconv if you want to obtain the micro-canonical density of states by removing the disturbance induced by the thermostat")
parser.add_argument("-iipi", "--input_ipi", nargs=1, type=str, default=None, help="the relative path to the i-PI inputfile")
parser.add_argument("-ivvac", "--input_vvac", nargs=1, type=str, default=None, help="the relative path to the input velocity-velocity autocorrelation function")
parser.add_argument("-k", "--input_kernel", nargs=1, type=str, default=[None], help="the relative path to the kernel function")
parser.add_argument("-mrows", "--maximum_rows", nargs=1, type=int, default=[-1], help="the maximum number of rows to be imported from INPUT_VVAC")
parser.add_argument("-s", "--stride", nargs=1, type=int, default=[1], help="the stride for computing the kernal")
parser.add_argument("-ovvac", "--output_vvac", nargs=1, type=str, default=["output-vvac.data"], help="the name of the output file containing the (de)convoluted spectrum")
parser.add_argument("-oflag","--output_flag", nargs=1, choices=["org","red"], default=["red"], help="choose orig_grid if you want OUTPUT_VVAC to have the same stride as INPUT_VVAC; choose reduced_grid if you want the OUTPUT_VVAC to have a stride of STRIDE")

# parses arguments.
if( len(sys.argv) > 1):
    args=parser.parse_args()
else:
    parser.print_help()
    sys.exit()

# stores the arguments
path2iipi=str(args.input_ipi[-1])
path2ivvac=str(args.input_vvac[-1])
path2ker=args.input_kernel[-1]
path2ovvac=str(args.output_vvac[-1])
oflag=str(args.output_flag[-1])
action=str(args.action[-1])
nrows=int(args.maximum_rows[-1])
stride=int(args.stride[-1])

print
print "#  __________________________________________ "
print "# |  __         __             _           _ |"
print "# | / /         \ \           (_)         | ||"
print "# || |_   _ _ __ | | __      ___ _ __   __| ||"
print "# || | | | | '_ \| | \ \ /\ / / | '_ \ / _` ||"
print "# || | |_| | | | | |  \ V  V /| | | | | (_| ||"
print "# || |\__,_|_| |_| |   \_/\_/ |_|_| |_|\__,_||"
print "# | \_\         /_/                          |"
print "# |__________________________________________|"
print
print

# imports input spectrum and the GLE parameters
print "# importing the parameters of the GLE."
Ap = input_A(path2iipi)
Cp = input_C(path2iipi)
Dp = np.dot(Ap,Cp) + np.dot(Cp,Ap.T)
print "# importing the input spectrum."
ivvac=input_vvac(path2ivvac, nrows, stride)
ix=ivvac[:,0]
iy=ivvac[:,1]

# computes the vvac kernel
if (path2ker == None):
    print "# computing the kernel."
    ker = gleKernel(ix, Ap, Dp)
    np.savetxt("ker.data", ker)
else:
    print "# importing the kernel."
    ker=np.loadtxt(path2ker)

# (de-)convolutes the spectrum
if(action == "conv"):
    print "# printing the output spectrum."
    output_vvac((ix, np.dot(iy,ker.T)), path2ovvac, input_vvac(path2ivvac, nrows, 1), oflag)
elif(action == "deconv"):
    print "# deconvoluting the input spectrum."
    oy, ocnvg = ISRA(ix, ker, iy)
    output_vvac((ix, oy), path2ovvac, input_vvac(path2ivvac, nrows, 1), oflag)
    np.savetxt("cnvg", ocnvg)

from os import sys
from tepy.entanglement.fidelity import *

filename = sys.argv[1]
total_dim = int(sys.argv[2])
outputname = sys.argv[3]
time_diff = int(sys.argv[4])
load_bool = sys.argv[5]
save_bool = sys.argv[6]
if sys.argv[7]:
   line_nums = int(sys.argv[7])
   lines = sys.argv[8:8+line_nums]
   vlines = [float(line) for line in lines]
if load_bool == 'T':
   load = True
else:
   load = False
if save_bool == 'T':
   save = True
else:
   save = False
   
plot_fidelity( total_dim=total_dim, datafile=filename, 
               savename=outputname, time_diff=time_diff, vlines=vlines, 
               belltype='phiplus', qutip_method=False, save=save, load=load)

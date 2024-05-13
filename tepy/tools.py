import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.fft import fft, fftfreq
from tepy.globalvars import *

def get_data( filename='POPULATIONS.dat' ):
   dat = np.loadtxt( filename )
   return dat
 
def get_py_indx_rhofile( istate, jstate, imag_or_real, total_states ):
   """
   Get the python index number from a TimeESR density matrix output
   based on the chosen i state and j state, i.e. the index that gets
   rho(i,j), both real and imaginary
   
   State number should be given in 1 base

   Example for 3 state system
   python index  0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
   column type   t Re 11 Re 21 Re 31 Re 12 Re 22 Re 32 Re 13 Re 23 Re 33 Im 11 Im 21 Im 31 Im 12 Im 22 Im 32 Im 13 Im 23 Im 33
   
   Converted to python index format (note the backwards look to rho_01 and how we access it -- it's due to Fortran2018's printing of arrays)
   """
   if imag_or_real == 'imag':
      indx = (total_states*total_states) + ( ( jstate - 1 ) * total_states) + istate
   elif imag_or_real == 'real':
      indx = ( ( jstate - 1 ) * total_states) + istate
   
   return indx
 
def plot_2d( ax_obj, x, y, xfact=1, yfact=1, swapxy=False, linewidth=1.5, marker='', color='k', label='', xlim=None, ylim=None, name=None, markersize=6 ):
   # limits are in the form [min, max, diff]   

   if swapxy:
      xplot, yplot = y * float(yfact), x * float(xfact)
   else:
      xplot, yplot = x * float(xfact), y * float(yfact)
   
   if xlim:
      ticks = np.arange(xlim[0],xlim[1]+(xlim[2]/100),xlim[2])
      ax_obj.set_xlim(xlim[0],xlim[1])
      ax_obj.set_xticks(ticks)

   if ylim:
      ticks = np.arange(ylim[0],ylim[1]+(ylim[2]/100),ylim[2])
      ax_obj.set_ylim(ylim[0],ylim[1])
      ax_obj.set_yticks(ticks)

   ax_obj.plot( xplot, yplot, linewidth=linewidth, marker=marker, label=name, markersize=markersize )
   
def get_envelope( signal, height=0.8, distance=10 ):
   from scipy.signal import find_peaks
   peaks, _ = find_peaks(signal, height=height, distance=distance)
   return peaks # peaks are locations

def get_neg_exp_curve( x, y, ymin ):
   from scipy.optimize import curve_fit
   
   xmin_indx = np.argmax(y)

   x_data, y_data = [], []
   for i, yval in enumerate(y[xmin_indx:]):
      if (yval >= ymin and not np.isnan(yval)):
         x_data.append(x[i])
         y_data.append(yval)

   x_dat, y_dat = np.array(x_data), np.array(y_data)

   # define type of function to search
   def model_func(x, a, k, b):
       return a * np.exp(-k*x) + b

   #p0 = (1.0,1./2000.,0.8) # starting search coefficients
   opt, pcov = curve_fit(model_func, x_dat, y_dat)#, p0)
   a, k, b = opt
   x_curve = np.linspace(x_dat[0],x_dat[-1],100)
   y_curve = model_func(x_curve, a, k, b)

   tau = 1./k
   
   residuals = y_dat - model_func(x_dat, *opt)
   ss_res = np.sum(residuals**2)
   ss_tot = np.sum((y_dat - np.mean(y_dat))**2)
   r_squared = 1.0 - (ss_res/ ss_tot)
   
   print("R squared value: ",r_squared)

   return x_curve, y_curve, tau, a, k, b
   
def get_neg_linear_curve( x, y, ymin ):
   from scipy.optimize import curve_fit
   
   xmin_indx = np.argmax(y)

   x_data, y_data = [], []
   for i, yval in enumerate(y[xmin_indx:]):
      if (yval >= ymin and not np.isnan(yval)):
         x_data.append(x[i])
         y_data.append(yval)

   x_dat, y_dat = np.array(x_data), np.array(y_data)

   # define type of function to search
   def model_func(x, k, b):
       return -k * x + b

   #p0 = (1.0,1./500.,0.2) # starting search coefficients
   opt, pcov = curve_fit(model_func, x_dat, y_dat)#, p0)
   k, b = opt
   x_curve = np.linspace(x_dat[0],x_dat[-1],100)
   y_curve = model_func(x_curve, k, b)

   residuals = y_dat - model_func(x_dat, *opt)
   ss_res = np.sum(residuals**2)
   ss_tot = np.sum((y_dat - np.mean(y_dat))**2)
   r_squared = 1.0 - (ss_res/ ss_tot)
   
   print("R squared value: ",r_squared)

   return x_curve, y_curve, k, b
 
def extract_peak_frequency(data, start, end, steps):
   """
   Get the peak frequency from a data set
   """
   sample_rate = (end-start)/steps
   fftfrequencies = fftfreq(steps,sample_rate)[:steps//2]
   fft_data = (2.0/steps) * np.abs(fft(data)[:steps//2])
   peak_coefficient = np.argmax(np.abs(fft_data))
   peak_freq = fftfrequencies[peak_coefficient]

   return abs(peak_freq)
    
def get_fft_spectra(data, start, end, steps):
   """
   Get the entire FFT spectra from a data set
   """
   sample_rate = (end-start)/steps
   fftfrequencies = fftfreq(steps,sample_rate)[:steps//2]
   fft_data = (2.0/steps) * np.abs(fft(data)[:steps//2])

   return fft_data, fftfrequencies
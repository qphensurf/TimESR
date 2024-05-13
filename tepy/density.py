import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tepy.globalvars import *
from tepy.tools import *

def get_dm( rho, dt, total_dim ):
   
   num_dim=4
   
   real_indices = np.array([ [get_py_indx_rhofile( i+1, j+1, 'real', total_dim ) for j in range(num_dim)] for i in range(num_dim) ])
   imag_indices = np.array([ [get_py_indx_rhofile( i+1, j+1, 'imag', total_dim ) for j in range(num_dim)] for i in range(num_dim) ])
   
   rho_data = np.zeros([num_dim,num_dim],dtype=complex)
   for i in range(num_dim):
      for j in range(num_dim):
         rho_data[i,j] = rho[dt,real_indices[i,j]] + 1j*rho[dt,imag_indices[i,j]]
   
   return rho_data

def get_dm_sub( rho, dt ):
   
   num_dim = 4
   
   rho_data = np.zeros([num_dim,num_dim],dtype=complex)
   for i in range(num_dim):
      for j in range(num_dim):
         rho_data[i,j] = rho[dt,i,j]   
   return rho_data 

def transform_dm( rho, u_matrix ):
   """
   Transform the density matrix using the unitary matrix
   """  
   rho_new = np.matmul(u_matrix,rho)
   return np.matmul(rho_new,np.transpose(np.conjugate(u_matrix)))
   
def transform_dm_rot( rho, time, rot_freq, freq_width, nfreq ):
   """
   Rotate the density matrix at a particular frequency about the z-axis
   """  
   import copy
   
   # The rotation frequency is not so easy to determine for TimeESR systems
   # so instead, we use the user-supplied frequency as the mid-point from which to test, 
   # plus or minus the range requested
   # We then perform a FFT of the off-diagonal element (the one that typically shows
   # the oscillations) and find the value for which the frequency is minimized to 0.
   if nfreq != 1:
      print( "Testing rotation frequency of: "+str(rot_freq)+" GHZ and width "+str(freq_width)+" GHZ for "+str(nfreq)+" points" )
      domega = 2.0*freq_width/float(nfreq-1)
      fft_peak_table = np.zeros((nfreq,2),dtype=np.cdouble)
      for ifreq in range(nfreq):
         omega = (rot_freq - freq_width) + (float(ifreq)*domega)
         print("   Working on frequency "+str(omega)+" GHZ")
         # Rho should be in the format of rho(ntimepoints,ndim,ndim)
         rho_sub = copy.deepcopy(rho)
         # Form unitary
         urot = np.zeros((len(time),2,2),dtype=np.cdouble)
         urot[:,0,0] = np.exp( -1.0j * 2.0*np.pi * 0.0 * time[:] )
         urot[:,1,1] = np.exp( -1.0j * 2.0*np.pi * omega * time[:] )
         # Perform rotation at each time step
         for i in range(len(time)):
            rho_sub[i,:,:] = transform_dm( rho_sub[i,:,:], urot[i,:,:] )
         # Tabulate frequencies obtained
         fftoffdiag_peak = extract_peak_frequency(np.real(rho_sub[:,0,1]), time[0], time[-1], len(time))
         fft_peak_table[ifreq,0] = omega
         fft_peak_table[ifreq,1] = fftoffdiag_peak
      # Return all indices with the smallest value
      min_freq_indices = np.flatnonzero(fft_peak_table[:,1] == fft_peak_table[:,1].min())
      # If the size of the array is larger than 1, pick the 'middle' or the 'middle-1' value
      if (len(min_freq_indices) == 1):
         index_loc = 0
      elif (len(min_freq_indices)%2 == 0):
         index_loc = int(np.floor(len(min_freq_indices)/2)) - 1
      else:
         index_loc = int(np.floor(len(min_freq_indices)/2))
      min_freq_index = min_freq_indices[index_loc]   
      min_freq = fft_peak_table[min_freq_index,0]
      min_freq_val = fft_peak_table[min_freq_index,1]
      print("Found minimal oscillations at rotation frequency "+str(np.real(min_freq))+" GHZ with peak frequency "+str(np.real(min_freq_val)))
      print("Table of frequencies and maximum peaks")
      print(fft_peak_table)
   else:
      min_freq = rot_freq
   # Perform the 'best' rotation found in the table
   rho_sub = copy.deepcopy(rho)
   # Form unitary
   urot = np.zeros((len(time),2,2),dtype=np.cdouble)
   urot[:,0,0] = np.exp( -1.0j * 2.0*np.pi * 0.0 * time[:] )
   urot[:,1,1] = np.exp( -1.0j * 2.0*np.pi * min_freq * time[:] )
   # Perform rotation at each time step
   for i in range(len(time)):
      rho_sub[i,:,:] = transform_dm( rho_sub[i,:,:], urot[i,:,:] )
      
   return rho_sub

def plot_dmelements( num_dim=4, total_dim=8, datafile='RHO.dat', time_diff=None, vlines=None, rot=False, 
                     rot_freq=0.0, state1=1, state2=2, freq_width=0.0, nfreq=100 ):

   xlabel = r'Time (ns)'
   xfact = 1
   yfact = 1
   
   matplotlib.rc('xtick', labelsize=labelsize)
   matplotlib.rc('ytick', labelsize=labelsize)
   plt.rcParams['figure.figsize'] = (fig_width_inches, fig_height_inches)
   plt.rcParams['figure.dpi'] = dpi

   prop_cycle = plt.rcParams['axes.prop_cycle']
   colors = prop_cycle.by_key()['color']
   
   rho = get_data( filename=datafile )
   time = rho[:,0]
   xmin = time[0] 
   xmax = time[-1]
   if time_diff:
      xdiff = time_diff
   else:
      xdiff = time[1] - time[0]

   if rot:
      rho_select = np.zeros((len(rho[:,0]),2,2),dtype=np.cdouble)
      list_states = [state1,state2]
      for i in range(2):
         for j in range(2):
            idx = list_states[i]
            idy = list_states[j]
            real_index = get_py_indx_rhofile( idx, idy, 'real', total_dim )
            imag_index = get_py_indx_rhofile( idx, idy, 'imag', total_dim )
            rho_select[:,i,j] = rho[:,real_index] + 1.0j*rho[:,imag_index]
      rho_plot = transform_dm_rot( rho_select, time, rot_freq, freq_width, nfreq )
      
      ylabel = r'Re[$\rho_{ij}$]'
      fig, ax = plt.subplots()
      fig.set_size_inches( fig_width_inches, fig_height_inches)
      fig.set_dpi(dpi)
      ax.tick_params(direction='in')
      ax.set_xlabel(xlabel, fontsize=fontsize)
      ax.set_ylabel(ylabel, fontsize=fontsize)

      for i in range(2):
         for j in range(2):
            state_a = list_states[i]
            state_b = list_states[j]
            color_idx = 2*i + j
            plot_2d( ax, time, np.real(rho_plot[:,i,j]), xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[color_idx], name=r'$\rho_{%d%d}$' % (state_a,state_b))
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
         
      ax.legend(frameon=False)
      plt.tight_layout()
      plt.savefig( 'rho_transformed_real.png', bbox_inches='tight', dpi=print_dpi)
      plt.close()
   
      ylabel = r'Im[$\rho_{ij}$]'
      fig, ax = plt.subplots()
      fig.set_size_inches( fig_width_inches, fig_height_inches)
      fig.set_dpi(dpi)
      ax.tick_params(direction='in')
      ax.set_xlabel(xlabel, fontsize=fontsize)
      ax.set_ylabel(ylabel, fontsize=fontsize)

      for i in range(2):
         for j in range(2):
            state_a = list_states[i]
            state_b = list_states[j]
            color_idx = 2*i + j
            plot_2d( ax, time, np.imag(rho_plot[:,i,j]), xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[color_idx], name=r'$\rho_{%d%d}$' % (state_a,state_b))
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )     
         
      ax.legend(frameon=False)
      plt.tight_layout()
      plt.savefig( 'rho_transformed_imag.png', bbox_inches='tight', dpi=print_dpi)
      plt.close()
   
      # Change rho_plot to a 2D matrix for saving
      rho_plot_save = rho_plot.reshape(len(time),4)
      rho_save = np.zeros((len(time),5),dtype=np.cdouble)
      rho_save[:,0] = time
      rho_save[:,1:5] = rho_plot_save
      np.savetxt( 'RHO_transformed.dat', rho_save )
   
   else:

      for i in range(num_dim):
         
         ylabel = r'Re[$\rho_{ij}$]'
         fig, ax = plt.subplots()
         fig.set_size_inches( fig_width_inches, fig_height_inches)
         fig.set_dpi(dpi)
         ax.tick_params(direction='in')
         ax.set_xlabel(xlabel, fontsize=fontsize)
         ax.set_ylabel(ylabel, fontsize=fontsize)

         real_indices = [get_py_indx_rhofile( i+1, j+1, 'real', total_dim ) for j in range(num_dim)]
         imag_indices = [get_py_indx_rhofile( i+1, j+1, 'imag', total_dim ) for j in range(num_dim)]

         for j in range(num_dim):
            plot_2d( ax, time, rho[:,real_indices[j]], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[j], name=r'$\rho_{%d%d}$' % (i+1,j+1))
         if vlines:
            for line in vlines:
               ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
         
         ax.legend(frameon=False)
         plt.tight_layout()
         plt.savefig( 'rho_row%d_real.png' % (i), bbox_inches='tight', dpi=print_dpi)
         plt.close()
         
         ylabel = r'Im[$\rho_{ij}$]' 
         fig, ax = plt.subplots()
         fig.set_size_inches( fig_width_inches, fig_height_inches)
         fig.set_dpi(dpi)
         ax.tick_params(direction='in')
         ax.set_xlabel(xlabel, fontsize=fontsize)
         ax.set_ylabel(ylabel, fontsize=fontsize)
         
         for j in range(num_dim):
            plot_2d( ax, time, rho[:,imag_indices[j]], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[j], name=r'$\rho_{%d%d}$' % (i+1,j+1))
         if vlines:
            for line in vlines:
               ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 ) 
              
         ax.legend(frameon=False)
         plt.tight_layout()
         plt.savefig( 'rho_row%d_imag.png' % (i), bbox_inches='tight', dpi=print_dpi)
         plt.close()

def plot_pops( num_populations=3, datafiles=['POPULATIONS.dat'], savename='Populations.png', names=[r'$\downarrow$',r'$\uparrow$',r'$\emptyset$'], time_diff=None, vlines=None ):

   xlabel = r'Time (ns)'
   ylabel = r'Population'
   xfact = 1
   yfact = 1
   
   matplotlib.rc('xtick', labelsize=labelsize)
   matplotlib.rc('ytick', labelsize=labelsize)
   plt.rcParams['figure.figsize'] = (fig_width_inches, fig_height_inches)
   plt.rcParams['figure.dpi'] = dpi

   prop_cycle = plt.rcParams['axes.prop_cycle']
   colors = prop_cycle.by_key()['color']

   fig, ax = plt.subplots()
   fig.set_size_inches( fig_width_inches, fig_height_inches)
   fig.set_dpi(dpi)
   ax.tick_params(direction='in')
   ax.set_xlabel(xlabel, fontsize=fontsize)
   ax.set_ylabel(ylabel, fontsize=fontsize)
   
   esr_data = {}
   time = {}
   for i, file in enumerate(datafiles):
      esr_data[file] = get_data( filename=file )
      time[file] = esr_data[file][:,0]
      xmin = time[file][0] 
      xmax = time[file][-1]
      if time_diff:
         xdiff = time_diff
      else:
         xdiff = time[file][1] - time[file][0]
      for j in range(num_populations):
         plot_2d( ax, time[file], esr_data[file][:,1+j], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], marker=None, linewidth=1.5, color=colors[j], name=names[j] )
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )
   ax.legend(frameon=False)
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)
      
def check_final_states( filename, times):
   rho = get_data( filename=filename )
   arraysize = len(rho[0,:])
   print(len(rho[:,0]))
   diag_states = [(1,1),(2,2),(3,3),(4,4)]
   coh_states = [(1,4),(4,1)]
   other_states = []
   for istate in range(1,9):
      for jstate in range(1,9):
         state = (istate, jstate)
         if state not in diag_states:
            if state not in coh_states:
               other_states.append((istate,jstate))
            
   for time in times:
      print("Time (ns): ",time)
      itime = int(time / (1000./500000))-1
      diag_val = 0.0
      coh_val = 0.0
      other_val = 0.0
      for state in diag_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         diag_val = diag_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      for state in coh_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         coh_val = coh_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      for state in other_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         other_val = other_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      print('Diag vals:',diag_val)
      print('Coh vals:',coh_val)
      print('Other vals:',other_val)
      #if abs(other_val) >= 0.1:
      #print('Other States %d %d:' % (istate,jstate),val)
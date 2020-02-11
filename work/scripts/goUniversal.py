#!/usr/bin/env python

##########################################################################
##  goUniversal.py
##
##  A python script to run or submit jobs for the common use cases
##  of the IMSRG++ code. We check whether there is a pbs or slurm
##  scheduler, assign the relevant input parameters, set names
##  for the output files, and run or submit.
##  						-Ragnar Stroberg
##  						TRIUMF Nov 2016
######################################################################

from os import path,environ,mkdir,remove
from sys import argv
from subprocess import call,PIPE
from time import time,sleep
from datetime import datetime

### Check to see what type of batch submission system we're dealing with
BATCHSYS = 'NONE'
if call('type '+'qsub', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'PBS'
elif call('type '+'srun', shell=True, stdout=PIPE, stderr=PIPE) == 0: BATCHSYS = 'SLURM'

### The code uses OpenMP and benefits from up to at least 24 threads
NTHREADS=32
exe = '%s/bin/imsrg++'%(environ['HOME'])

### Flag to swith between submitting to the scheduler or running in the current shell
#batch_mode=False
batch_mode=True
if 'terminal' in argv[1:]: batch_mode=False

### This comes in handy if you want to loop over Z
ELEM = ['n','H','He','Li','Be','B','C','N',
       'O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K',
       'Ca','Sc','Ti','V','Cr','Mn','Fe','Co',  'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y',
       'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',  'Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
       'Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb']# ,'Bi','Po','At','Rn','Fr','Ra','Ac','Th','U','Np','Pu']

### ARGS is a (string => string) dictionary of input variables that are passed to the main program
ARGS  =  {}

### Maximum value of s, and maximum step size ds
ARGS['smax'] = '500'
ARGS['dsmax'] = '0.5'

#ARGS['lmax3'] = '10' # for comparing with Heiko

### Norm of Omega at which we split off and start a new transformation
ARGS['omega_norm_max'] = '0.25'

### Model space parameters used for reading Darmstadt-style interaction files
ARGS['file2e1max'] = '14 file2e2max=28 file2lmax=10'
ARGS['file3e1max'] = '14 file3e2max=28 file3e3max=14'

### Name of a directory to write Omega operators so they don't need to be stored in memory. If not given, they'll just be stored in memory.
#ARGS['scratch'] = 'SCRATCH'

### Generator for core decoupling, can be atan, white, imaginary-time.  (atan is default)
#ARGS['core_generator'] = 'imaginary-time'
### Generator for valence deoupling, can be shell-model, shell-model-atan, shell-model-npnh, shell-model-imaginary-time (shell-model-atan is default)
#ARGS['valence_generator'] = 'shell-model-imaginary-time'

### Solution method
ARGS['method'] = 'magnus'
#ARGS['method'] = 'brueckner'
#ARGS['method'] = 'flow'
#ARGS['method'] = 'HF'
#ARGS['method'] = 'MP3'

### Tolerance for ODE solver if using flow solution method
#ARGS['ode_tolerance'] = '1e-5'

if BATCHSYS == 'PBS':
  FILECONTENT = """#!/bin/bash
#PBS -N %s
#PBS -q oak
#PBS -d %s
#PBS -l walltime=192:00:00
#PBS -l nodes=1:ppn=%d
#PBS -l vmem=250gb
#PBS -o imsrg_log/%s.o
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=%d
%s
  """

elif BATCHSYS == 'SLURM':
  FILECONTENT = """#!/bin/bash
#SBATCH --account=rrg-holt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=%d
#SBATCH --mem=125G
#SBATCH --output=imsrg_log/%s.%%j
#SBATCH --time=%s
time srun %s
"""

### Make a directory for the log files, if it doesn't already exist
if not path.exists('imsrg_log'): mkdir('imsrg_log')

### Loop over multiple jobs to submit
for Z in range(10,11):
 A=Z
 for reference in ['%s%d'%(ELEM[Z],A)]:
  ARGS['reference'] = reference
  for e in [6]:
   for hw in [2]:

     ARGS['emax'] = '%d'%e
     ARGS['2bme'] = '/global/scratch/exch/AtomicHamiltonian/Coulomb_Hamil_LO%d_emax10_e2max20.snt'%(hw)
     ARGS["2bme"] = "/project/6006601/tmiyagi/MtxElmnt/AtomicHamiltonian/Coulomb_Hamil_LO%d_emax10_e2max20.snt"%(hw)
     ARGS['LECs'] = 'Coulomb'
     ARGS["physical_system"]="atomic"
     ARGS["basis_type"]="AO"
     ARGS["fmt2"]="tokyo"
     ARGS["hw"]=str(hw)
     ARGS["A"]=str(A)
     ARGS["valence_space"]="0hw-shell"

     ARGS['Operators'] = ''    # Operators to consistenly transform, separated by commas.

    ### Make an estimate of how much time to request. Only used for slurm at the moment.
     time_request = '0-1:00'
     if   e <  5 : time_request = '0-1:00'
     elif e <  8 : time_request = '0-1:00'
     elif e < 10 : time_request = '0-1:00'
     elif e < 12 : time_request = '0-1:00'

     jobname  = '%s_%s_%s_%s_e%s_s%s_hw%s_A%s' %(ARGS['valence_space'], ARGS['LECs'],ARGS['method'],ARGS['reference'],ARGS['emax'],ARGS['smax'],ARGS['hw'],ARGS['A'])
     logname = jobname + datetime.fromtimestamp(time()).strftime('_%y%m%d%H%M.log')

  ### Some optional parameters that we probably want in the output name if we're using them
     if 'lmax3' in ARGS:  jobname  += '_l%d'%(ARGS['lmax3'])
     if 'eta_criterion' in ARGS: jobname += '_eta%s'%(ARGS['eta_criterion'])
     if 'core_generator' in ARGS: jobname += '_' + ARGS['core_generator']
     if 'BetaCM' in ARGS: jobname += '_' + ARGS['BetaCM']
     ARGS['flowfile'] = 'output/BCH_' + jobname + '.dat'
     ARGS['intfile']  = 'output/' + jobname

     cmd = ' '.join([exe] + ['%s=%s'%(x,ARGS[x]) for x in ARGS])

  ### Submit the job if we're running in batch mode, otherwise just run in the current shell
     if batch_mode==True:
       sfile = open(jobname+'.batch','w')
       if BATCHSYS == 'PBS':
          sfile.write(FILECONTENT%(jobname,environ['PWD'],NTHREADS,logname,NTHREADS,cmd))
          sfile.close()
          call(['qsub', jobname+'.batch'])
       elif BATCHSYS == 'SLURM':
          sfile.write(FILECONTENT%(NTHREADS,jobname,time_request,cmd))
          sfile.close()
          call(['sbatch', jobname+'.batch'])
       remove(jobname+'.batch') # delete the file
       sleep(0.1)
     else:
       call(cmd.split())  # Run in the terminal, rather than submitting


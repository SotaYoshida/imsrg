#!/usr/local/bin/python3
import os, sys
HOME = os.path.expanduser('~')
path_ui = HOME+'/kshell/kshell-20180517/bin/'
import numpy as np
import subprocess
elist = [2,4,6,8]
hwlist = [10,15,20,25,30,35,40,45,50]
Nucl = 'O18'
for e in elist:
    for hw in hwlist:
        snt = 'magnus_O18_custom_NAT_hw' + str(hw) +\
                '_e' + str(e) + '_A18.snt'
        snt_base = os.path.splitext(os.path.basename(snt))[0]
        if(not os.path.isfile(snt)): continue
        f = open('input.txt','w')
        f.write('\n')  # single / MPI: default is single
        f.write(snt+'\n') # interaction file
        f.write('0,2\n')  # # of nucleons in valence space
        f.write('\n') # J, Parity, number of lowst states
        f.write('\n') # truncation for "+" parity
        f.write('\n') # truncation for "-" parity
        f.write('beta_cm = 0.0\n') # modifying parameters
        f.write('\n') # no transition probabilities
        f.write('\n') # number of valence protons and neutrons
        f.write('\n')
        f.close()
        cmd = 'python ' + path_ui + 'kshell_ui.py < input.txt'
        subprocess.call(cmd,shell=True)
        cmd = './'+str(Nucl)+'_'+snt_base+'.sh'
        subprocess.call(cmd,shell=True)

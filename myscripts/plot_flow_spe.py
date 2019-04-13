#!/usr/local/bin/python3
import os, sys
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt

llist = ["s","p","d","f","g","h","i","j","k","l","m","n","o"]

def get_label(i):
    # get label in imsrg code
    n,l,j,z = get_nljz(i)
    if(z == -1): pn = "p"
    if(z ==  1): pn = "n"
    return pn+str(n)+llist[l]+str(j)


def get_nljz(i):
    idx = 0
    for N in range(100):
        for l in range(N,-2,-2):
            if(l < 0): continue
            n = (N-l)//2
            for j2 in range(2*l+1,2*l-3,-2):
                if(j2 < 0): continue
                for tz in [-1,1]:
                    if(i == idx): return n, l, j2, tz
                    idx += 1

f = sys.argv[1]
fout = 'flow_spe_'+ os.path.splitext(os.path.basename(f))[0] + '.pdf'

plt.rcParams = myplt.set_style(roman=True)
fig,axs = myplt.set_canvas(1,1,12,16)
colors = myplt.colors()
proton_spe = {}
neutron_spe = {}

cnt = 0
is_spe = True
i = 1
while is_spe == True:
    i+=1
    if(i%2 == 0): continue
    try:
        X,Y = myplt.data2d(f,x=1,y=i)
        label = get_label(i-2)
        axs.plot(X[:],Y[:],label=label)
    except:
        is_spe = False

#axs.set_ylim(90,110)
### comment ###
axs.set_xlabel('$s$')
axs.set_ylabel('Energy (MeV)')
### comment ###
plt.legend(ncol=2)
plt.tight_layout()
#plt.show()
plt.savefig(fout)

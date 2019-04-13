#!/usr/local/bin/python3
import os, sys
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt

llist = ["s","p","d","f","g","h","i","j","k","l","m","n","o"]
emax = 6
norbs = (emax+1) * (emax+2)
def get_label_2b(idx):
    i,j = get_idx_2b(idx)
    return get_label(i) + " " + get_label(j) + " " + \
            get_label(i) + " " + get_label(j)

def get_idx_2b(idx):
    cnt = 0
    for i in range(norbs):
        for j in range(i+1):
            if(idx == cnt): return i, j
            cnt += 1

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
fout = 'flow_gamma_'+ os.path.splitext(os.path.basename(f))[0] + '.pdf'

plt.rcParams = myplt.set_style(roman=True)
fig,axs = myplt.set_canvas(1,1,12,16)
colors = myplt.colors()
ls = ["-","--",":","-."]
proton_spe = {}
neutron_spe = {}

cnt = 0
is_spe = True
label_ref = ["n0s1","n0p1","n0p3","n1s1","n0d3","n0d5", "n1p1", "n1p3", "n0f5", "n0f7"
        "p0s1","p0p1","p0p3","p1s1","p0d3","p0d5", "p1p1", "p1p3", "p0f5", "p0f7"]
label_ref = ["n0s1","n0p1","n0p3","n1s1","n0d3","n0d5", "p0s1","p0p1","p0p3","p1s1","p0d3","p0d5"]
#label_ref = ["n0s1","n0p1","n0p3", "p0s1","p0p1","p0p3"]
#label_ref = ["n0p1","n0p3","n1s1","n0d3","n0d5", "p0p1","p0p3","p1s1","p0d3","p0d5"]
#label_ref = ["n1s1","n0d3","n0d5", "p1s1","p0d3","p0d5"]
#label_ref = ["n0p1","n0p3","p0p1","p0p3"]
i = 1
while is_spe == True:
    i+=1
    try:
        c = colors[i%len(colors)]
        l = ls[i%len(ls)]
        label = get_label_2b(i-2)
        label_each = label.split()
        #print(label_each)
        if(not label_each[0] in label_ref): continue
        if(not label_each[1] in label_ref): continue
        if(not label_each[2] in label_ref): continue
        if(not label_each[3] in label_ref): continue
        X,Y = myplt.data2d(f,x=1,y=i)
        axs.plot(X[:],Y[:],c=c,ls=l,label=label)
        #axs.plot(X[:],Y[:],c=c,ls=l)
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

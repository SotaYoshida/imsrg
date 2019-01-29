#!/usr/local/bin/python3
import os, sys
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt

hwlist = range(10,55,5)
elist = range(2,10,2)
colors = myplt.colors()
marks = myplt.markers()
cnt = 0
for emax in elist:
    c = colors[cnt]
    m = marks[cnt]
    x = []
    y = []
    for hw in hwlist:
        f = HOME + '/Dropbox/Results_IMSRG_2018/data/EM500_srg2/Summary_magnus_O16_custom_NAT_hw' + str(hw) + \
                '_e' + str(emax) + '_A16.dat'
        _y = myplt.data_point(f,r=5, c=2)
        if(_y == None): continue
        x.append(hw)
        y.append(_y)
    plt.plot(x,y,c=c,marker=m)
    cnt += 1
plt.show()

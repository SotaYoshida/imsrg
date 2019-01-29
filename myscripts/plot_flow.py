#!/usr/local/bin/python3
import os, sys
import matplotlib.pyplot as plt
HOME=os.path.expanduser('~')
sys.path.append(HOME)
from python_myset.myplot import myplt

f = sys.argv[1]
X, Y = myplt.data2d(f,x=1,y=2)
plt.plot(X,Y)

f = sys.argv[2]
X, Y = myplt.data2d(f,x=1,y=2)
plt.plot(X,Y)
plt.show()

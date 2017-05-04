# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 12:47:36 2015

@author: wolfganggottwald
"""

import math
import matplotlib.pyplot as py
print "Bitte geben sie die Zahl ein, deren Wurzel sie ausrechnen wollen"
a=input()
reala=math.sqrt(a)
xn=1.0
xnn=0.5*(xn+(a/xn))
l=[]
k=0
while abs(xnn-reala)>0.00001:
   xnn=0.5*(xnn+(a/xnn))
   k=k+1
   print abs(xnn-reala),"Differenz"
   l.append(abs(xnn-reala))
print xnn,"Ergebnis der Iteration"
print reala,"tatsächlicher Wert"
print l
py.plot(l)
py.show(l)



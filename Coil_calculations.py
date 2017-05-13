# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as co
import time

import csv
import pandas as pd
import scipy.optimize 

    #Rcoil=inner diameter of coil [m]
    #Dwire=diameter of cable[m]
    #Nturns= windings per layer  (counter for that is n)
    #Nlayers=number of layers (counter for that is m)
    #L=length of coil =N_turns*D_wire [m]
    #I=current[A]
    #pr= r coordinate of the point, where the field is calculated
    #phi= angle coordinate of the point P, where field is calc.
    #pz= z coordinate of -||-
    
    #rm= R_coil+((m+0.5)*D_wire)....distance from z-axis to center 
    #of cable in m-th layer
    
    #theta=angle coordinate of point Q in m-th layer, goes from
    #0 to 2pi in 1000 steps
    
    #zn= z coordinate of nth winding
    #mu=magnetic field constant
    #dtheta=step distance at theta summation
   


thetasteps=1000
dtheta=2*math.pi/thetasteps

class coil:
    
    def __init__(self,Rcoil,Dwire,Nturns,Nlayers,I,z0):
        """
    Rcoil=inner diameter of coil [m]
    Dwire=diameter of cable[m]
    Nturns= windings per layer  (counter for that is n)
    Nlayers=number of layers (counter for that is m)
    I=current[A]
    z0=z coordinate of the middle of the coil    
        """
        self.Rcoil=Rcoil
        self.Dwire=Dwire
        self.Nturns=Nturns
        self.Nlayers=Nlayers
        self.I=I
        self.z0=z0



    
    def rm(self,m):
        return (self.Rcoil+m*self.Dwire-0.5*self.Dwire)   

    
    def zn(self,n):
        return (n*self.Dwire-0.5*self.Dwire-self.Nturns*0.5*self.Dwire) 
        
    
    
    def radialfield(self,r,z):
        """
        gives the radialfieldcomponent Br at a point (r,z) in gauß
        
        """
        B_r=0.
 
        for m in range(1,self.Nlayers+1,1):
            for n in range(1,self.Nturns+1,1):
                for i in np.arange(0,2*math.pi,dtheta):
                    B_r+=(self.rm(m)*(z-self.z0-self.zn(n))*dtheta*math.cos(i))/((r**2+((z-self.z0)-self.zn(n))**2+(self.rm(m))**2-2.*r*self.rm(m)*math.cos(i))**(1.5))
        return (self.I*co.mu_0/(4.*math.pi))*B_r*10000

    
    def longitudinalfield(self,r,z):
        """
        gives the longitudinalfieldcomponent Bz at a point (r,z) in gauß
        """
        B_z=0.
        phi=0.
        for m in range(1,self.Nlayers+1,1):
            for n in range(1,self.Nturns+1,1):
                for i in np.arange(0,2*math.pi,dtheta):
                    B_z+=(-self.rm(m)*(r*math.cos(i-phi)-self.rm(m))*dtheta)/((r**2+((z-self.z0)-self.zn(n))**2+(self.rm(m))**2-2.*r*self.rm(m)*math.cos(i))**(1.5))
        return ((self.I*co.mu_0)/(4.*math.pi))*B_z*10000

        
        
        
        
        
        
        
        
        
        #longfieldcalc_r und radialfieldcalc_r calculate the field for 
        #one r value and different z values between the origin
        # which is in the center of the coil and the set limit
        


        

        
         #longfieldcalc_z und radialfieldcalc_z calculate the field for one z value
               
    def longfieldcalc_z(self,r,dr,z):
        """
        calls longitudinalfield
        calculates Bz for a fixed z value and r values between -r to r in dr steps
        returns a list with the entries
        
        
        """
        l=[]
        for b in np.arange(-r,r+dr,dr):
            l.append(self.longitudinalfield(b,z))
        return l

        
        
    def radialfieldcalc_z(self,r,dr,z):
        """
        calls radialfield
        sets a value for z  and calculates Br at r values between -r to r with stepping distance dr
        returns list
        """
        l=[]
        for b in np.arange(-r,r+dr,dr):
                l.append(self.radialfield(b,z))
        return l

        
        
    def plot_longfield_z(self,r,dr,z):
        """
        calls longfieldcalc_z
        plots Bz at a set z value and r values between -r and r with stepping distance dr
        image can be saved, by changing plt.show() to plt.savefig()
        
        """
        xz=[]
        for i in np.arange(-r,r+dr,dr):
            xz.append(i)
        
        plt.plot(xz,self.longfieldcalc_z(r,dr,z))
        plt.ylabel("Bz [Gauss]")
        plt.xlabel("distance from center [m]")
        plt.show()
        #plt.savefig("Bz_fixed_z.eps")

        
        
    def plot_radialfield_z(self,r,dr,z):
        
        """
        calls radialfieldcalc_z
        plots Br at a set z value and r values between -r and r with stepping distance dr
        image can be saved, by changing plt.show() to plt.savefig()
        """
        xr=[]
        for i in np.arange(-r,r+dr,dr):
            xr.append(i)
        plt.plot(xr,self.radialfieldcalc_z(r,dr,z))
        plt.ylabel("Br [Gauss]")
        plt.xlabel("distance from center [m]")
        plt.show()
        #plt.savefig("Br_fixed_z.eps")

            
          
        
        
        
        
        

    def radialfieldcalc_r(self,r,zstart,zlimit,zstep):
        
        """
        calls radialfield 
        sets a value for r and calculates Br at z values from zstart to zlimit in stepping distance zsteps 
        returns list 
        """
        l=[]

        for a in np.arange(zstart,zlimit+zstep,zstep):
            l.append(self.radialfield(r,a))
        return l        
            
    def longfieldcalc_r(self,r,zstart,zlimit,zstep):
        
        """
        calls longitudinalfield
        sets a value for r and calculates Bz at z values from zstart to zlimit in stepping distance zsteps
        returns list
        """
        l=[]
        for b in np.arange(zstart,zlimit+zstep,zstep):
            l.append(self.longitudinalfield(r,b))
        return l
    
    def plot_longfield_r(self,r,zstart,zlimit,zstep):
        """
        calls longfieldcalc_r
        plots Bz at set r value with z values from zstart to zlimit with stepping distance zsteps
        """
        xz=[]

        for i in np.arange(zstart,zlimit+zstep,zstep):
            xz.append(i)
        
        plt.plot(xz,self.longfieldcalc_r(r,zstart,zlimit,zstep))
        plt.ylabel("Bz [Gauss]")
        plt.xlabel("distance from begin of coil[m]")
        plt.show()
        #plt.savefig("Bz_fixed_r.eps")
    
    def plot_radialfield_r(self,r,zstart,zlimit,zstep):
        
        """
        calls radialfieldcalc_r
        
        plots Br at set r value with z values from zstart to zlimit with stepping distance zsteps
        
        
        """
        
        xr=[]

        for i in np.arange(zstart,zlimit+zstep,zstep):
            xr.append(i)
        
        plt.plot(xr,self.radialfieldcalc_r(r,zstart,zlimit,zstep))
        plt.ylabel("Br [Gauss]")
        plt.xlabel("distance from begin of coil[m]")
        plt.show()
        #plt.savefig("Br_fixed_r.eps")
    
            


            
        

         #now plot field fo different values of z
    def zfielddifferentvalues(self,r,dr,z1,z2,z3,z4):
        
        """
        sets a z value from 4 different, then calculates and plots Bz from -r to r there
        then calls the next z value
        plots all 4 in one graph and saves it
    
        """
        t=[]
        for i in np.arange(-r,r+dr,dr):
            t.append(i)
        
        #Length of the coil and windings
        N=float(self.Nturns)*float(self.Nlayers) 
        L=self.Dwire*float(self.Nturns)
          #make strings where the legend description is saved
        stringa=str("z=")+str(z1)+str("m")
        stringb=str("z=")+str(z2)+str("m")
        stringc=str("z=")+str(z3)+str("m")
        stringd=str("z=")+str(z4)+str("m")
        stringLegende=str("L=")+str(L)+str("m")+str(",R=")+str(self.Rcoil)+str("m")+str(",N=")+str(N)
        #calc field at z values
        a=self.longfieldcalc_z(r,dr,z1)
        b=self.longfieldcalc_z(r,dr,z2)
        c=self.longfieldcalc_z(r,dr,z3)
        d=self.longfieldcalc_z(r,dr,z4)
          #makes a string with the legend description
        stringLegende=str("L=")+str(L)+str("m")+str(",R=")+str(self.Rcoil)+str("m")

        #value at the end gives the z coordinate where the field is calculated, e.g. 0,04m
        plt.plot(t,a,'k-',label=stringa)
        plt.plot(t,b,'g-',label=stringb)        
        plt.plot(t,c,'b-',label=stringc)
        plt.plot(t,d,'r-',label=stringd)
        plt.legend(loc="best",title=stringLegende, ncol=1, shadow=True, fancybox=True)
        plt.ylabel("Bz [Gauss]")
        plt.xlabel("distance from center [mm]")
        #plt.show()
        Dateiname=str("Bz_4values_")+stringLegende
        plt.savefig(Dateiname+".eps")
        #saves plot
        
    def rfielddifferentvalues(self,r,dr,z1,z2,z3,z4):
        """
        sets a z value from 4 different, then calculates and plots Br from -r to r there
        then calls the next z value
        plots all 4 in one graph and saves it
    
        """
        t=[]
        for i in np.arange(-r,r+dr,dr):
            t.append(i)

        L=self.Dwire*float(self.Nturns)
        N=float(self.Nturns)*float(self.Nlayers)       
         #make strings where the description of the legend is saved
        stringa=str("z=")+str(z1)+str("m")
        stringb=str("z=")+str(z2)+str("m")
        stringc=str("z=")+str(z3)+str("m")
        stringd=str("z=")+str(z4)+str("m")
        stringLegende=str("L=")+str(L)+str("m")+str(",R=")+str(self.Rcoil)+str("m")+str(",N=")+str(N)
        #calculates field at z values
        a=self.radialfieldcalc_z(r,dr,z1)
        b=self.radialfieldcalc_z(r,dr,z2)
        c=self.radialfieldcalc_z(r,dr,z3)
        d=self.radialfieldcalc_z(r,dr,z4)
        plt.plot(t,a,'k-',label=stringa) 
        plt.plot(t,b,'g.',label=stringb)        
        plt.plot(t,c,'b.',label=stringc)
        plt.plot(t,d,'r.',label=stringd)
        plt.legend(loc="best",title=stringLegende, ncol=1, shadow=True, fancybox=True)
        plt.ylabel("Br [Gauss]")
        plt.xlabel("distance from center [mm]")
        Dateiname=str("Br_4values_")+stringLegende
        #plt.show()
        plt.savefig(Dateiname+".eps")

        
        

    def totalfield(self,r,z):
        """
        calls functions to calculate Br and Bz
        calculates totalfield from them
        returns value
        """
        return math.sqrt((self.longitudinalfield(r,z))**2+(self.radialfield(r,z))**2)
        
        
 
    def totalfieldcalc(self,r,dr,z):
        """
        saves values of totalfield for set z value in a list
        """
        tot=[]
        for i in np.arange(-r,r+dr,dr):
            tot.append(math.sqrt(abs(self.longitudinalfield(i,z))**2+abs(self.radialfield(i,z))**2))
        return tot
    
 


    #total field at one z value
    def plottotalfield(self,r,dr,z):
        """
        plots totalfield for set z value and r values between -r and r, stepping distance dr
        """
        xtot=[]
        for i in np.arange(-r,r+dr,dr):
            xtot.append(i)
        plt.plot(xtot,self.totalfieldcalc(r,dr,z))
        plt.ylabel("B [Gauss]")
        plt.xlabel("distance from center [mm]")
        plt.show()
        #plt.savefig("B.eps")
       
        
    def totalfielddifferentvalues(self,r,dr,z1,z2,z3,z4):
        """
        plots totalfield for 4 different z values and r varying from -r to r each time
        
        """
        t=[]
        for i in np.arange(-r,r+dr,dr):
            t.append(i)

        L=self.Dwire*float(self.Nturns)
        N=float(self.Nturns)*float(self.Nlayers)
        #make strings where the description of the legend is saved
        stringa=str("z=")+str(z1)+str("m")
        stringb=str("z=")+str(z2)+str("m")
        stringc=str("z=")+str(z3)+str("m")
        stringd=str("z=")+str(z4)+str("m")
        stringLegende=str("L=")+str(L)+str("m")+str(",R=")+str(self.Rcoil)+str("m")+str(",N=")+str(N)
        #calculates field at z values
        a=self.totalfieldcalc(r,dr,z1)
        b=self.totalfieldcalc(r,dr,z2)
        c=self.totalfieldcalc(r,dr,z3)
        d=self.totalfieldcalc(r,dr,z4)
        plt.plot(t,a,'k-',label=stringa) 
        plt.plot(t,b,'g.',label=stringb)        
        plt.plot(t,c,'b.',label=stringc)
        plt.plot(t,d,'r.',label=stringd)
        plt.legend(loc="best",title=stringLegende, ncol=1, shadow=True, fancybox=True)
        plt.ylabel("B_total [Gauss]")
        plt.xlabel("distance from center [mm]")
        Dateiname=str("B_tot_4values_")+stringLegende
        #plt.show()
        plt.savefig(Dateiname+".eps")

     
    
    def plot_OFS(self):
        """plots the optimal field shape for B0 in Gauß"""
        B0=self.tot_field_one_point(0.,0.)
        L=self.Dwire*self.Nturns
        z=np.arange(-L/2,L/2,0.005)
        plt.plot(z,B0*(np.cos(math.pi*z/L))**2,label='OFS')
        plt.legend(loc='best')
        plt.xlabel("z [m]")
        plt.ylabel("Optimal field shape [G]")
        plt.show()
    
        
        
    def approx_eq(self):
        """
        just a approximationequation for the magnetic field in the center r=0
        """
        L=self.Dwire*float(self.Nturns)
        N=float(self.Nturns)*float(self.Nlayers)
        return 10000*mu*N*self.I/L
        
    def Inductivity_Wheeler(self):
        """calculates the inductivity of one coil, where the radius is bigger  than the length, applied for the guidance coils"""
        Dwindings=self.Nlayers*self.Dwire
        Coillength=self.Nturns*self.Dwire
        N=self.Nturns*self.Nlayers
        return (0.8*self.Rcoil**2*N**2/(6*self.Rcoil+9*Coillength+10*Dwindings))*39.3701/10**6
    #conversion from mycroHenry into Henry, the 39. factor comes from the conversion into meters, as the formula was found in inches
    #found the formula as Wheelers formula on: www.home.earthlink.net/~jimlux/hv/wheeler.htm

    def Inductivity_longcoil(self):
        """calculates the inductivity for one long coil (l/R=0.6), ergo longer than the radius"""
        #formula from wikipedia
        return (self.Nturns*self.Nlayers)**2*co.mu_0*2*math.pi*self.Rcoil**2/(self.Dwire*self.Nturns+2*self.Dwire*self.Nlayers/(2.2))
    
    def Inductiv_resistance(self,frequency):
        """calculates the inductiv resistance of the coil"""
        return 2*math.pi*frequency*self.Inductivity_Wheeler()
    
    
        
class coilarray:

    def __init__(self,coils):
        self.coils=coils
        """
        reads a list of coils and performs operation with this coilarray
        """
        
    def r_field_one_point(self,r,z):
        """
        calls coil.radialfield for each coil in coilarray to calculate Br at (r,z) and sums up the all influences
        """
        B_r=0
        for coil in self.coils:
            B_r+=coil.radialfield(r,z)
            
        return B_r
            

        
        
    def r_field_calc_fixed_r(self,r,zstart,zlimit,zstep):
        """
        calls coilarray.r_field_one_point to calculate Br from all coils
        
        sets and r value and then goes from zstart to zlimit and
        calculates the field at points zstep apart from each other
        """
        B_r=[]
        for i in np.arange(zstart,zlimit+zstep,zstep):
            B_r.append(self.r_field_one_point(r,i))
        return B_r


    def r_field_plot_fixed_r(self,r,zstart,zlimit,zstep):
        """
        calls coilarray. r_field_calc_fixed_r
        plots Br for set r and z values from zstart to zlimit
        
        """
        l=[]
        for i in np.arange(zstart,zlimit+zstep,zstep):
            l.append(i)
            
        plt.plot(l,self.r_field_calc_fixed_r(r,zstart,zlimit,zstep),label="B_r-field of coils")
        plt.ylabel("$B_{rad}$ [Gauss]")
        plt.xlabel("z [m]")
        plt.legend(loc="best")
        #plt.show()
        plt.savefig("computed_Br_field_helmholtzcoils.png")
        
        
        

    def z_field_one_point(self,r,z):
        """
        calls coil.longitudinalfield to calculate Bz from all coils
        returns value
        """
        B_z=0
        for coil in self.coils:
            B_z+=coil.longitudinalfield(r,z)
            
        return B_z
        
        

    def z_field_calc_fixed_r(self,r,zstart,zlimit,zstep):
        """
        calls coilarray.z_field_one_point 
        sets an r value and calculates Bz at points from zstart to zlimit, which are zsteps apart
        """
        B_z=[]
        for i in np.arange(zstart,zlimit+zstep,zstep):
            B_z.append(self.z_field_one_point(r,i))
        return B_z
        
        
    def z_field_plot_fixed_r(self,r,zstart,zlimit,zstep):
        """
        calls coilarray.z_field_calc_fixed_r 
        plots Bz from zstart to zlimit for set r value
        """
        l=[]

        for i in np.arange(zstart,zlimit+zstep,zstep):
            l.append(i)
            
        plt.plot(l,self.z_field_calc_fixed_r(r,zstart,zlimit,zstep),label="Bz-field of  coils")
        plt.ylabel("$B_{long}$ [Gauss]")
        plt.xlabel("z [m]")
        plt.legend(loc="best")
        #plt.show()
        plt.savefig("computed_Bz_field_helmholtzcoils.png")


            

        
    #calculates the list of values for B for fixed z and all coils     
    def r_field_calc_fixed_z(self,z,r,dr):
        """
        calls coilarray.r_field_one_point
        sets z value
        calulates Br from -r to r with steps that are dr apart
        returns Br as a list
        """
        B_r=[]
        for i in np.arange(-r,r+dr,dr):
            B_r.append(self.r_field_one_point(i,z))
        return B_r


    def r_field_plot_fixed_z(self,z,r,dr):
        """
        plots Br for set z value and r between -r to r, with dr apart from each other
        """
        l=[]
        for i in np.arange(-r,r+dr,dr):
            l.append(i)
            
        plt.plot(l,self.r_field_calc_fixed_z(z,r,dr))
        plt.ylabel("Br [Gauss]")
        plt.xlabel("distance from origin [mm]")
        plt.show()
        
        
        


        
        
#plots the z field for one r value
    def z_field_calc_fixed_z(self,z,r,dr):
        """
        sets z value
        calls coilarray.z_field_one_point, calculates Bz at r coordinates from -r to r with dr apart
        returns Bz as list
       
        """
        B_z=[]
        for i in np.arange(-r,r+dr,dr):
            B_z.append(self.z_field_one_point(i,z))
        return B_z
        
        

    def z_field_plot_fixed_z(self,z,r,dr):
        """
        sets z value
        calls coilarray.z_field_calc_fixed_z
        plots Bz for fixed z value
        """
        l=[]

        for i in np.arange(-r,r+dr,dr):
            l.append(i)
            
        plt.plot(l,self.z_field_calc_fixed_z(z,r,dr))
        plt.ylabel("Bz [Gauss]")
        plt.xlabel("distance from origin [mm]")
        plt.show()
        


    #calculate total field 
    
    def tot_field_one_point(self,r,z):
        """
        calculates total field by calling coilarray.z_field_one_point and coilarray.r_field_one_point
        """
        return math.sqrt(abs((self.z_field_one_point(r,z)))**2+abs((self.r_field_one_point(r,z))**2))

    
    #now want to calculate integral of B over the length of a coilarray
    
    def J(self,r,zstart,zlimit,zstep):
        """
        calculates integral of B over the length of the coil:
            1. calls coilarray.tot_field_one_point to calculate tot_field at given point
            2.multiplies it with the stepping distance zstep
            3.weights it over the number of steps taken
            4.adds it to list
            5.sums list to value
        returns value for J
        """
        J=[]
        for i in np.arange(zstart,zlimit+zstep,zstep):
            J.append((self.tot_field_one_point(r,i)*zstep)/float(((zstart-zlimit)/zstep)))
        return np.sum(J)
    
    
    def dJ_over_J(self,r,zstart,zlimit,zstep):
        """
        normalizes J by dividing it by J(at center ,z=0)
        subtracts 1, so that dj/j is 0 at the center and multiplies with 1000, for the scale to be better visible
        because changes are from the order of 0.001 
        """
        return 1000*(self.J(r,zstart,zlimit,zstep)/self.J(0.,zstart,zlimit,zstep)-1)
 
        
    def plot_dJ_over_J(self,rads,zstart,zlimit,zstep,name):
        """
        plots dJ_J, needs a string of a name, in order for the graphics to
        be saved.
        """
        l=[]
        for rad in rads:
            l.append(self.dJ_over_J(rad,zstart,zlimit,zstep))


        plt.plot(rads,l,label=name)
        plt.legend()
        plt.ylabel(r'$1000*\left(\frac{\Delta J}{J}-1\right)$')
        plt.xlabel("distance from center [m]")
    
        #plt.show()
        plt.savefig("dJ_J_"+name+".eps")
        
            
#now some functions to write the data into csv files, in order to save it


    def B_fixed_r_to_csv(self,r,zstart,zlimit,zstep,coilname):
        """
        writes B for fixed r into rows of a csv file
        """
        
        Br=self.r_field_calc_fixed_r(r,zstart,zlimit,zstep)
        Bz=self.z_field_calc_fixed_r(r,zstart,zlimit,zstep)
        z_coordinate=[]
        for j in np.arange(zstart,zlimit+zstep,zstep):
            z_coordinate.append(j)
            
        table={'z_coordinate':z_coordinate,'Br':Br,'Bz':Bz}
        
        data=pd.DataFrame(table,columns=['z_coordinate','Br','Bz'],dtype=float)
        data.to_csv("B_fixed_r_"+str(coilname)+"_r="+str(r)+"m.csv")
    
    
    def B_fixed_z_to_csv(self,z,r,dr,coilname):
        """
        writes B for fixed z into rows of a csv file
        """
        
        Br=self.r_field_calc_fixed_z(z,r,dr)
        Bz=self.z_field_calc_fixed_z(z,r,dr)
        r_coordinate=[]
        for j in np.arange(-r,r+dr,dr):
            r_coordinate.append(j)
            
        table={'r_coordinate':r_coordinate,'Br':Br,'Bz':Bz}
        
        data=pd.DataFrame(table,columns=['r_coordinate','Br','Bz'],dtype=float)
        data.to_csv("B_fixed_z_"+str(coilname)+"_z="+str(z)+"m.csv")
        
#now we want to read those csv files and create plots, which saves us time, if we use the same coils all the time

    def B_fixed_r_read_plot_csv(self,coilname,OFS,coillength):
        """reads a csv file and plots it
        Enter csvname/picname as strings (with "")
        also plots the optimal field shape for that coil if wanted, enter Yes/No at OFS field for that"""

        if OFS=="Yes":
            data=pd.read_csv("B_fixed_r_"+coilname+".csv")
        
            
            middle=int(len(data['z_coordinate'])/2.)
            end=int(len(data['z_coordinate']))
            #setting the length of NSE coils to 0.724 m cause it has 724 windings of 0.001 m Dwire
            L=coillength
            B0=float(math.sqrt((data['Br'][middle])**2+(data['Bz'][middle])**2))

            xaxis=np.arange(-L/2,L/2,0.005)
            
                      
            plt.plot(data['z_coordinate'],data['Br'],label='$B_r$')
            plt.plot(data['z_coordinate'],data['Bz'],label='$B_z$')
            plt.plot(xaxis,B0*(np.cos(math.pi*xaxis/L))**2,label='OFS')
            plt.legend(loc='best')
            plt.xlabel("z coordinate[m]")
            plt.ylabel("magnetic field [gauss]")
        
            #plt.show()
            plt.savefig(coilname+".eps")
        elif OFS=="No":           
            data=pd.read_csv("B_fixed_r_"+coilname+".csv")
            plt.plot(data['z_coordinate'],data['Br'],label='$B_r$')
            plt.plot(data['z_coordinate'],data['Bz'],label='$B_z$')
            plt.legend(loc='best')
            plt.xlabel("z coordinate[m]")
            plt.ylabel("magnetic field [gauss]")

            #plt.show()
            plt.savefig(coilname+".eps")
        else:
            print "Try again, enter Yes for Optimal field shape or No for only the field"
        
    def B_fixed_z_read_plot_csv(self,coilname):
        """reads a csv file and plots it"""
        data=pd.read_csv("B_fixed_z_"+coilname+".csv")
        plt.plot(data['r_coordinate'],data['Br'],label='$B_r$')
        plt.plot(data['r_coordinate'],data['Bz'],label='$B_z$')
        plt.legend(loc='best')
        plt.xlabel("r coordinate[m]")
        plt.ylabel("magnetic field [gauss]")
        #plt.show()
        plt.savefig(coilname+".eps")
        
        
    def plot_OFS(self,length):
        """plots the optimal field shape for B0 in Gauß"""
        B0=self.tot_field_one_point(0.,0.)
        print B0, "B0"
        L=length
        z=np.arange(-L/2,L/2,0.005)
        plt.plot(z,B0*(np.cos(math.pi*z/L))**2,label='OFS')
        plt.legend(loc='best')
        plt.xlabel("z [m]")
        plt.ylabel("Optimal field shape [G]")
        plt.show()
   


   
def csv_read_dJ_J(rads,zstart,zlimit,zstep,csvname,plotname):
        
    """reads the B field of reseda from the csv file and then plots dJ over J
    rads is a list of radia, for which the caluclations can be done and zstart, zlimit and zstep refer 
    to the area where the field is calculated        
    csvname and plotname need to be strings, so put in in the form of "xxxx"
        """

    J=[]
    for r in rads:
        #ready computed RESEDA-files are named 
        #the program now reads several csv files ,which have been created in advance
        data=pd.read_csv(csvname+"_r="+str(r)+"m.csv")
        tot_field=[]
        for i in range(len(data)):
            tot_field.append(math.sqrt((data['Br'][i])**2+(data['Bz'][i])**2))
        J_r=[]
        for i in range(len(tot_field)):
            J_r.append(float((tot_field[i]*zstep)/float(((zstart-zlimit)/zstep))))
        J_r=np.sum(J_r)
        J.append(J_r)
        
    dJ_J=[]
    for i in range(len(rads)):
          dJ_J.append(1000000*((J[i]/J[0])-1))

    #now plotting it

    plt.plot(rads,dJ_J,label=plotname)
    plt.legend(loc='best')
    plt.ylabel(r'$\frac{\Delta J}{J}-1$[ppm]')
    plt.xlabel("distance from center [m]")
    
    #plt.show()
    plt.savefig("dJ_J_"+plotname+".eps")


def difference_between_fields(zstart,zlimit,zstep,csvname1,csvname2):
    data1=pd.read_csv(csvname1)
    data2=pd.read_csv(csvname2)
    r_diff=[]
    z_diff=[]
    t=[]
    for i in np.arange(len(data1)):
        r_diff.append(data1['Br'][i]-data2['Br'][i])
    for i in np.arange(len(data1)):
        z_diff.append(data1['Bz'][i]-data2['Bz'][i])
    for i in np.arange(zstart,zlimit+zstep,zstep):
        t.append(i)

    plt.plot(t,z_diff,label="Bz-difference")
    plt.plot(t,r_diff,label="Br-difference")
    plt.legend(loc="best")
    plt.xlabel("z coordinate, along coil axis[m]")
    plt.ylabel("difference in magnetic fields [G]")
    
    plt.show()
    
#following function is changeable, in order to generate plots, either only Br or only Bz, with OFS shape or not
#ideally one would write a list of radia, for which the plots should be generated, such as:
#rads=[0.0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02]
#for rad in rads:
#    B_fixed_r_read_plot_two_csv_files(rad,"B_fixed_r_NSE_arm1_old_r="+str(rad)+"m.csv","old_NSE","B_fixed_r_NSE_arm1_new_r="+str(rad)+"m.csv","new_NSE")
#    plt.clf()

def B_fixed_r_read_plot_two_csv_files(rad,csvname1,coilname1,csvname2,coilname2):
    data1=pd.read_csv(csvname1)
    data2=pd.read_csv(csvname2)
    #plt.plot(data1['z_coordinate'],data1['Br'],label="$B_r$"+str(coilname1))
    plt.plot(data1['z_coordinate'],data1['Bz'],label="$B_z$,"+str(coilname1))
    #plt.plot(data2['z_coordinate'],data2['Br'],label="$B_r$"+str(coilname2))
    plt.plot(data2['z_coordinate'],data2['Bz'],label="$B_z$,"+str(coilname2))
    #now the OFS shape
    middle1=int(len(data1['z_coordinate'])/2.)
    middle2=int(len(data2['z_coordinate'])/2.)
    end=int(len(data1['z_coordinate']))
    L=1.
    B0_old=float(math.sqrt((data1['Br'][middle1])**2+(data1['Bz'][middle1])**2))
    B0_v2=float(math.sqrt((data2['Br'][middle2])**2+(data2['Bz'][middle2])**2))
    print B0_v2
    xaxis=np.arange(-L/2,L/2,0.005)
    plt.plot(xaxis,B0_old*(np.cos(math.pi*xaxis/L))**2,label='OFS old')
    plt.plot(xaxis,B0_v2*(np.cos(math.pi*xaxis/L))**2,label='OFS v2')
    plt.legend(loc="best")
    plt.xlabel("z coordinate [m]")
    plt.ylabel("magnetic field [G]")
    #plt.show()
    plt.savefig("B_fields_r="+str(rad)+"m_with_OFS.eps")
    
def Bz_fixed_r_read_plot(coilname,rads):
    """enter name of coilarray, e.g. "helmholtzcoilarray", for which you wrote the magnetic field into csv files
    enter list of radia, for which you calculated the magnetic field
    plots the longitudinal field for fixed r values"""
    for rad in rads:
        filename=str("B_fixed_z_"+coilname+"_r="+str(rad)+"m.csv")
        data=pd.read_csv(filename)
        plt.plot(data['z_coordinate'],data['Bz'],label="r="+str(rad)+"m"+str(coilname))
    plt.legend(loc='best')
    plt.xlabel("distance from beamcenter [m]")
    plt.ylabel("$B_z$")
    #plt.show()
    plt.savefig("Bz_fixed_r_"+coilname+"r="+str(rads)+".eps")


def Br_fixed_r_read_plot(coilname,rads):
    """enter name of coilarray, e.g. "helmholtzcoilarray", for which you wrote the magnetic field into csv files
    enter list of radia, for which you calculated the magnetic field
    plots the radial field for fixed r values
        """
    for rad in rads:
        filename=str("B_fixed_r_"+coilname+"_r="+str(rad)+"m.csv")
        data=pd.read_csv(filename)
        plt.plot(data['z_coordinate'],data['Br'],label="r="+str(rad)+"m"+str(coilname))
    plt.legend(loc='best')
    plt.xlabel("z coordinate [m]")
    plt.ylabel("$B_r$[G]")
    #plt.show()
    plt.savefig("Br_fixed_r_"+coilname+"r="+str(rads)+".eps")
    
def Bz_fixed_z_read_plot(coilname,zvalues):
    """enter name of coilarray, e.g. "helmholtzcoilarray", for which you wrote the magnetic field into csv files
    enter list of radia, for which you calculated the magnetic field
    plots the longitudinal field for fixed z values
    """
    for zvalue in zvalues:
        filename=str("B_fixed_z_"+coilname+"_z="+str(zvalue)+"m.csv")
        data=pd.read_csv(filename)
        plt.plot(data['r_coordinate'],data['Bz'],label="z="+str(zvalue)+"m"+str(coilname))
    plt.legend(loc='best')
    plt.xlabel("distance from beamcenter [m]")
    plt.ylabel("$B_z$[G]")
    plt.savefig("Bz_fixed_z_"+coilname+"_z="+str(zvalues)+".eps")

def Br_fixed_z_read_plot(coilname,zvalues):
    """enter name of coilarray, e.g. "helmholtzcoilarray", for which you wrote the magnetic field into csv files
    enter list of radia, for which you calculated the magnetic field
    plots the radial field for fixed z values
    """
    for zvalue in zvalues:
        filename=str("B_fixed_z_"+coilname+"_z="+str(zvalue)+"m.csv")
        data=pd.read_csv(filename)
        plt.plot(data['r_coordinate'],data['Br'],label="z="+str(zvalue)+"m"+str(coilname))
    plt.legend(loc='best')
    plt.xlabel("distance from beamcenter [m]")
    plt.ylabel("$B_r$[G]")
    plt.savefig("Br_fixed_z_"+coilname+"_z="+str(zvalues)+".eps")

    

    
    #now for building OFS coils

def build_OFS_coils_csv(coilname,N_coils,d_insulator,d_wire,R_inner_coil,L_inner_coil,dist):
    """constructs a coilsetup, whose geometric form ist like the OFS, with an arbitrary amplitude A, which 
    corrresponds to how high the setup should be biult, when constructed
    N_coils is the number of coils 
    d_insulator could be everything between 0.001[m] and 0.02[m]
    d_wire is usually 0.0014 [m]
    R_innercoil is usually 0.05 [m], cause the neutron beam is 4 cm in diameter
    L_inner_coil can be everything between 1m and 1,4 m which is ideal if, 1,4 mm cable is used, cause then 
    we have exactly 1000 windings in the innermost coil
    dist corresponds to how much difference there should be between the hypothetical maximum of the
    cos^2 shape of the coilarray and the last coil constructed (usually some centimeters)
    we don't want the coils to be too short, because then the radial field gets more inhomogenuos
    function returns dwire
    saves the input as a txt file
    """
    radius=lambda n: R_inner_coil+n*(d_wire+d_insulator)
    #the nth coil gets radius(n) as its radius
    A=radius(N_coils)+dist-R_inner_coil
    #A is the maximum of the
    OFS=lambda x: A*(np.cos(np.pi*x/L_inner_coil)**2)+R_inner_coil
    x_coordinate=lambda x: np.arccos(np.sqrt((radius(x)-R_inner_coil)/A))*L_inner_coil/np.pi
    l=[]
    r=[]
    t=[]
    cl=[]
    
    for n in np.arange(0,N_coils+1,1):
        l.append(round(x_coordinate(n),4))
        r.append(radius(n))
        cl.append(round(2*x_coordinate(n),4))
        #rounding the length to millimeters
        
    for n in np.arange(0,N_coils+1,1):
        t.append(int(round(cl[n]/d_wire,0)))

                 
    #take into account here, that we round the number of turns to the next float with .0, the conversion to integers happens, when
    #you call build_OFS_coilarray!!
    #had some problems with computing this...
    table={'radius':r,'turns':t,"coillength":cl}
    data=pd.DataFrame(table,columns=['radius','turns','coillength'],dtype=float)

    data.to_csv("coilparameters_for_"+str(coilname)+".csv")
    #now plot the result
    #ax=np.arange(-L_inner_coil/2,L_inner_coil/2,0.001)
    #plt.xlabel("z coordinate alongside coil [m]")
    #plt.ylabel("OFS field shape,with coils as hlines")
    #plt.plot(ax,OFS(ax))
    #plt.scatter(l,r)
    #for i in np.arange(0,N_coils+1,1):
    #    plt.hlines(OFS(l[i]),-l[i],l[i])
        
    #plt.savefig(str(coilname)+"_coil_distribution.eps")
    #print data
    #want to save the starting parameters as a txt file, so that this info doesnt get lost in the process
    #starting_parameters=open("starting_parameters_for_"+str(coilname)+".txt","w")
    #starting_parameters.write("coilname="+str(coilname)+" ------ Number of coils="+str(N_coils)+" ------ Thickness of insulator[m]="+str(d_insulator)
                             #+" ------ diameter of wire[m]="+str(d_wire)+" ------ radius inner coil[m]= "+str(R_inner_coil)+" ------ length inner coil[m]="
                             #+str(L_inner_coil)+" ------ distance between top winding of Nth coil and max of hypothetical cos^2 shape[m]="+str(dist))

def build_OFS_coilarray(coilname,current,z0,D_wire,rads,zstart,zlimit,zsteps):
    """
    enter coilname as str
    reads csv file from build_OFS_coils_csv
    makes OFS coilarray from it
    calculates field and makes csv file from it which is then saved, just like in B_fixed_r_to_csv
    saves the input as a txt file
    """
    N_layers=1
    csvname="coilparameters_for_"+coilname+".csv"
    #saving the field_calc_parameters as a textfile, so that the info is not lost
    #field_params=open("field_calc_parameters_for_"+str(coilname)+".txt","w")
    #field_params.write("current[A]="+str(current)+" ------ z0 [m]="+str(z0)+" ------ radia where calculated [m]="+
                         #str(rads)+" ------ z boundaries where field is calculated [m]="+str(zstart)+","+str(zlimit)+" ------ z stepping range[m]="+str(zsteps))

    
    data=pd.read_csv(csvname)
    built_coilarray=[]
    
    for i in np.arange(len(data['radius'])):
        built_coilarray.append(coil(data['radius'][i],D_wire,int(data['turns'][i]),N_layers,current,z0))
    for rad in rads:
        coilarray(built_coilarray).B_fixed_r_to_csv(rad,zstart,zlimit,zsteps,coilname)
        



    

start_time = time.time()
rads=[0.0,0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02]
 
#build_OFS_coils_csv("NSE_v5",15,0.001,0.0014,0.05,1.4,0.01)
#build_OFS_coilarray("NSE_v5",5.,0.,0.0014,rads,-0.6,0.6,0.05)
build_OFS_coils_csv("NSE_v9",12,0.001,0.0014,0.12,1.4,0.01)
build_OFS_coilarray("NSE_v9",2.5,0.,0.0014,rads,-0.6,0.6,0.05)
end_time=time.time()
calc_time=open("calculation_time_NSE_v9.txt",'w')
calc_time.write(str((end_time-start_time)/60.)+" minutes")
#new array
start_time2=time.time()
build_OFS_coils_csv("NSE_v10",15,0.001,0.0014,0.05,1.4,0.005)
build_OFS_coilarray("NSE_v10",2.,0.,0.0014,rads,-0.6,0.6,0.05)
#NSE_v3.B_fixed_r_read_plot_csv("NSE_v3_r=0.0m","Yes",1.)

end_time2=time.time()
calc_time2=open("calculation_time_NSE_v10.txt",'w')
calc_time2.write(str((end_time2-start_time2)/60.)+" minutes")


       





from __future__ import print_function

import numpy
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

plt.close()


# load the image
image = Image.open('C:/Users\sury_\Downloads\Atto\Task11_holey_mirror.png')
# convert image to numpy array
data = np.array(image)

#Calibration abd center
xpix=data.shape[1]
ypix=data.shape[0]
x=xpix*7/169.5
y=ypix*7/169.5
xx=np.linspace(0,x,xpix)
yy=np.linspace(0,y,ypix)
cx=1095*7/169.5
cy=887*7/169.5


#Plot data
plt.figure(1)
plt.imshow(data, cmap='jet',extent=[0,x,y,0])
plt.colorbar(label = 'Intensity, a.u.')
plt.ylabel('y[mm]')
plt.xlabel('x[mm]')
plt.axvline(x = cx, color = 'k', linestyle = '--')
plt.axhline(y =  cy, color = 'k', linestyle = '--')


#Selecting row an colum of the center of mass and max.
intx=data[887,:]
inty=data[:,1095]

intx = intx.astype(float)
intx[1007:1185,]=np.nan
np.savetxt('C:/Users\sury_\Downloads\Atto\intx.txt',intx)


a1 =       204.94
b1 =       45.17
c1 =       6.26
d =       15.04
fitx=a1*np.exp(-2*((xx-b1)/c1)**2)+d

plt.figure(2)
plt.plot(xx,intx,'*b')
plt.plot(xx,fitx)
plt.xlabel('x [mm]')
plt.ylabel('Intensity[a.u]')
plt.legend(['Data','Fit'])


#Plot for y

inty = inty.astype(float)
inty[798:977,]=np.nan
np.savetxt('C:/Users\sury_\Downloads\Atto\inty.txt',inty)


a1 =       129.81
b1 =       36.387
c1 =        6.26
d =       13.74
fity = a1*np.exp(-2*((yy-b1)/c1)**2)+d

plt.figure(3)
plt.plot(yy,inty,'*b')
plt.plot(yy,fity)
plt.legend(['Data','Fit'])
plt.xlabel('y [mm]')
plt.ylabel('Intensity [a.u]')


#%%
lamda=1030*10**-9 #m
winx=6.386*10**-3
winy=5.34*10**-3
en=1*10**-3 #J
tau=40*10**-15 #fs
f=90 #cm
wx=lamda*f/(np.pi*winx)#[cm]
Ipeakx=(2*en/(tau*wx**2*np.pi))*np.sqrt((4*np.log(2))/np.pi)#W/cm²
wy=lamda*f/(np.pi*winy)#[cm]
Ipeaky=(2*en/(tau*wy**2*np.pi))*np.sqrt((4*np.log(2))/np.pi)#W/cm²
Eholex=0.48
Eholey=0.60
Ipeakxh=Ipeakx-(Ipeakx*Eholex)
Ipeakyh=Ipeaky-(Ipeaky*Eholey)
#%%
lum=1.03 #lambda um
Ip=15.76 #eV
Ein=1.2398/lum
Up=9.33*10**-14*Ipeakx*lum**2 #eV
Ecoff=Ip+3.2*Up #eV
coffx=Ecoff/Ein

Up=9.33*10**-14*Ipeaky*lum**2 #eV
Ecoff=Ip+3.2*Up #eV
coffy=Ecoff/Ein

Up=9.33*10**-14*Ipeakxh*lum**2 #eV
Ecoff=Ip+3.2*Up #eV
coffxh=Ecoff/Ein

Up=9.33*10**-14*Ipeakyh*lum**2 #eV
Ecoff=Ip+3.2*Up #eV
coffyh=Ecoff/Ein




imagenf = Image.open('C:/Users\sury_\Downloads\Atto\Task13_nofilter.png')
# convert image to numpy array

datanf = np.array(imagenf)
plt.figure(4)
plt.title('With no filter')
plt.imshow(datanf, cmap='jet')
plt.colorbar(label = 'Intensity, a.u.')
plt.ylabel('y [pix]')
plt.xlabel('x [pix]')

wlc=np.sum(datanf, axis=0) #wavelength components
sp=np.sum(datanf,axis=1) #spatial profile

plt.figure(5)
plt.plot(wlc)
plt.ylabel('Intensity[a.u]')
plt.xlabel('x [pix]')

plt.figure(6)
plt.plot(np.linspace(0,len(datanf),len(datanf)),sp)
np.savetxt('C:/Users\sury_\Downloads\Atto\sp.txt',sp)
plt.xlabel('y [pix]')
plt.ylabel('Intensity[a.u]')


#For Mg


imageMg = Image.open('C:/Users\sury_\Downloads\Atto\Task13_Mgfilter.png')
dataMg = np.array(imageMg)

plt.figure(7)
plt.title('Mg with filter')
plt.imshow(dataMg, cmap='jet')
plt.colorbar(label = 'Intensity, a.u.')
plt.ylabel('y [pix]')
plt.xlabel('x [pix]')
wlcMg=np.sum(dataMg, axis=0) #wavelength components
spMg=np.sum(dataMg,axis=1) #spatial profile

plt.figure(8)
plt.plot(wlcMg)
plt.ylabel('Intensity [a.u]')
plt.xlabel('x [pix]')

plt.figure(9)
plt.plot(np.linspace(0,len(dataMg),len(dataMg)),spMg)
plt.xlabel('y [pix]')
plt.ylabel('Intensity[a.u]')

#Last image

a1 =   4.9e+05
b1 =       16.5
c1 =       5.1
d =   1.74e+04
mm=np.linspace(0,40,len(sp))
fitsm=a1*np.exp(-2*((mm-b1)/c1)**2)+d

plt.figure(13)
plt.plot(mm,sp/max(sp))
plt.plot(mm,fitsm/max(fitsm),'-*')
plt.xlabel('x [mm]')
plt.ylabel('Intensity [a.u]')
plt.legend(['Data','Fitting'])
plt.xlim(2,25)
plt.show()


#For Te

imageTe = Image.open('C:/Users\sury_\Downloads\Atto\Task13_Tefilter.png')
dataTe = np.array(imageTe)

plt.figure(10)

plt.title('Te Filter')
plt.imshow(dataTe, cmap='jet')
plt.colorbar(label = 'Intensity, a.u.')
plt.ylabel('y [pix]')
plt.xlabel('x [pix]')
wlcTe=np.sum(dataTe, axis=0) #wavelength components
spTe=np.sum(dataTe,axis=1) #spatial profile

plt.figure(11)
plt.plot(wlcTe)
plt.ylabel('Intensity[a.u]')
plt.xlabel('x [pix]')

plt.figure(12)
plt.plot(np.linspace(0,len(dataTe),len(dataTe)),spTe)
plt.xlabel('y [pix]')
plt.ylabel('Intensity [a.u]')
























#read no filter image

imagenf = Image.open('Task13_nofilter.png')
# convert image to numpy array
datanf = np.array(imagenf)
x2_sum = np.sum(datanf,axis=0)
y2_sum= np.sum(datanf,axis=1)



#Te image
Tefilter=Image.open('Task13_Tefilter.png')
Data4=np.array(Tefilter)
x4_sum=np.sum(Data4,axis=0)
y4_sum=np.sum(Data4,axis=1)

#Mg image

Mgfilter=Image.open('Task13_Mgfilter.png')
Data3=np.array(Mgfilter)
x3_sum=np.sum(Data3,axis=0)
y3_sum=np.sum(Data3,axis=1)

#Calibration

Mg= (x3_sum)/ (x2_sum) # dividing with filter/ without filter
Te= (x4_sum)/ (x2_sum)



#Plot

plt.figure()
plt.plot(Mg)
plt.plot(Te)
plt.legend(['Mg-ref', 'Te-ref'])
plt.ylabel('Transmision')
plt.xlabel('Pixels')




#Fitting parameters

a= 4.592e+04   #fitting parameters
b=  -11.78
unnamed = np.arange(1,1937)

xaxis= (a/unnamed)+b  #Converting pixels into eV

Tefilter_eV= x4_sum
Mgfilter_eV= x3_sum
nofilter_eV= x2_sum

Mg_norm= ((Mgfilter_eV)/np.max(Mgfilter_eV))
no_norm= ((nofilter_eV)/np.max(nofilter_eV))
Te_norm= ((Tefilter_eV)/np.max(Tefilter_eV))

#Real data
trans_Mg= np.loadtxt('Mg1.txt')
y2=Mg_norm/no_norm





# create figure and axis objects with subplots()
fig,ax = plt.subplots()
# make a plot
ax.plot(xaxis,Mg_norm,color="red",label='Mg filter')
ax.plot(xaxis,no_norm,color = 'orange',label = 'No filter')
# set x-axis label
ax.set_xlabel('Photon Energy (eV)')
# set y-axis label
ax.set_ylabel('Intensity')
ax.legend()
plt.xlim(15,80)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(trans_Mg[:,0],trans_Mg[:,1],color ='green',label = 'Theoretical T')
ax2.plot(xaxis,y2,color = 'blue', label = 'Experimental T')
ax2.set_ylabel('Transmission')
ax2.legend()



#Te
a1=  5.07e+04
b1=-23.57
xaxis1= a1/unnamed+b1
trans_Te=np.loadtxt('TeData1.txt')
y2=Te_norm/no_norm


# create figure and axis objects with subplots()
fig1,ax = plt.subplots()
# make a plot
ax.plot(xaxis,Te_norm,color="red",label = 'Te filter')
ax.plot(xaxis,no_norm,color = 'orange',label = 'No filter')
# set x-axis label
ax.set_xlabel('Photon Energy (eV)')
# set y-axis label
ax.set_ylabel('Intensity')
ax.legend()
plt.xlim(10,80)

# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(trans_Te[:,0],trans_Te[:,1],color = 'green',label = 'Theoretical T')
ax2.plot(xaxis,y2,color = 'blue',label = 'Experimental T')
ax2.set_ylabel('Transmission')
ax2.legend()

plt.show()



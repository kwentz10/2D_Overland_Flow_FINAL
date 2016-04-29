# -*- coding: utf-8 -*-
"""
Created on Wed Mar  13 2006

@author: Katherine
"""

from __future__ import division
import numpy as np
from landlab import RasterModelGrid
from matplotlib import pyplot as plt
from landlab.plot.imshow import imshow_node_grid


#__________INITIALIZE___________#
#Initial Parameters
zmax=30 #max elevation
valley_s=0.005 #valley slope
side_s=0.02 #side slope

#Space Array
x_min=0 #m
x_max=100 #m
dx=10 #m
spacearray=np.arange(x_min, x_max+dx, dx)
num_rows=len(spacearray)
num_cols=len(spacearray)

#Time Array
t_min=0 #s
t_max=18000 #s
dt=0.3 #s
timearray=np.arange(t_min, t_max+dt, dt)
time_interval=[timearray[i:i+int((round(len(timearray)/4)))] for i in range(0,len(timearray),int(round(len(timearray)/4)))]

#2D grid dimensions
x_ncells=num_cols
y_ncells=num_rows
mg=RasterModelGrid(x_ncells,y_ncells,dx) #make 2D grid

#Create data fields
zb=mg.add_empty('node','land_surface_elevation')
h=mg.add_zeros('node','water_thickness')
zw=mg.add_empty('node','water_elevation')
w_slope=mg.add_zeros('link','water_slope') #make all links have zero slope
q=mg.add_zeros('link','flux')

#Initialize elevation
zb[:]=zmax-valley_s*mg.node_x
zb+=side_s * np.abs(mg.node_y - mg.dx * ((num_rows - 1) / 2.0))
zw[:]=zb+h

#Define constants 
n=0.045 #roughness of bed surface 
r=mg.add_ones('node','rain_rate')
inf=mg.add_ones('node','infiltration_rate')

#Boundary conditions (ENWS)
mg.set_closed_boundaries_at_grid_edges(False, True, True, True) #water can flow out of right side

#Plot setup
nplots=10 #number of plots
tplot=int(t_max/nplots)

#Number of runs
imax=len(timearray)


#__________RUN___________#
for t in range(imax):
    
    #Parameters for Manning Eq
    w_slope[mg.active_links]=mg.calculate_gradients_at_active_links(zw)  
#    h_edge=mg.map_value_at_max_node_to_link(mg,'water_elevation','water_thickness',out=None)
    h_edge=mg.map_max_of_link_nodes_to_link('water_thickness')

    #Calculate ice flux
    q=-np.sign(w_slope)*(1/n)*(h_edge**(5/3))*((abs(w_slope))**(1/2))
    dqdx=mg.calculate_flux_divergence_at_nodes(q[mg.active_links])
    
    #Change in rain rate over time
    if timearray[t] in time_interval[1]:
        r[:]=r*(10/1000/3600)   
    elif timearray[t] in time_interval[2]:
        r[:]=r*(10/1000/3600)*3
    elif timearray[t] in time_interval[3]:
        r[:]=r*(10/1000/3600)*6
    elif timearray[t] in time_interval[4]:
        r[:]=r*(10/1000/3600)*9 
    
    #Change in infiltration rate over time
    if timearray[t] in time_interval[1]:
        inf[:]=inf*(5/1000/3600)  
    elif timearray[t] in time_interval[2]:
        inf[:]=inf*(5/1000/3600)*20
    elif timearray[t] in time_interval[3]:
        inf[:]=inf*(5/1000/3600)*10
    elif timearray[t] in time_interval[4]:
        inf[:]=inf*(5/1000/3600)
    
    #Calculate h and zw
    dhdt=-dqdx+r-inf
    h[mg.core_nodes]+=dhdt[mg.core_nodes]*dt
    for i in range(0, len(h)): #make sure bottom limit of z does not go below zb
        if h[i]<0:
            h[i]=0
    zw[:]=h+zb

    #Plot figures 
    if timearray[t]%tplot==0:
        
        water_thickness_raster=mg.node_vector_to_raster(h)
        water_elevation_raster=mg.node_vector_to_raster(zw)
        water_thickness_raster_in_mm=water_thickness_raster*1000
        
        fb1=plt.figure(1) #figure blueprint
        figure1 = plt.subplot()
        figure1.set_title('Water Thickness')
        figure1.set_xlabel('Distance (m)')
        figure1.set_ylabel('Distance (m)')
        l1=figure1.imshow(water_thickness_raster_in_mm) #add in ticks=v for stable min/max of colorbar
        u1=plt.colorbar(mappable=l1)
        u1.set_label('Millimeters')
        plt.pause(0.00001)
        fb1.clear()
        
        fb2=plt.figure(2)
        figure2=plt.subplot()       
        figure2.clear()
        figure2.set_title('Water Elevation')
        figure2.set_xlabel('Distance (m)')
        figure2.set_ylabel('Distance (m)')
        l2=figure2.imshow(water_elevation_raster)
        u2=plt.colorbar(mappable=l2)
        u2.set_label('Meters')
        plt.pause(0.00001)
        fb2.clear()
        
   
#__________FINALIZE___________# 
water_thickness_raster=mg.node_vector_to_raster(h)
water_elevation_raster=mg.node_vector_to_raster(zw)
water_thickness_raster_in_mm=water_thickness_raster*1000

fb1=plt.figure(1) #figure blueprint
figure1 = plt.subplot()
figure1.set_title('Water Thickness')
figure1.set_xlabel('Distance (m)')
figure1.set_ylabel('Distance (m)')
l1=figure1.imshow(water_thickness_raster_in_mm)
u1=plt.colorbar(mappable=l1) #add in ticks=v for stable min/max of colorbar
u1.set_label('Millimeters')
        
fb2=plt.figure(2)
figure2=plt.subplot()       
figure2.clear()
figure2.set_title('Water Elevation')
figure2.set_xlabel('Distance (m)')
figure2.set_ylabel('Distance (m)')
l2=figure2.imshow(water_elevation_raster)
u2=plt.colorbar(mappable=l2)
u2.set_label('Meters')

          
    

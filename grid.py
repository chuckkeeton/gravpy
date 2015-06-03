#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay 
import random

import trinterior as trint

b = 1 # radius of magnification function

def f(x):
    '''our magnification function'''
    return np.sqrt(b**2-x**2)
    
def nf(x):
    '''our magnification function (negative y-axis values)'''
    return -1* np.sqrt(b**2-x**2)
    
def carmapping(x,y):
    '''mapping of cartesian coordinates from image to source plane'''
    x = np.array(x)
    y = np.array(y)
    return np.transpose([x - b*x/np.sqrt(x**2 + y**2), y - b*y/np.sqrt(x**2 + y **2)])

def polmapping(r,th):
    '''mapping of polar coordinates from image to source plane'''
    r  = np.array(r)
    th = np.array(th)
    
    return np.transpose( [np.cos(th)(r-b), np.sin(th)(r-b)] )

def polartocar(r,th):
    '''convert polar coordinates to cartesian coordinates'''
    r  = np.array(r)
    th = np.array(th)
            
    return np.transpose([r*np.cos(th),r*np.sin(th)])

def condition(x,y,coord):
    '''returns the distance between the critical curve and the point '''
    if coord == 'car':
        return np.sqrt(x**2 + y**2) - b
    elif coord == 'pol':
        return x - b
    else:
        return NaN

def relation(x,y,coord):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    dif = condition(x,y,coord)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

def buildrelations(xran,yran,coord):
    '''applies relation() on 2D range specified by xran and yran. Returns a 2D array.'''
    xx,yy = np.meshgrid(xran,yran,sparse=True)
    
    return np.transpose(relation(xx,yy,coord))
    
def changep(fir,sec,thr,frt):
    '''[[fir, sec],[thr,frt]] - returns true if there is a change in magnification in the box specified'''
    return (fir*sec < 0 or sec*frt < 0 or frt*thr < 0 or thr*fir < 0)

RECURSEDEPTH = 4 #just here for future purpose when we want to set the max depth

def points(stack,xran,yran,coord,n=0):
    '''collects the points around our magnification function and appends them to the 'stack'. Subdivides grid by 2 in each recursion.'''
    mat = buildrelations(xran,yran,coord)
    xlen = len(xran)
    ylen = len(yran)
    for i in range(xlen-1):
        for j in range(ylen-1):
            if n != RECURSEDEPTH and changep(mat[i][j],mat[i][j+1],mat[i+1][j],mat[i+1][j+1]):
                points(stack,
                       np.linspace(xran[i],xran[i+1],3,endpoint=True),
                       np.linspace(yran[j],yran[j+1],3,endpoint=True),
                       coord,n+1)
            else:
                stack.append([xran[i+1],yran[j]])
                stack.append([xran[i+1],yran[j+1]])
                stack.append([xran[i],yran[j]])
                stack.append([xran[i],yran[j+1]])

##### SCRIPT BEGINS ######
    
# plot our function, 200 points
funx = np.linspace(-1,1,200,endpoint=False)
funy1 = f(funx)
funy2 = nf(funx) 

## critical curves
critx = np.hstack((funx,funx[::-1])) 
crity = np.hstack((funy1,funy2[::-1])) #here and above the second pairs are put in reverse order so the line connecting the points makes a nice circe rather than crossing the origin to start from the other side

critpairs = zip(critx,crity)

## caustics?

transcrit = carmapping(critx,crity)

# initial cartesian grid, coarse,
spacing = 0.5
upperend = 2.5
lowerend = -upperend
x = np.arange(lowerend,upperend+spacing,spacing)
y1 = np.arange(0,upperend+spacing,spacing)
y2 = np.arange(lowerend,spacing,spacing)

#intial polar grid
rspacing = 0.03
rupper = 1.2
thetadivisions = 42
r = np.arange(0,rupper+rspacing,rspacing)
theta = np.linspace(0,2*np.pi,thetadivisions)

## the rs and thetas formed from a cartesian product. Used for vector operations.
rs = np.tile(r,len(theta))
thetas = np.repeat(theta,len(r))

#polstack = [] # stack for holding the polar points #needed if we subgrid on polar grids
carstack = [] # stack for holding the cartesian points

# SUBGRIDDING!!!!
#points(polstack,r,theta,'pol') # if we wanted to subgrid on the polar grid, uncomment
points(carstack,x,y1,'car') #generate subgrid on cartesian grid on for the +y axis range
points(carstack,x,y2,'car') #generate subgrid on cartesian grid on for the -y axis range

carstack = np.array(carstack) 
polstack = np.array(polartocar(rs,thetas))
stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
transformed = np.array(carmapping(*np.transpose(stack))) #transform pairs from image to source plane

dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

# finding source positions
image = (random.random()*2-1,random.random()*2-1) #a random coordinate pair with x,y chosen from -1:1

lenstri = transformed[dpoints.simplices.copy()] #triangles on the source plane
imagetri= stack[dpoints.simplices.copy()] #triangles on the image plane

indices = trint.find(image,lenstri) #list of which triangles contain point on source plane
sourcetri = imagetri[indices] 
sourcepos = np.sum(sourcetri,axis=1)/3.0 #list of the centroid coordinates for the triangles which contain the point 'image'

# plot them all

plt.subplot(1,2,1) # image plane
plt.title('Image Plane')
plt.axis([lowerend,upperend,lowerend,upperend])
plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

plt.plot(critx,crity,zorder=2) # plot of critical curve(s)

plt.triplot(stack[:,0],stack[:,1], dpoints.simplices.copy(),zorder=1) # plot of the Delaunay Triangulization

plt.scatter(*zip(*sourcepos), marker='*', color='black',s=100, zorder=2) # plot of the (approximate/centroid) positions of the image

plt.subplot(1,2,2) # source plane
plt.title('Source Plane')
plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

plt.triplot(transformed[:,0],transformed[:,1], dpoints.simplices.copy(), zorder=1) # plot of the transformed Delaunay Triangulization

plt.scatter(*zip(image), marker='*', color='red', s= 200, zorder=2 ) # plot of (observed) image position

plt.plot(transcrit[:,0],transcrit[:,1], color ='green', zorder=2) # plot of caustics?
plt.show()


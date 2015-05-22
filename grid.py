#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay 


import trinterior as trint

import random

b = 1 # radius of magnification function

def f(x):
    '''our magnification function'''
    return np.sqrt(b**2-x**2)
    
def nf(x):
    '''our magnification function (negative y-axis values)'''
    return -1* np.sqrt(b**2-x**2)
    
def condition(x,y,coord):
    '''returns the distance between the critical curve and the point '''
    if coord == 'car':
        return np.sqrt(x**2 + y**2) - b
    elif coord == 'pol':
        return x - b
    else:
        return NaN

def carmapping(stack):
    newstack = []
    for i in stack:
        (x,y) = i
        newstack.append([x - b*x/np.sqrt(x**2 + y**2), y - b*y/np.sqrt(x**2 + y **2)])
    return newstack
        
def polmapping(stack):
    newstack = []
    for i in stack:
        (r,th) = i
        newstack.append( [np.cos(th)(r-b),np.sin(th)(r-b)] )
    return newstack

def polartocar(stack):
    newstack = []
    for i in stack:
        (r,th) = i
        newstack.append([r*np.cos(th),r*np.sin(th)])
    return newstack

def relation(x,y,coord):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    dif = condition(x,y,coord)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

#    if dif < 0: 
#        return -1
#    elif dif > 0:
#        return 1
#    else: return 0


def buildrelations(xran,yran,coord):
    '''applies relation() on 2D range specified by xran and yran. Returns a 2D array.'''
    return [[relation(j,i,coord) for i in yran] for j in xran]
    
def changep(fir,sec,thr,frt):
    '''[[fir, sec],[thr,frt]] - returns true if there is a change in magnification in the box specified'''
    return (fir*sec < 0 or sec*frt < 0 or frt*thr < 0 or thr*fir < 0)

RECURSEDEPTH = 4 #just here for future purpose when we want to set the max depth

def points(stack,xran,yran,coord,n=0):
    '''collects the points around our magnification function and appends them to the 'stack'. Subdivides grid by 2 in each recursion.'''
    mat = buildrelations(xran,yran,coord)
    xindex = len(xran)
    yindex = len(yran)
    for i in range(xindex-1):
        for j in range(yindex-1):
            if changep(mat[i][j],mat[i][j+1],mat[i+1][j],mat[i+1][j+1]):
                if n == RECURSEDEPTH:
                    stack.append([xran[i+1],yran[j]])
                    stack.append([xran[i+1],yran[j+1]])
                    stack.append([xran[i],yran[j]])
                    stack.append([xran[i],yran[j+1]])
                else: 
                    points(stack,
                           np.linspace(xran[i],xran[i+1],3,endpoint=True),
                           np.linspace(yran[j],yran[j+1],3,endpoint=True),
                           coord,n+1)
            else: #if n == 0: 
                stack.append([xran[i+1],yran[j]])
                stack.append([xran[i+1],yran[j+1]])
                stack.append([xran[i],yran[j]])
                stack.append([xran[i],yran[j+1]])

                    

    
# plot our function, 200 points
funx = np.linspace(-1,1,200,endpoint=False)
funy1 = [f(i) for i in funx]
funy2 = [nf(i) for i in funx]

## critical curves
critx = np.hstack((funx,funx[::-1])) 
crity = np.hstack((funy1,funy2[::-1])) #here and above the second pairs are put in reverse order so the line connecting the points makes a nice circe rather than crossing the origin to start from the other side

critpairs = zip(critx,crity)

## caustics?

transcrit = np.array(carmapping(critpairs))

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


polstack = [] # stack for holding the polar points
[[polstack.append((ra,th)) for th in theta] for ra in r]
carstack = [] # stack for holding the cartesian points

# SUBGRIDDING!!!!
#points(polstack,r,theta,'pol')
points(carstack,x,y1,'car') #generate subgrid on cartesian grid on for the +y axis range
points(carstack,x,y2,'car') #generate subgrid on cartesian grid on for the -y axis range

carstack = np.array(carstack) 
polstack = np.array(polartocar(polstack))
stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
transformed = np.array(carmapping(stack)) #transform pairs from image to source plane

dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

# finding source positions
image = (random.random()*2-1,random.random()*2-1) #a random coordinate pair with x,y chosen from -1:1

lenstri = transformed[dpoints.simplices.copy()] #triangles on the source plane
imagetri= stack[dpoints.simplices.copy()] #triangles on the image plane

indices = trint.find(image,lenstri) #list of which triangles contain point on source plane
sourcetri = imagetri[indices] 
sourcepos = [np.sum(t,axis=0)/3.0 for t in sourcetri] #list of the centroid coordinates for the triangles which contain the point 'image'

# plot them all

plt.subplot(1,2,1) # image plane
plt.title('Image Plane')
plt.axis([lowerend,upperend,lowerend,upperend])
plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

plt.plot(critx,crity,zorder=2) # plot of critical curve(s)

#plt.plot(*np.transpose(stack),marker=',',color='m',ls='') # plot of the points on grid
plt.triplot(stack[:,0],stack[:,1], dpoints.simplices.copy(),zorder=1) # plot of the Delaunay Triangulization

plt.scatter(*zip(*sourcepos), marker='*', color='black',s=100, zorder=2) # plot of the (approximate) positions of the image

plt.subplot(1,2,2) # source plane
plt.title('Source Plane')
plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

plt.triplot(transformed[:,0],transformed[:,1], dpoints.simplices.copy(), zorder=1) # plot of the transformed Delaunay Triangulization

plt.scatter(*zip(image), marker='*', color='red', s= 200, zorder=2 ) # plot of (observed) image position

plt.plot(transcrit[:,0],transcrit[:,1], color ='green', zorder=2) # plot of caustics?
plt.show()


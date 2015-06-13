#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op
import random

import trinterior as trint

def f(x,b):
    '''our magnification function'''
    return np.sqrt(b**2-x**2)
    
def nf(x,b):
    '''our magnification function (negative y-axis values)'''
    return -1* np.sqrt(b**2-x**2)
    
def mapping(v,u,b):
    w,z = u #image
    x,y = v #actual position we want to find
    return [x - b*x/np.sqrt(x**2 + y**2) - w, y - b*y/np.sqrt(x**2 + y **2) - z]

def carmapping(x,y,b):
    '''mapping of cartesian coordinates from image to source plane'''
    x = np.array(x)
    y = np.array(y)
    return np.transpose([x - b*x/np.sqrt(x**2 + y**2), y - b*y/np.sqrt(x**2 + y **2)])

def polmapping(r,th,b):
    '''mapping of polar coordinates from image to source plane'''
    r  = np.array(r)
    th = np.array(th)
    
    return np.transpose( [np.cos(th)(r-b), np.sin(th)(r-b)] )

def polartocar(r,th):
    '''convert polar coordinates to cartesian coordinates'''
    r  = np.array(r)
    th = np.array(th)
            
    return np.transpose([r*np.cos(th),r*np.sin(th)])

def distance(x,y,coord,b):
    '''returns the distance between the critical curve and the point '''
    if coord == 'car':
        return np.sqrt(x**2 + y**2) - b
    elif coord == 'pol':
        return x - b
    else:
        raise NameError('Unrecogonized Coordinate System')

def relation(x,y,coord,b):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    dif = distance(x,y,coord,b)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

def buildrelations(xran,yran,coord,b):
    '''applies relation() on 2D range specified by xran and yran. Returns a 2D array.'''
    xx,yy = np.meshgrid(xran,yran,sparse=True)
    return np.transpose(relation(xx,yy,coord,b))
    
def mag_change(fir,sec,thr,frt):
    '''[[fir, sec],[thr,frt]] - returns true if there is a change in magnification in the box specified'''
    return (fir*sec < 0 or sec*frt < 0 or frt*thr < 0 or thr*fir < 0)

RECURSEDEPTH = 4 #just here for future purpose when we want to set the max depth

def points(stack,xran,yran,coord,b,n=0):
    '''collects the points around our magnification function and appends them to the 'stack'. Subdivides grid by 2 in each recursion.'''
    mat = buildrelations(xran,yran,coord,b)
    xlen = len(xran)
    ylen = len(yran)
    for i in range(xlen-1):
        for j in range(ylen-1):
            if n != RECURSEDEPTH and mag_change(mat[i][j],mat[i][j+1],mat[i+1][j],mat[i+1][j+1]):
                points(stack,
                       np.linspace(xran[i],xran[i+1],3,endpoint=True),
                       np.linspace(yran[j],yran[j+1],3,endpoint=True),
                       coord,b,n+1)
            else:
                stack.append([xran[i+1],yran[j]])
                stack.append([xran[i+1],yran[j+1]])
                stack.append([xran[i],yran[j]])
                stack.append([xran[i],yran[j+1]])

def neg_mag(x,y):
    r = np.sqrt(x**2+y**2)
    mag = 1/(1-1/r)
    return mag <= 0

    
def mag_limit(x,y,b):
    r = np.sqrt(x**2+y**2)
    mag = r/(r-b)
    return mag > 3 or mag < -2
    
def sprinkle_points(stack,xran,yran,num):
    xmin, xmax = xran
    ymin, ymax = yran
    for i in range(num):
        stack.append(
            (random.uniform(xmin,xmax), random.uniform(ymin,ymax)))


    
                         
def points2(stack,xran,yran,spacing,b):
    xs = np.tile(xran,len(yran))
    ys = np.repeat(yran,len(xran))
    pairs = zip(xs,ys)
    sides,neg_sides = 3,7
    width,neg_width = spacing/sides,spacing/neg_sides
    
    for i in pairs:
        x,y = i
        #if neg_mag(x,y):
         #   sprinkle_points(stack,(x-spacing/2.,x+spacing/2.),(y-spacing/2.,y+spacing/2.),200)
        if mag_limit(x,y,b):
            sprinkle_points(stack,(x-spacing/2.,x+spacing/2.),(y-spacing/2.,y+spacing/2.),100)
            
    stack += pairs
        
        
def generate_ranges(critargs,carargs,polargs):
    # generate our critical curve function, 200 points
    b, numpoints = critargs
    funx = np.linspace(-b,b,numpoints,endpoint=False)
    funy1 = f(funx,b)
    funy2 = nf(funx,b) 

    ## critical curves
    critx = np.hstack((funx,funx[::-1])) 
    crity = np.hstack((funy1,funy2[::-1])) #here and above the second pairs are put in reverse order so the line connecting the points makes a nice circe rather than crossing the origin to start from the other side

    ## caustics?

    caustics = carmapping(critx,crity,b)

    # initial cartesian grid, coarse,
    lowerend, upperend, spacing = carargs
    #spacing = 0.5
    #upperend = 2.5
    #lowerend = -upperend
    x = np.arange(lowerend,upperend+spacing,spacing)
    y = np.arange(lowerend,upperend+spacing,spacing)

    #intial polar grid
    
    rupper, rspacing, thetadivision = polargs
    #rspacing = 0.03
    #rupper = 0.8 #rspacing*10
    #thetadivision = 42
    thetaspacing = 2*np.pi/thetadivision
    r = np.arange(0,rupper+rspacing,rspacing)
    theta = np.arange(0,2*np.pi,thetaspacing)

    ## the rs and thetas formed from a cartesian product. Used for vector operations.
    rs = np.tile(r,len(theta))
    thetas = np.repeat(theta,len(r))

    return [ [critx, crity], caustics, [x,y], [r, theta], [rs, thetas]]



def transformations(car_ranges, pol_ranges, spacing, b):

    x,y = car_ranges
    r,theta = pol_ranges
    
    polstack = [] # stack for holding the polar points #needed if we subgrid on polar grids
    carstack = [] # stack for holding the cartesian points
    stack = []
    #points(polstack,r,theta,'pol') # if we wanted to subgrid on the polar grid, uncomment
    #points(carstack,x,y,'car',b) #generate subgrid on cartesian grid
    points2(carstack,x,y,spacing,b)
    stack = np.array(stack)
    carstack = np.array(carstack) 
    polstack = np.array(polartocar(r,theta)) #comment if we subgrid on polar grids
    stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
    transformed = np.array(carmapping(stack[:,0],stack[:,1],b)) #transform pairs from image to source plane

    dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

    return [stack, transformed, dpoints]



def find_source(stack, transformed, simplices, image_loc,b):
    lenstri = transformed[simplices] #triangles on the source plane
    imagetri= stack[simplices] #triangles on the image plane

    indices = trint.find2(image_loc,lenstri.copy()) #list of which triangles contain point on source plane

    sourcetri = imagetri[indices] 
    sourcepos = np.sum(sourcetri.copy(),axis=1)/3.0 #list of the centroid coordinates for the triangles which contain the point 'image'
    realpos = np.array([(op.root(lambda x0: mapping(x0,image_loc,b),v)).x for v in sourcepos]) # use centroid coordinates as guesses for the actual root finding algorithm

    return [sourcepos, realpos]

def plot_graphs(critx,crity,caustics,stack,transformed,simplices,sourcepos,realpos,image,lowerend,upperend):
    plt.subplot(1,2,1) # image plane
    plt.title('Image Plane')
    plt.axis([lowerend,upperend,lowerend,upperend])
    plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis
    
    plt.plot(critx,crity,zorder=2) # plot of critical curve(s)

    plt.triplot(stack[:,0],stack[:,1], simplices,zorder=1) # plot of the Delaunay Triangulization

    plt.scatter(*zip(*sourcepos), marker='*', color='black',s=100, zorder=2) # plot of the (approximate/centroid) positions of the image
    plt.scatter(*zip(*realpos), marker='*', color='purple', s=100, zorder=2)

    plt.subplot(1,2,2) # source plane
    plt.title('Source Plane')
    plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

    plt.triplot(transformed[:,0],transformed[:,1], simplices, zorder=1) # plot of the transformed Delaunay Triangulization

    plt.scatter(*zip(image), marker='*', color='red', s= 200, zorder=2 ) # plot of (observed) image position

    plt.plot(caustics[:,0],caustics[:,1], color ='green', zorder=2) # plot of caustics?
    plt.show()


def run(critargs,carargs,polargs,show_plot=True,image=(random.random()*2-1,random.random()*2-1)):
    b = critargs[0]
    lowerend, upperend, spacing = carargs
    
    args = generate_ranges(critargs,carargs,polargs)
    critx, crity = args[0]
    caustics = args[1]
    x,y = args[2]
    r,theta = args[3]
    rs,thetas = args[4]

    stack, transformed, dpoints = transformations((x,y),(rs,thetas),spacing, b)
    
    #image = (random.random()*2-1,random.random()*2-1) #a random coordinate pair with x,y chosen from -1:1

    sourcepos, realpos = find_source(stack, transformed, dpoints.simplices.copy(), image, b)



    if show_plot:
        plot_graphs(critx,crity,caustics,stack,transformed,dpoints.simplices.copy(),sourcepos,realpos,image,lowerend,upperend)

critargs = [1.,200] #critical curve radius, number of points on critical curve function
carargs = [-2.5,2.5,0.5] # lower bound, upper bound, initial spacing (all 3 quantities apply the same to x and y axes)
polargs = [0.3,0.03,42] # outer radius, radius spacing, number of divisions in angle (for 360 degrees)
image = [0.6,-0.1] #image location if we want to specify

run(critargs,carargs,polargs)

# [1,200], [-2.5,2.5,0.5], [0.3,0.03,42] # sample runtime parameters


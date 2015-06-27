#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op
import random

import trinterior as trint

# lens model modules: 'models_list[modelargs[0]]' gives us the module to use
import sis
models_list = {'SIS':sis}   

def relation(x,y,coord,modelargs):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    model = models_list[modelargs[0]]
    dif = model.distance(x,y,coord,modelargs)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

def buildrelations(xran,yran,coord,modelargs):
    '''applies relation() on 2D range specified by xran and yran. Returns a 2D array.'''
    xx,yy = np.meshgrid(xran,yran, sparse=True) #sparse =True for sparse matrices, gives a row (x) and column (y) vector. 
    return np.transpose(relation(xx,yy,coord,modelargs))

def polartocar(r,th):
    '''convert polar coordinates to cartesian coordinates'''
    r  = np.array(r)
    th = np.array(th)
            
    return np.transpose([r*np.cos(th),r*np.sin(th)])
    
def mag_change(fir,sec,thr,frt):
    '''[[fir, sec],
        [thr, frt]] - returns true if there is a change in magnification in the box specified. 
    Note: the order of the points inputted does not matter as long as rotation symmetry is preserved. (i.e. the right points are adjacent to each other)'''
    return (fir*sec < 1 or sec*frt < 1 or frt*thr < 1 or thr*fir < 1)

RECURSEDEPTH = 4 #just here for future purpose when we want to set the max depth

def points(stack,xran,yran,coord,modelargs,n=0):
    '''collects the points around our magnification function and appends them to the 'stack'. Subdivides grid by 2 in each recursion.'''
    mat = buildrelations(xran,yran,coord,modelargs)
    xlen = len(xran)
    ylen = len(yran)
    for i in range(xlen-1):
        for j in range(ylen-1):
            if n != RECURSEDEPTH and mag_change(mat[i][j],mat[i][j+1],mat[i+1][j],mat[i+1][j+1]):
                points(stack,
                       np.linspace(xran[i],xran[i+1],3,endpoint=True),
                       np.linspace(yran[j],yran[j+1],3,endpoint=True),
                       coord,modelargs,n+1)
            else:
                stack.append([xran[i+1],yran[j]])
                stack.append([xran[i+1],yran[j+1]])
                stack.append([xran[i],yran[j]])
                stack.append([xran[i],yran[j+1]])

def points5(xran,yran,spacing,modelargs,recurse_depth=3,caustics=False):
    '''A vectorized approach to bulding a 'recursive' subgrid without recursion. Algorithm works by vectorizing each level of cell-size, handling each level in one complete calculation before proceeding to the next. '''
    x = xran[0:-1]
    y = yran[0:-1]
    xs, ys = np.tile(x,len(y)),np.repeat(y,len(x))

    grid_pairs = np.column_stack((np.tile(xran,len(yran)),np.repeat(yran,len(xran))))

    gridm_x_y = np.vstack((np.dstack((xs,xs,xs+spacing,xs+spacing)),np.dstack((ys,ys+spacing,ys,ys+spacing))))
    cells = np.transpose(gridm_x_y,[1,2,0])

    temp_cells = cells.copy()
    cells_sel = points5_wrapper(temp_cells,modelargs)
    if not caustics:
        output_pairs = grid_pairs.copy()

    for i in range(recurse_depth):
        temp_cells = subdivide_cells(cells_sel,spacing,i+1)
        cells_sel = points5_wrapper(temp_cells,modelargs)
        if not caustics:
            output_pairs = np.vstack((output_pairs,np.reshape(cells_sel,(-1,2))))

    if not caustics:
        return output_pairs
    else:
        return cells_sel
    
    
def points5_wrapper(cells,modelargs):
    '''Takes a list of cells and returns the cells for which the magnification changes sign. Function itself is a condensed command for three lines of code, which I did not want to write over and over again.'''
    cells_mag = mag_of_cells(cells,modelargs)
    mag_change_mask = cell_mag_change(cells_mag)
    return cells[mag_change_mask]
    
def mag_of_cells(cells,modelargs):
    '''Takes a list of cells and returns the magnification values for each point in the cells. Retains shape and order of the original list of cells.'''
    cells_x_y = np.reshape(cells,(-1,2))
    cells_x = cells_x_y[:,0]
    cells_y = cells_x_y[:,1]

    mag_x_y = relation(cells_x,cells_y,'car',modelargs)
    
    return np.reshape(mag_x_y,(-1,4))
    
def cell_mag_change(cells_mag):
    '''Takes a list of magnification values of a list of cells and returns a boolean mask for which the magnification changes across a cell. Uses numpy vectorization.'''
    fir,sec,thr,frt = cells_mag.T
    
    less1 = np.vstack((fir*sec,sec*frt,frt*thr,thr*fir)) < 1
    
    return np.any(less1,axis=0)

def subdivide_cells(cells,grid_spacing,cell_depth):
    '''Divides each cell in a list of cells (given by [[p1, p2],[p3, p4]]) into four smaller cells. On the cartesian grid, the points are ordered as so:
    p2--------p4         |
    |         |      II  |  I
    |         |     -----|-----
    |         |      III |  IV
    p1--------p3         |
    '''
    #cry/reimplement if we need to subdivide cells by any different than 4:1
    
    spacing = grid_spacing/ 2**cell_depth
    
    #below code uses broadcasting to 'shrink' each cell into a fourth of its size, but still retaining one of its vertices. This happens four times, shrinking each time towards one of the four vertices, leaving four quarter cells that together make up the original cell.
    quadrant1 = cells + [[spacing,spacing],[spacing,0],[0,spacing],[0,0]]
    quadrant2 = cells + [[0,spacing],[0,0],[-spacing,spacing],[-spacing,0]]
    quadrant3 = cells + [[0,0],[0,-spacing],[-spacing,0],[-spacing,-spacing]]
    quadrant4 = cells + [[spacing,0],[spacing,-spacing],[0,0],[0,-spacing]]

    return np.vstack((quadrant1,quadrant2,quadrant3,quadrant4))


def generate_ranges(carargs,polargs,modelargs):
    # generate our critical curve function, 200 points
    
    # initial cartesian grid, coarse,
    lowerend, upperend, spacing = carargs
    
    x = np.arange(lowerend,upperend+spacing,spacing)
    y = np.arange(lowerend,upperend+spacing,spacing)

    #intial polar grid
    rupper, rspacing, thetadivisions = polargs
    
    r = np.arange(0,rupper+rspacing,rspacing)
    theta = np.linspace(0,2*np.pi,thetadivisions,endpoint=False)

    ## the rs and thetas formed from a cartesian product. Used for vector operations.
    rs = np.tile(r,len(theta))
    thetas = np.repeat(theta,len(r))
    
    ## critical curves
    critx, crity = np.transpose(points5(x,y,spacing,modelargs,recurse_depth=8,caustics=True))

    ## caustics?
    model = models_list[modelargs[0]]
    causticsx,causticsy = np.transpose(model.carmapping(critx,crity,modelargs))

    return [ [critx, crity], [causticsx, causticsy], [x,y], [r, theta], [rs, thetas] ]


def transformations(car_ranges, pol_ranges, spacing, modelargs):

    x,y = car_ranges
    r,theta = pol_ranges
    model = models_list[modelargs[0]]
    
    polstack = [] # stack for holding the polar points #needed if we subgrid on polar grids
    carstack = [] # stack for holding the cartesian points
    stack = []
    #points(polstack,r,theta,'pol') # if we wanted to subgrid on the polar grid, uncomment
    #points2(carstack,x,y,'car',b) #generate subgrid on cartesian grid
    carstack = points5(x,y,spacing,modelargs)
    stack = np.array(stack)
    carstack = np.array(carstack) 
    polstack = np.array(polartocar(r,theta)) #comment if we subgrid on polar grids
    stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
    transformed = np.array(model.carmapping(stack[:,0],stack[:,1],modelargs)) #transform pairs from image to source plane

    dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

    return [stack, transformed, dpoints]



def find_source(stack, transformed, simplices, image_loc, modelargs):
    model = models_list[modelargs[0]]

    lenstri = transformed[simplices] #triangles on the source plane
    imagetri= stack[simplices] #triangles on the image plane

    indices = trint.find2(image_loc,lenstri.copy()) #list of which triangles contain point on source plane

    sourcetri = imagetri[indices] 
    sourcepos = np.sum(sourcetri.copy(),axis=1)/3.0 #list of the centroid coordinates for the triangles which contain the point 'image'
    realpos = np.array(
        [(op.root(lambda x0: model.mapping(x0,image_loc,modelargs),v)).x 
         for v in sourcepos]) # use centroid coordinates as guesses for the actual root finding algorithm

    return realpos

def plot_graphs(critx,crity,
                causticsx,causticsy,
                stack,transformed,
                simplices,
                realpos,image,
                lowerend,upperend):

    plt.subplot(1,2,1) # image plane
    plt.title('Image Plane')
    plt.axis([lowerend,upperend,lowerend,upperend])
    plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis
    
    plt.scatter(critx,crity, color='red', s=1, zorder=2) # plot of critical curve(s)

    plt.triplot(stack[:,0],stack[:,1], simplices, color='blue', zorder=1) # plot of the Delaunay Triangulization

    #plt.scatter(*zip(*sourcepos), marker='*', color='black',s=100, zorder=2) # plot of the (approximate/centroid) positions of the image
    plt.scatter(*zip(*realpos), marker='*', color='purple', s=100, zorder=2)

    plt.subplot(1,2,2) # source plane
    plt.title('Source Plane')
    plt.gca().set_aspect('equal', adjustable='box') # equal ratios on x and y axis

    plt.triplot(transformed[:,0],transformed[:,1], simplices, color='blue', zorder=1) # plot of the transformed Delaunay Triangulization

    plt.scatter(*zip(image), marker='*', color='red', s= 200, zorder=2 ) # plot of (observed) image position

    plt.scatter(causticsx, causticsy, color ='green', s=1, zorder=2) # plot of caustics
    plt.show()


def run(carargs,polargs,modelargs,show_plot=True,image=np.random.uniform(-1,1,2)):
    b = modelargs[1]

    lowerend, upperend, spacing = carargs
    
    args = generate_ranges(carargs,polargs,modelargs)
    critx, crity = args[0]
    causticsx,causticsy = args[1]
    x,y = args[2]
    r,theta = args[3]
    rs,thetas = args[4]

    stack, transformed, dpoints = transformations((x,y),(rs,thetas),spacing, modelargs)
    
    realpos = find_source(stack, transformed, dpoints.simplices.copy(), image, modelargs)

    if show_plot:
        plot_graphs(
            critx,crity,
            causticsx,causticsy,
            stack,transformed,dpoints.simplices.copy(),
            realpos,image,
            lowerend,upperend)


#### sample parameter run ####

pmodelargs = ['SIS',1.0] #model string, - for SIS: radius
pcarargs = [-2.5,2.5,0.5] # lower bound, upper bound, initial spacing (all 3 quantities apply the same to x and y axes)
ppolargs = [0.3,0.03,42] # outer radius, radius spacing, number of divisions in angle (for 360 degrees)
pimage = [0.6,-0.1] #image location -if- we want to specify

run(pcarargs,ppolargs,pmodelargs)

# [1,200], [-2.5,2.5,0.5], [0.3,0.03,42] # sample runtime parameters
# cProfile.run('run([1,200], [-2.5,2.5,0.5], [0.3,0.03,42],show_plot=False)')

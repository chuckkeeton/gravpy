#!/usr/bin/python

import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op

import trinterior as trint
import plots

# lens model modules: 'models_list[modelargs[0]]' gives us the module to use
import sie
models_list = {'SIE':sie}   

def process_modelargs(modelargs):
    '''Groups same-model mass components into 2D arrays for quicker(?) vectorized computation.'''

    models = [mass_component[0] for mass_component in modelargs]
    params = [mass_component[1:] for mass_component in modelargs]

    sortedmodels, invInd = np.unique(models,return_inverse=True)
    num_models = len(sortedmodels)
    
    ind = [np.argwhere(invInd==i).flatten() for i in range(num_models)]
    result = [ [sortedmodels[i], np.array(map(params.__getitem__,ind[i])).T.reshape((-1,1,len(ind[i])))] 
        for i in range(num_models)] #hairy... but it works....

    return result

def relation(x,y,modelargs):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    
    dif = distance(x,y,modelargs)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

def buildrelations(xran,yran,coord,modelargs):
    '''applies relation() on 2D range specified by xran and yran. Returns a 2D array.'''
    xx,yy = np.meshgrid(xran,yran, sparse=True) #sparse =True for sparse matrices, gives a row (x) and column (y) vector. 
    return np.transpose(relation(xx,yy,coord,modelargs))

def polartocar(r,th):
    '''convert polar coordinates to cartesian coordinates'''
    r  = np.array(r)
    th = np.array(th)
            
    return np.transpose(r*[np.cos(th),np.sin(th)])

def mapping(v,u,modelargs):
    '''Function used for root finding. u is image (non-variable), and v is the true position in the image plane (variable). Not intended for vectorized inputs. Returns the deflection vector along with the inverse-magnification matrix (aka the Jacobian).'''

    w,z = u #image
    x,y = v #actual position we want to find

    phiarr = potdefmag(x,y,modelargs).ravel()
    phix,phiy,phixx,phiyy,phixy = phiarr[1:6]

    return [[x-w,y-z] - np.array((phix,phiy)),
            [[1-phixx, -phixy],[-phixy,1-phiyy]]]

def carmapping(x,y,modelargs):
    '''mapping of cartesian coordinates from image to source plane'''
    
    phiarr = potdefmag(x,y,modelargs)
    phix,phiy = phiarr[1:3]
    
    return np.transpose([x,y] - np.array((phix,phiy)))

def distance(x,y,modelargs):
    '''returns the distance between the critical curve and the point '''

    phiarr = potdefmag(x,y,modelargs)
    phixx,phiyy,phixy = phiarr[3:6]
    
    return (1-phixx)*(1-phiyy)-phixy**2

def potdefmag(x,y,modelargs):
    '''The wrapper used to find the phi values given models' parameters. The output is (6,x) where x is the length of the x,y arguments given in the invocation. This command seeks out the correct module to contact for each model calculation.'''
    phi2Darray = []
    x = x.reshape((-1,1)) #broadcasting-ready
    y = y.reshape((-1,1)) #--^

    for mass_component in modelargs:
        model = models_list[mass_component[0]]
        args = mass_component[1]
        x0,y0 = args[1:3]
        te = args[4]
        
        #transformation from actual coords to natural coords (the frame/axes are rotated so there is no ellipticity angle in the calculation). This makes the expressions in the model modules simpler to calculate. 
        c = np.cos(te*np.pi/180);c2=c**2
        s = np.sin(te*np.pi/180);s2=s**2;sc=s*c
        xp = -s*(x-x0) + c*(y-y0)
        yp = -c*(x-x0) - s*(y-y0)

        pot,px,py,pxx,pyy,pxy = model.phiarray(xp,yp,args)

        # Inverse transformation back into desired coordinates. 
        newphix = -s*px-c*py
        newphiy = c*px-s*py 
        newphixx= s2*pxx+c2*pyy+2*sc*pxy
        newphiyy= c2*pxx+s2*pyy-2*sc*pxy
        newphixy= sc*(pyy-pxx)+(s2-c2)*pxy
        newarr = (pot,newphix,newphiy,newphixx,newphiyy,newphixy)

        phi2Darray.append(np.sum(newarr,axis=2)) #add up contributions from different parameters (from same model) into one result phiarray. Equivalent to a final sum over all different models simliar to the one after this for-loop, but this sum makes appending easier (due to dimensionality issues).
        
    return np.sum(phi2Darray,axis=0)

def cond_break(x,y,conds,function_calls,modelargs):
    n = np.max(x.shape)
    empty_shape = np.zeros((6,n,0),dtype='float64')

    allmodels = empty_shape.copy()

    zipped = np.dstack(conds,function_calls)

    for cond,function in zipped:
        temp_ind = np.flatnonzero(cond) 
        
        temp_args = np.take(modelargs,temp_ind,axis=2)
        temp_models = function(x[:,temp_args],y[:,temp_args],temp_args) if temp_args.size!= 0 else empty_shape
    
        allmodels = np.concatenate((allmodels,temp_models),axis=2)

    sorted_indices = np.argsort(conds)
    sorted_models = np.take(allmodels,sorted_indices,axis=2)

    return sorted_models
       

def mag_of_cells(cells,modelargs):
    '''Takes a list of cells and returns the magnification values for each point in the cells. Retains shape and order of the original list of cells.'''
    cells_x_y = cells.reshape((-1,2))
    cells_x = cells_x_y[:,0]
    cells_y = cells_x_y[:,1]

    # possible area to optimize later, here, phiarray calculates on same values twice.
    mag_x_y = relation(cells_x,cells_y,modelargs)
    
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

def points5_wrapper(cells,modelargs):
    '''Takes a list of cells and returns the cells for which the magnification changes sign. Function itself is a condensed command for three lines of code, which I did not want to write over and over again.'''
    cells_mag = mag_of_cells(cells,modelargs)
    mag_change_mask = cell_mag_change(cells_mag)
    return np.compress(mag_change_mask,cells,axis=0) #equivalent to cells[mag_change_mask] but faster

def points5(xran,yran,spacing,modelargs,recurse_depth=3,caustics_mode=False):
    '''A vectorized approach to bulding a 'recursive' subgrid without recursion. Algorithm works by vectorizing each level of cell-size, handling each level in one complete calculation before proceeding to the next. '''
    x = xran[0:-1]
    y = yran[0:-1]
    xs, ys = np.tile(x,len(y)),np.repeat(y,len(x))

    grid_pairs = np.column_stack((np.tile(xran,len(yran)),np.repeat(yran,len(xran))))

    gridm_x_y = np.vstack((np.dstack((xs,xs,xs+spacing,xs+spacing)),np.dstack((ys,ys+spacing,ys,ys+spacing))))
    cells = np.transpose(gridm_x_y,[1,2,0])

    temp_cells = cells.copy()
    cells_sel = points5_wrapper(temp_cells,modelargs)

    if not caustics_mode:
        output_pairs = grid_pairs.copy()

    for i in range(recurse_depth):
        temp_cells = subdivide_cells(cells_sel,spacing,i+1)
        cells_sel = points5_wrapper(temp_cells,modelargs)
        
        if not caustics_mode:
            output_pairs = np.vstack( (output_pairs,cells_sel.reshape((-1,2))) )

    if not caustics_mode:
        return output_pairs
    else:
        return np.mean(cells_sel,axis=1) # don't want the vertices of each cell; just the (center) of each cell


def generate_ranges(carargs,polargs,modelargs,caustics=True):
    '''Generates the sequences used for the other core & gridding functions.
    Returns an array containing:
    [ [x,y], [r,theta], [critx, crity], [causticsx, causticsy]]
    where critx,crity refers to the x and y values for the critical curves,
    causticsx, causticsy refers to the x and y values for the caustics,
    x,y are the cartesian ranges for the originial mesh-grid,
    r,theta are the cartesian x and y values for the supplementary polar grid(s).
    If caustics is set to False, then [critx, crity] and [causticsx, causticsy] are not returned-- Instead just [x,y]  and [r,theta] are returned.
    '''
        
    # initial cartesian grid, coarse,
    lowerend, upperend, spacing = carargs
    
    x = np.arange(lowerend,upperend+spacing,spacing)
    y = np.arange(lowerend,upperend+spacing,spacing)

    #supplemental polar grid(s)
    polargrids = np.reshape([],(0,2))
    
    for grid in polargs:
        center, rupper, rdivisions, thetadivisions = grid
        #rdivisions = int(rupper/rspacing)
        logr1 = np.log10(rupper+1)

        r = np.logspace(0,logr1,num=rdivisions)-1
        theta = np.linspace(0,2*np.pi,thetadivisions,endpoint=False)
        
        ## the rs and thetas formed from a cartesian product. Used for vector operations.
        temprs = np.tile(r,len(theta))
        tempthetas = np.repeat(theta,len(r))

        shiftedpairs = polartocar(temprs,tempthetas) + center

        polargrids = np.vstack((polargrids,shiftedpairs))
    
    #caustics and critical curves, if caustics == True
    if caustics:
        ## critical curves
        critx, crity = np.transpose(points5(x,y,spacing,modelargs,recurse_depth=8,caustics_mode=True))
        ## caustics
        causticsx,causticsy = np.transpose(carmapping(critx,crity,modelargs))
        return [ [x,y], polargrids, [critx, crity], [causticsx, causticsy] ]
    #otherwise just return x,y grid and polar grids
    else:
        return  [[x,y], polargrids]


def transformations(car_ranges, pol_ranges, spacing, modelargs, recurse_depth=3):
    '''Generates the subgridding points (more points around the critical curve), the transformations from image plane to source plane, and the Delaunay Triangulization object for plotting.'''
    x,y = car_ranges
      
    carstack = points5(x,y,spacing,modelargs,recurse_depth=recurse_depth)  #generate subgrid on cartesian grid
    polstack = np.array(pol_ranges)
    stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
    transformed = np.array(carmapping(stack[:,0],stack[:,1],modelargs)) #transform pairs from image to source plane

    dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

    return [stack, transformed, dpoints]


def find_source(stack, transformed, simplices, image_loc, modelargs):
    '''Employs the algorithm in the 'trinterior' module to find the positions of the image in the image plane. Returns the coordinate pair(s) in an array.'''
    

    lenstri = np.take(transformed,simplices,axis=0) #triangles on the source plane
    imagetri= np.take(stack,simplices,axis=0) #triangles on the image plane

    indices = trint.find2(image_loc,lenstri) #list of which triangles contain point on source plane

    sourcetri = imagetri[indices] 
    sourcepos = np.sum(sourcetri,axis=1)/3.0 #list of the centroid coordinates for the triangles which contain the point 'image'
    realpos = np.array(
        [(op.root(mapping,v,args=(image_loc,modelargs),jac=True)).x 
         for v in sourcepos]) # use centroid coordinates as guesses for the actual root finding algorithm

    return realpos


    


def run(carargs,polargs,modelargs,
        show_plot=True,caustics=True,image=np.random.uniform(-1,1,2),recurse_depth=3):
    '''The master command that wraps and executes all the commands to run a gridding example. Use this function (excusively) when using this module.'''
    
    lowerend, upperend, spacing = carargs
    bettermodelargs = process_modelargs(modelargs)

    args = generate_ranges(carargs,polargs,bettermodelargs,caustics=caustics)
    
    x,y = args[0]
    polargrids = args[1]

    if caustics:
        critx, crity = args[2]
        causticsx,causticsy = args[3]

    stack, transformed, dpoints = transformations((x,y),polargrids,spacing, bettermodelargs, recurse_depth=recurse_depth)
    
    realpos = find_source(stack, transformed, dpoints.simplices, image, bettermodelargs)

    if show_plot:
        plots.source_image_planes(
            stack,transformed,dpoints.simplices,
            realpos,image,
            lowerend,upperend,
            caustics=[[critx,crity],[causticsx,causticsy]] if caustics else False)



# cProfile.run('run([-2.5,2.5,0.5], [0.3,0.03,42], [\'SIS\',1.0],show_plot=False)')

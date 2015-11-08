#!/usr/bin/python
import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op

import trinterior as trint
import plots

# lens model modules: 'models_list[modelargs[0]]' gives us the module to use
import sie, alpha
models_list = {'SIE':sie,'alpha':alpha}   

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

def relation(x,y,modelargs,vec=False):
    '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
    
    dif = magnification(x,y,modelargs,vec=vec)
    return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

def polartocar(r,th):
    '''convert polar coordinates to cartesian coordinates'''
    r  = np.array(r)
    th = np.array(th)
            
    return np.transpose(r*[np.cos(th),np.sin(th)])

def mapping(v,u,modelargs,vec=False):
    '''Function used for root finding. u is image (non-variable), and v is the true position in the image plane (variable). Not intended for vectorized inputs. Returns the deflection vector along with the inverse-magnification matrix (aka the Jacobian).'''

    w,z = u #image
    x,y = v #actual position we want to find

    phiarr = potdefmag(x,y,modelargs,vec=vec,numexpr=False).ravel()
    phix,phiy,phixx,phiyy,phixy = phiarr[1:6]

    return [[x-w,y-z] - np.array((phix,phiy)),
            [[1-phixx, -phixy],[-phixy,1-phiyy]]]

def carmapping(x,y,modelargs,vec=False):
    '''mapping of cartesian coordinates from image to source plane'''
    
    phiarr = potdefmag(x,y,modelargs,vec=vec)
    phix,phiy = phiarr[1:3]
    
    return np.transpose([x,y] - np.array((phix,phiy)))

def magnification(x,y,modelargs,vec=False):
    '''returns the magnification of the points '''

    phiarr = potdefmag(x,y,modelargs,vec=vec)
    phixx,phiyy,phixy = phiarr[3:6]
    
    return (1-phixx)*(1-phiyy)-phixy**2

def ellipticity_calculation(x,y,mass_component,vec=False,numexpr=True):
    model = models_list[mass_component[0]]
    args = mass_component[1] if vec else np.array(mass_component[1:])
    x0,y0 = args[1:3]
    te = args[4]
    
    #transformation from actual coords to natural coords (the frame/axes are rotated so there is no ellipticity angle in the calculation). This makes the expressions in the model modules simpler to calculate. 
    c = np.cos(te*np.pi/180);c2=c*c
    s = np.sin(te*np.pi/180);s2=s*s;sc=s*c
    xp = -s*(x-x0) + c*(y-y0)
    yp = -c*(x-x0) - s*(y-y0)

    pot,px,py,pxx,pyy,pxy = model.phiarray(xp,yp,args,vec=vec,numexpr=numexpr)

    # Inverse transformation back into desired coordinates. 
    new_phix = -s*px-c*py
    new_phiy = c*px-s*py 
    new_phixx= s2*pxx+c2*pyy+2*sc*pxy
    new_phiyy= c2*pxx+s2*pyy-2*sc*pxy
    new_phixy= sc*(pyy-pxx)+(s2-c2)*pxy

    return  np.array((pot,new_phix,new_phiy,new_phixx,new_phiyy,new_phixy))

def potdefmag(xi,yi,modelargs,vec=False,numexpr=True):
    '''The wrapper used to find the phi values given models' parameters. The output is (6,x) where x is the length of the x,y arguments given in the invocation. This command seeks out the correct module to contact for each model calculation.'''
    phi2Darray = []

    #broadcasting-ready if vectorizing else (check for and) turn scalars into vectors of length 1
    x = np.expand_dims(xi,axis=1) if vec else np.atleast_1d(xi) 
    y = np.expand_dims(yi,axis=1) if vec else np.atleast_1d(yi) 
                

    for mass_component in modelargs:
        
        phiarray = ellipticity_calculation(x,y,mass_component,vec=vec,numexpr=numexpr)

        if vec:
            phi2Darray.extend(phiarray.transpose([2,0,1]))
        else:
            phi2Darray.append(phiarray)
            
    
    return np.sum(phi2Darray,axis=0)
    
    

def cond_break(x,y,modelargs,conds,function_calls,break_num=10):
    '''Experimental function that decides whether to vectorize multiple mass components or just iterate over a for loop. Current break is at 10 mass components'''
    num = x.shape[1]
    
    if num<break_num:
        return cond_break_for(x,y,modelargs,conds,function_calls)
    else:
        return cond_break_vec(x,y,modelargs,conds,function_calls)
    
    
def cond_break_for(x,y,modelargs,conds,function_calls):
    n = x.shape[1]

    ind = np.nonzero(conds)[1]
    ordered_funcs = np.array(function_calls)[ind]

    phiarrays = [ordered_funcs[i](x[:,i],y[:,i],modelargs[:,:,i]) for i in range(n)]
    
    return np.transpose(phiarrays,[1,2,0])
    
    
def cond_break_vec(x,y,modelargs,conds,function_calls):
    '''This code is essentially some filtering and ordering to seperate cases based on conditions and then recombine resulting phiarrays back into the same order they came in argument-wise. For example, in the 'sie' module, there is a split in calculations for elliptical and spherical cases. This wrapper will split the incoming arguments into their respective function calls and merge the outputs back together in the order they came in.'''

    n = x.shape[0]
    empty_shape = np.zeros((6,n,0),dtype='float64')

    allmodels = []
    allinds = [] 

    zipped = zip(conds,function_calls)

    for cond,function in zipped:
        temp_ind = np.flatnonzero(cond) 
        
        temp_args = np.take(modelargs,temp_ind,axis=2)
        temp_x,temp_y = np.take(x,temp_ind,axis=1),np.take(y,temp_ind,axis=1)
        temp_models = function(temp_x,temp_y,temp_args) if temp_args.size!= 0 else empty_shape
        
        allinds.append(temp_ind)
        allmodels.append(temp_models)
        
    npallinds = np.concatenate(allinds,axis=1)
    npallmodels = np.concatenate(allmodels,axis=2)
    sorted_indices = np.argsort(npallinds)
    sorted_models = np.take(npallmodels,sorted_indices,axis=2)
    
    return sorted_models
       

def mag_of_cells(cells,modelargs,vec=False):
    '''Takes a list of cells and returns the magnification values for each point in the cells. Retains shape and order of the original list of cells.'''
    num_points_in_cell = cells.shape[1] #usually 4, unless we're caching mag values (then it's 3)
    cells_x_y = cells.reshape((-1,2))
    cells_x = cells_x_y[:,0]
    cells_y = cells_x_y[:,1]

    mag_x_y = relation(cells_x,cells_y,modelargs,vec=vec)
    
    return np.reshape(mag_x_y,(-1,num_points_in_cell))
    

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

def points5_wrapper(cells,modelargs,vec=False,subdivide=False):
    '''Takes a list of cells and returns the cells for which the magnification changes sign. Function itself is a condensed command for three lines of code, which I did not want to write over and over again.'''
    temp_cells = cells if not subdivide else subdivide_cells(cells,subdivide[0],subdivide[1])
    cells_mag = mag_of_cells(temp_cells,modelargs,vec=vec)
    mag_change_mask = cell_mag_change(cells_mag)
    return np.compress(mag_change_mask,temp_cells,axis=0) #equivalent to cells[mag_change_mask] but faster

def for_points5_wrapper(cells,grid_spacing,cell_depth,modelargs,vec=False):
    '''Takes a list of cells and returns the cells for which the magnification changes sign. Function itself is a condensed command for three lines of code, which I did not want to write over and over again.'''
    subdivided_cells = subdivide_cells(cells,grid_spacing,cell_depth)
    cells_mag = mag_of_cells(subdivided_cells,modelargs,vec=vec)
    mag_change_mask = cell_mag_change(cells_mag)
    return np.compress(mag_change_mask,subdivided_cells,axis=0) #equivalent to cells[mag_change_mask] but faster

def for_points5_wrapper_cached(cells,mag_cells,grid_spacing,cell_depth,modelargs,vec=False):
    subdivided_cells = subdivide_cells(cells,grid_spacing,cell_depth)
    
    
    q1,q2,q3,q4 = np.split(subdivided_cells,4)
    q1,q2,q3,q4 = [np.delete(q,i,axis=1) for q,i in zip([q1,q2,q3,q4],[3,1,0,2])]
    
    m1,m2,m3,m4 = [mag_of_cells(q,modelargs,vec=vec) for q in [q1,q2,q3,q4]]
    
    c1,c2,c3,c4 = mag_cells.T #cache of old mag values (did not change)
    
    mag_combined = np.vstack([np.insert(m,i,c,axis=1) for m,i,c in zip([m1,m2,m3,m4],[3,1,0,2],[c4,c2,c1,c3])])
    
    mag_change_mask = cell_mag_change(mag_combined)
    return [np.compress(mag_change_mask,mag_combined,axis=0),np.compress(mag_change_mask,subdivided_cells,axis=0)]
    
def points5(xran,yran,spacing,modelargs,recurse_depth=3,caustics_mode=False,vec=False):
    '''A vectorized approach to bulding a 'recursive' subgrid without recursion. Algorithm works by vectorizing each level of cell-size, handling each level in one complete calculation before proceeding to the next. '''
    x = xran[0:-1]
    y = yran[0:-1]
    xs, ys = np.tile(x,len(y)),np.repeat(y,len(x))

    grid_pairs = np.column_stack((np.tile(xran,len(yran)),np.repeat(yran,len(xran))))

    gridm_x_y = np.vstack((np.dstack((xs,xs,xs+spacing,xs+spacing)),np.dstack((ys,ys+spacing,ys,ys+spacing))))
    cells = np.transpose(gridm_x_y,[1,2,0])

    # we don't want to subdivide the first iteration
    cells_mag = mag_of_cells(cells,modelargs,vec=vec)
    mag_change_mask = cell_mag_change(cells_mag)
    
    cells_sel = np.compress(mag_change_mask,cells,axis=0) #equivalent to cells[mag_change_mask] but faster
    cells_mag = np.compress(mag_change_mask,cells_mag,axis=0) # = cells_mag[mag_change_mask]
    output_pairs = []
    
    for i in range(recurse_depth):
        
        cells_mag,cells_sel = for_points5_wrapper_cached(cells_sel,cells_mag,spacing,i+1,modelargs,vec=vec)
        
        if not caustics_mode:
            output_pairs.append(cells_sel)

    if not caustics_mode:
        output_pairs = np.vstack(output_pairs).reshape((-1,2))
        return np.vstack((grid_pairs,output_pairs))
    else:
        return np.mean(cells_sel,axis=1) # don't want the vertices of each cell; just the (center) of each cell


def generate_ranges(carargs,polargs,modelargs,caustics=True,vec=False):
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
        critx, crity = np.transpose(points5(x,y,spacing,modelargs,recurse_depth=8,caustics_mode=True,vec=vec))
        ## caustics
        causticsx,causticsy = np.transpose(carmapping(critx,crity,modelargs,vec=vec))
        return [ [x,y], polargrids, [critx, crity], [causticsx, causticsy] ]
    #otherwise just return x,y grid and polar grids
    else:
        return  [[x,y], polargrids]


def transformations(car_ranges, pol_ranges, spacing, modelargs, recurse_depth=3, vec=False):
    '''Generates the subgridding points (more points around the critical curve), the transformations from image plane to source plane, and the Delaunay Triangulization object for plotting.'''
    x,y = car_ranges
      
    carstack = points5(x,y,spacing,modelargs,recurse_depth=recurse_depth,vec=vec)  #generate subgrid on cartesian grid
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
    sourcepos = np.mean(sourcetri,axis=1) #list of the centroid coordinates for the triangles which contain the point 'image'
    realpos = np.array(
        [(op.root(mapping,v,args=(image_loc,modelargs),jac=True,tol=1e-4)).x
         for v in sourcepos]) # use centroid coordinates as guesses for the actual root finding algorithm

    return realpos

def run(carargs,polargs,modelargs,
        show_plot=True,caustics=True,image=np.random.uniform(-1,1,2),recurse_depth=3,vec=False):
    '''The master command that wraps and executes all the commands to run a gridding example. Use this function (excusively) when using this module.'''
    
    lowerend, upperend, spacing = carargs
    
    tempmodelargs = process_modelargs(modelargs) if vec else modelargs
    
    args = generate_ranges(carargs,polargs,tempmodelargs,caustics=caustics,vec=vec)
    
    x,y = args[0]
    polargrids = args[1]

    if caustics:
        critx, crity = args[2]
        causticsx,causticsy = args[3]

    stack, transformed, dpoints = transformations((x,y),polargrids,spacing, tempmodelargs, recurse_depth=recurse_depth)
    
    realpos = find_source(stack, transformed, dpoints.simplices, image, tempmodelargs)

    stackx,stacky=np.transpose(stack)
    tranx,trany = np.transpose(transformed)
    mag = magnification(stackx,stacky,tempmodelargs)
    
    
    if show_plot:
        plots.source_image_planes(
            stack,transformed,dpoints.simplices,
            realpos,image,
            lowerend,upperend,
            caustics=[[critx,crity],[causticsx,causticsy]] if caustics else False)
#        plots.mag_map(tranx,trany,mag,dpoints.simplices)

   




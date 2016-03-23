#!/usr/bin/python
import numpy as np
from scipy.spatial import Delaunay 
import scipy.optimize as op
import numexpr as ne
import trinterior as trint
import plots

class gravlens:
        
    def __init__(self,carargs,polargs,modelargs,show_plot=True,include_caustics=True,image=np.random.uniform(-1,1,2),recurse_depth=3,caustics_depth=8):
        self.carargs = carargs
        self.xspacing = carargs[0][2]
        self.yspacing = carargs[1][2]
        self.polargs = polargs
        self.modelargs = modelargs
        self.show_plot = show_plot
        self.include_caustics = include_caustics
        self.image = image
        self.recurse_depth = recurse_depth
        self.caustics_depth= caustics_depth
        
        self.num_eval = 0
        
    def relation(self,x,y):
        '''tells us if the point pair (x,y) is outside, inside, or on the critical curve'''
        
        dif = self.magnification(x,y)
        return np.sign(dif) # -1 for inside, +1 for outside, 0 for exactly on

    @staticmethod
    def polartocar(r,th):
        '''convert polar coordinates to cartesian coordinates'''
        r  = np.array(r)
        th = np.array(th)
        
        return np.transpose(r*[np.cos(th),np.sin(th)])
    
    def mapping(self,v):
        '''Function used for root finding. \'v\' is the true position in the image plane (what we want to find). Not intended for vectorized inputs. Returns the deflection vector along with the inverse-magnification matrix (aka the Jacobian).'''
    
        w,z = self.image
        x,y = v #actual position we want to find
    
        phiarr = self.potdefmag(x,y,numexpr=False).ravel()
        phix,phiy,phixx,phiyy,phixy = phiarr[1:6]
    
        return [[x-w-phix,y-z-phiy],
                [[1-phixx, -phixy],
                 [-phixy ,1-phiyy]]]
    
    
    def carmapping(self,x,y):
        '''mapping of cartesian coordinates from image to source plane'''
        print "--Mapping Call--"
        phiarr = self.potdefmag(x,y)
        phix,phiy = phiarr[1:3]
        
        return np.transpose([x-phix,y-phiy])
    
    def magnification(self,x,y):
        '''returns the magnification of the points '''
        print "--Magnification Call--"
        phiarr = self.potdefmag(x,y)
        phixx,phiyy,phixy = phiarr[3:6]
        
        return (1-phixx)*(1-phiyy)-phixy**2
    
    def potdefmag(self,xi,yi,numexpr=True):
        '''The wrapper used to find the phi values given models' parameters. The output is (6,x) where x is the length of the x,y arguments given in the invocation. This command seeks out the correct module to contact for each model calculation.'''
        phi2Darray = []
    
        #turn scalars into vectors of length 1
        x = np.atleast_1d(xi) 
        y = np.atleast_1d(yi) 

        self.num_eval += x.size

        print "Evaluating %d points..." % x.size

        for mass_component in self.modelargs:
            
            phiarray = mass_component.phiarray(x,y,numexpr=numexpr) 
            phi2Darray.append(phiarray)
                
        return np.sum(phi2Darray,axis=0)
    
    
    def mag_of_cells(self,cells):
        '''Takes a list of cells and returns the magnification values for each point in the cells. Retains shape and order of the original list of cells.'''
        num_points_in_cell = cells.shape[1] #usually 4, unless we're caching mag values (then it's 3)
        cells_x_y = cells.reshape((-1,2))
        cells_x = cells_x_y[:,0]
        cells_y = cells_x_y[:,1]
    
        mag_x_y = self.relation(cells_x,cells_y)
        
        return np.reshape(mag_x_y,(-1,num_points_in_cell))
    
    @staticmethod
    def cell_mag_change(cells_mag):
        '''Takes a list of magnification values of a list of cells and returns a boolean mask for which the magnification changes across a cell. Uses numpy vectorization.'''
        fir,sec,thr,frt = cells_mag.T

        with np.errstate(invalid='ignore'): 
            less1 = np.vstack((fir*sec,sec*frt,frt*thr,thr*fir)) < 1
        
        output = np.any(less1,axis=0)
                
        if np.count_nonzero(output) == 0: #np.all doesn't catch the error
            raise ValueError("Magnification does not change across grid, change grid parameters so that critical curves are seen.")
        else:
            return output

    
    def subdivide_cells(self,cells,cell_depth):
        '''Divides each cell in a list of cells (given by [[p1, p2],[p3, p4]]) into four smaller cells. On the cartesian grid, the points are ordered as so:
        p2--------p4         |
        |         |      II  |  I
        |         |     -----|-----
        |         |      III |  IV
        p1--------p3         |
        '''
        #cry/reimplement if we need to subdivide cells by any different than 4:1

        dx = self.xspacing/ 2**cell_depth
        dy = self.yspacing/ 2**cell_depth
        
        #below code uses broadcasting to 'shrink' each cell into a fourth of its size, but still retaining one of its vertices. This happens four times, shrinking each time towards one of the four vertices, leaving four quarter cells that together make up the original cell.
        quadrant1 = cells + [[dx, dy],[dx,  0],[0  , dy],[0  ,  0]]
        quadrant2 = cells + [[0 , dy],[0 ,  0],[-dx, dy],[-dx,  0]]
        quadrant3 = cells + [[0 ,  0],[0 ,-dy],[-dx,  0],[-dx,-dy]]
        quadrant4 = cells + [[dx,  0],[dx,-dy],[0  ,  0],[0  ,-dy]]
    
        return np.vstack((quadrant1,quadrant2,quadrant3,quadrant4))
    
    def for_points5_wrapper_cached(self,cells,mag_cells,cell_depth):
        '''Function that subdivdes given cells, computes the magnification values for new points, merges the magnification values from \'mag_cells\' with the newly computed magnifcation values, and returns the magnification values and cells where a critical curve was detected.'''
                
        subdivided_cells = self.subdivide_cells(cells,cell_depth)
        
        q1,q2,q3,q4 = np.split(subdivided_cells,4) #need the points for each quadrant
        
        q1,q2,q3,q4 = [np.delete(q,i,axis=1) for q,i in zip([q1,q2,q3,q4],[3,1,0,2])] # remove the cells that stay the same
        
        m1,m2,m3,m4 = [self.mag_of_cells(q) for q in [q1,q2,q3,q4]] # mag of points above
        
        c1,c2,c3,c4 = mag_cells.T #cache of old mag values (did not change)
        
        mag_combined = np.vstack([np.insert(m,i,c,axis=1) for m,i,c in zip([m1,m2,m3,m4],[3,1,0,2],[c4,c2,c1,c3])]) #combine old values with newly calculated values
        
        mag_change_mask = self.cell_mag_change(mag_combined) #which cells had a change in magnification
        
        # return mag values of selected cells (used for the next level in gridding) and the cells themselves that had a mag change
        return [np.compress(mag_change_mask,mag_combined,axis=0),np.compress(mag_change_mask,subdivided_cells,axis=0)]
    
    def points5(self,xran,yran):#,caustics_mode=False):
        '''A vectorized approach to bulding a 'recursive' subgrid without recursion. Algorithm works by vectorizing each level of cell-size, handling each level in one complete calculation before proceeding to the next. '''
                
        x = xran[0:-1]
        y = yran[0:-1]
        xs, ys = np.tile(x,len(y)),np.repeat(y,len(x))
    
        grid_pairs = np.column_stack((np.tile(xran,len(yran)),np.repeat(yran,len(xran))))
    
        gridm_x_y = np.vstack((np.dstack((xs,xs,xs+self.xspacing,xs+self.xspacing)),np.dstack((ys,ys+self.yspacing,ys,ys+self.yspacing))))
        cells = np.transpose(gridm_x_y,[1,2,0])
    
        # we don't want to subdivide the first iteration
        cells_mag = self.mag_of_cells(cells)
        mag_change_mask = self.cell_mag_change(cells_mag)
        
        cells_sel = np.compress(mag_change_mask,cells,axis=0) #equivalent to cells[mag_change_mask] but faster
                
        cells_mag = np.compress(mag_change_mask,cells_mag,axis=0) # = cells_mag[mag_change_mask]
        output_pairs = []
        
        depth = self.caustics_depth if self.include_caustics else self.recurse_depth
        for i in range(depth):
            
            cells_mag,cells_sel = self.for_points5_wrapper_cached(cells_sel,cells_mag,i+1)
            
            #if not caustics_mode:
            output_pairs.append(cells_sel)
    
        #if not caustics_mode:
        output_pairs = np.vstack(output_pairs).reshape((-1,2))
        if self.include_caustics:
            self.critical_lines = np.transpose(np.mean(cells_sel,axis=1)) # don't want the vertices of each cell; just the (center) of each cell
            
        return np.vstack((grid_pairs,output_pairs))
        
    
    
    def generate_ranges(self):
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
        [xlowerend, xupperend, xspacing],[ylowerend,yupperend,yspacing] = self.carargs

        x = np.arange(xlowerend,xupperend+xspacing,xspacing)
        y = np.arange(ylowerend,yupperend+yspacing,yspacing)
    
        #supplemental polar grid(s)
        polargrids = np.reshape([],(0,2))
        
        for grid in self.polargs:
            center, rupper, rdivisions, thetadivisions = grid
            
            logr1 = np.log10(rupper+1)
    
            r = np.logspace(0,logr1,num=rdivisions)-1
            theta = np.linspace(0,2*np.pi,thetadivisions,endpoint=False)
            
            ## the rs and thetas formed from a cartesian product. Used for vector operations.
            temprs = np.tile(r,len(theta))
            tempthetas = np.repeat(theta,len(r))
    
            shiftedpairs = self.polartocar(temprs,tempthetas) + center
    
            polargrids = np.vstack((polargrids,shiftedpairs))
            
        #caustics and critical curves, if include_caustics == True
#        if self.include_caustics:
#            ## critical curves
#            critx, crity = np.transpose(self.points5(x,y,caustics_mode=True))
#            ## caustics
#            causticsx,causticsy = np.transpose(self.carmapping(critx,crity))
#            return [ [x,y], polargrids, [critx, crity], [causticsx, causticsy] ]
        #otherwise just return x,y grid and polar grids
#        else:
        return  [[x,y], polargrids]
    
    
    def transformations(self,car_ranges, pol_ranges):
        '''Generates the subgridding points (more points around the critical curve), the transformations from image plane to source plane, and the Delaunay Triangulization object for plotting.'''
        
        x,y = car_ranges
        
        carstack = self.points5(x,y)  #generate subgrid on cartesian grid
        polstack = np.array(pol_ranges)
        stack = np.concatenate((carstack,polstack),axis=0) #combine list of cartesian and polar pairs
        transformed = self.carmapping(stack[:,0],stack[:,1]) #transform pairs from image to source plane
        if self.include_caustics:
            criticalx, criticaly = self.critical_lines
            causticsx, causticsy = np.transpose(self.carmapping(criticalx,criticaly))
            self.caustics = [[criticalx,criticaly],[causticsx,causticsy]]
        else:
            self.caustics = None

            
        dpoints = Delaunay(stack) # generate Delaunay object for triangulization/triangles

        self.stack = stack
        self.transformed = transformed
        self.dpoints = dpoints
        
        
    def find_source(self):
        '''Employs the algorithm in the 'trinterior' module to find the positions of the image in the image plane. Returns the coordinate pair(s) in an array.'''

        simplices = self.dpoints.simplices
    
        lenstri = np.take(self.transformed,simplices,axis=0) #triangles on the source plane
        imagetri= np.take(self.stack,simplices,axis=0) #triangles on the image plane
    
        indices = trint.find2(self.image,lenstri) #list of which triangles contain point on source plane
    
        sourcetri = imagetri[indices] 
        sourcepos = np.mean(sourcetri,axis=1) #list of the centroid coordinates for the triangles which contain the point 'image'
        realpos = np.array(
            [(op.root(self.mapping,v,jac=True,tol=1e-4)).x
             for v in sourcepos]) # use centroid coordinates as guesses for the actual root finding algorithm
    
        self.realpos = realpos
    
    def run(self):
        '''The master command that wraps and executes all the commands to run a gridding example. Use this function (excusively) when using this module.'''
                        
        self.validate_arguments()
        
        args = self.generate_ranges()
        
        x,y = args[0]
        polargrids = args[1]
    
        #if self.include_caustics:
            #critx, crity = args[2]
            #causticsx,causticsy = args[3]
            #self.caustics = [[critx,crity],[causticsx,causticsy]]
        #else:
            
        self.transformations((x,y),polargrids)
        
        self.find_source()
        
        if self.show_plot:
            self.plot()

    def validate_arguments(self):
        if not np.array(self.carargs).shape==(2,3):
            raise AssertionError("'carargs' are not of the correct shape")
        if not type(self.polargs[0]) is list:
            raise AssertionError("'polargs' is not a list of list(s)")
        if not type(self.modelargs) is list:
            raise AssertionError("'modelargs' is not a list of models")

        
    def plot(self):

        [xlowerend, xupperend, xspacing],[ylowerend,yupperend,yspacing] = self.carargs
        
        plots.source_image_planes(
                self.stack,self.transformed,self.dpoints.simplices,
                self.realpos,self.image,
                xlowerend,xupperend,ylowerend,yupperend,
                caustics=self.caustics)
    
    
    
    

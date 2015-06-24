#!/usr/bin/python

import numpy as np
import numpy.linalg as la

def sub(u,v):
    '''Subtract two vectors '''
    return (u[0]-v[0],u[1]-v[1])

def det(u,v):
    '''Return the determinant of the 2-D matrix formed by u and v'''
    return u[0]*v[1]-u[1]*v[0]

def inside_triangle(P,tri):
    '''Explicit expression of the components in the algebra of finding barycentric coordiantes. Plugging in an x and y value for v0,v1,v2 reduces the computation to three determinants from before which computed six dot products.'''
    # Code derived from http://www.blackpawn.com/texts/pointinpoly/
    # Uses barycentric coordinates to find whether the point is inside the triangle

    A,B,C = tri
    
    v0 = sub(C,A)
    v1 = sub(B,A)
    v2 = sub(P,A)

    invDenom = 1.0/det(v1,v0)
    
    u = det(v2,v1) * (-invDenom)
    v = det(v2,v0) *   invDenom

    return u >= 0 and v >= 0 and u + v < 1

def find(v,triangles):
    '''Given a point v and a list of triangles defined by three verticies, this function will return the indicies of the triangles in which the point is inside the triangle.'''

    indices = [i for i,triangle in enumerate(triangles) if inside_triangle(v,triangle)]
            
    return indices

def np_inside_triangle(A,B,C,P):
    '''Explicit expression of the components in the algebra of finding barycentric coordiantes. Plugging in an x and y value for v0,v1,v2 reduces the computation to three determinants from before which computed six dot products. This particular function uses numpy vector operations, where the inputs are vectors of pairs. 'a,b,c' are the lists of pairs of the triangle vertices, and p is still a single point for which we are searching for.'''
    # Code derived from http://www.blackpawn.com/texts/pointinpoly/
    # Uses barycentric coordinates to find whether the point is inside the triangle
    # Uses numpy vector operations

    v0 = C-A
    v1 = B-A
    v2 = P-A
    
    with np.errstate(divide='ignore',invalid='ignore'): 
        #my understanding is that invalid values arise from floating point errors, so nothing to worry about? still gives good results...
        invDenom = 1.0/la.det(np.transpose(np.dstack((v1,v0)),[0,2,1])) #divide/0 
    
        u = la.det(np.transpose(np.dstack((v2,v1)),[0,2,1])) * (-invDenom) #invalid value in det()
        v = la.det(np.transpose(np.dstack((v2,v0)),[0,2,1])) *   invDenom  #invalid value in det() 

        coords = np.vstack((u,v,1-u-v))
        indices = np.nonzero(np.all(coords>=0,axis=0)) #invalid value in greater_equal
    
    return indices

def find2(p,triangles):
    '''Wrapper for the algorithm that uses numpy vector operations. Use this one for better speed!'''
    a,b,c = np.transpose(triangles,[1,0,2])
    return np_inside_triangle(a,b,c,p)

def np_inside_triangle_dots(A,B,C,P):
    v0 = C-A
    v1 = B-A
    v2 = P-A

    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom

    coords = np.vstack((u,v,1-u-v))
    indices = np.nonzero(np.all(coords>=0,axis=0)) #invalid value in greater_equal
    
    return indices
    
### the following code is for testing purposes ###


sampletriangles = np.load("triangles.npy") #load a working example of a list of triangles formed by Delaunay Triangulization, for use in the test function below
def test():
    '''a function used to time the 'inside_triangle' algorithm'''
    find2([0.5,0.5],sampletriangles)

# run from the command line to time the function:
## python -m timeit -v -s'import trinterior' 'trinterior.test()'

# Speedups
# dot products -> determinants in inside_triangle(): 80 msec to 60 msec
# for loop -> list comprehension in find(): 60 msec to 56 msec
# find() -> find2() -- full transition to numpy vector operations: 56 msec to 15 msec

    

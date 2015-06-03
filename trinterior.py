#!/usr/bin/python

import numpy as np


def dot(u,v):
    '''Dot Product of two vectors'''
    return u[0]*v[0] + u[1]*v[1]

def sub(u,v):
    '''Subtract two vectors '''
    return (u[0]-v[0],u[1]-v[1])

def condition(P, tri):
    '''Returns a Boolean whether P is inside the triangle defined by vertices A,B,C'''
    # Code shamelessly copied from http://www.blackpawn.com/texts/pointinpoly/
    # Uses barycentric coordinates to find whether the point is inside the triangle
    (A,B,C) = tri
    

    v0 = sub(C,A)
    v1 = sub(B,A)
    v2 = sub(P,A)
    
    dot00 = dot(v0,v0)
    dot01 = dot(v0,v1)
    dot02 = dot(v0,v2)
    dot11 = dot(v1,v1)
    dot12 = dot(v1,v2)
    
    denom = (dot00 * dot11 - dot01 * dot01)
    if denom == 0: return False
    
    u = (dot11 * dot02 - dot01 * dot12) / denom
    v = (dot00 * dot12 - dot01 * dot02) / denom
    
    return u >= 0 and v >= 0 and u + v < 1


def find(v,triangles):
    '''Given a point v and a list of triangles defined by three verticies, this function will return the indicies of the triangles in which the point is inside the triangle.'''
    length = len(triangles)
    indices = []
    for i in range(length):
        if condition(v,triangles[i]):
            indices.append(i)
    
    return indices
        
### the following code is for testing purposes ###


sampletriangles = np.load("triangles.npy") #load a working example of a list of triangles formed by Delaunay Triangulization, for use in the test function below
def test():
    '''a function used to time the 'condition' algorithm'''
    find([0.5,0.5],sampletriangles)

# run from the command line to time the function:
## python -m timeit -v -s'import trinterior' 'trinterior.test()'


    

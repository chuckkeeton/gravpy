#!/usr/bin/python

import numpy as np


def dot(u,v):
    return u[0]*v[0] + u[1]*v[1]

def sub(u,v):
    return (u[0]-v[0],u[1]-v[1])

def condition(P, A, B, C):
    '''Returns a Boolean whether P is inside the triangle defined by vertices A,B,C'''
    # Code shamelessly copied from http://www.blackpawn.com/texts/pointinpoly/

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
    length = len(triangles)
    indices = []
    for i in range(length):
        if condition(v,*triangles[i]):
            indices.append(i)
    
    return indices
        

sampletriangles = np.load("triangles.npy")
def test():
    find([0.5,0.5],sampletriangles)


    
#u = (0.3,0.3)
#x = (0.1,0.1)
#y = (0.5,5)
#z = (3,0)

#print tot(u,x,y,z)

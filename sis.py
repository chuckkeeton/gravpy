import numpy as np

#param list: radius,

def mapping(v,u,modelargs):
    b = modelargs[1]
    
    w,z = u #image
    x,y = v #actual position we want to find
    return [x - b*x/np.sqrt(x**2 + y**2) - w, y - b*y/np.sqrt(x**2 + y **2) - z]

def carmapping(x,y,modelargs):
    '''mapping of cartesian coordinates from image to source plane'''
    b = modelargs[1]

    x = np.array(x)
    y = np.array(y)
    return np.transpose([x - b*x/np.sqrt(x**2 + y**2), y - b*y/np.sqrt(x**2 + y **2)])

def polmapping(r,th,modelargs):
    '''mapping of polar coordinates from image to source plane'''
    b = modelargs[1]

    r  = np.array(r)
    th = np.array(th)
    
    return np.transpose( [np.cos(th)*(r-b), np.sin(th)*(r-b)] )

def distance(x,y,coord,modelargs):
    '''returns the distance between the critical curve and the point '''
    b = modelargs[1]
    if coord == 'car':
        return np.sqrt(x**2 + y**2) - b
    elif coord == 'pol':
        return x - b
    else:
        raise NameError('Unrecogonized Coordinate System')

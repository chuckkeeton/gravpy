from gravpy import gravlens
from gravpy.models import SIE

carargs = [[-2.5,2.5,0.5],[-2.5,2.5,0.5]] # lower bound, upper bound, initial spacing (two sets--for x and y axes)
polarargs = [[(0,0),1,10,42]] # center position (coordinate pair), outer radius, number of divisions in radius, number of divisions in angle (for 360 degrees) 
modelargs = [SIE(1,0,0,0,0,0)]

example = gravlens(carargs,polarargs,modelargs)
example.run()

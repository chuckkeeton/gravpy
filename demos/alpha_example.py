from gravpy import gravlens
from gravpy.models import Alpha

carargs = [[-2.5,2.5,0.5],[-2.5,2.5,0.5]] # lower bound, upper bound, initial spacing (two sets--for x and y axes)
polarargs = [[(-0.5,-0.5),1,10,42],[(0.5,0.5),1,10,42]] # center position (coordinate pair), outer radius, number of divisions in radius, number of divisions in angle (for 360 degrees) 
modelargs = [Alpha(1, -0.5, -0.5, 0, 0, 0, 1.5),
             Alpha(0.5, 0.5, 0.5, 0.5, 45, 0.25, -1), ]

example = gravlens(carargs,polarargs,modelargs)
example.run()

from gravpy import Gravlens
from gravpy.models import Alpha

# lower bound, upper bound, initial spacing (two sets--for x and y axes)
carargs = [[-2.5, 2.5, 0.5], [-2.5, 2.5, 0.5]]
# center position (coordinate pair), outer radius, number of divisions in radius,
# number of divisions in angle (for 360 degrees)
polarargs = [
        [(-0.5, -0.5), 1, 10, 42],
        [(0.5, 0.5), 1, 10, 42]
    ]
modelargs = [
        Alpha(1, -0.5, -0.5, 0, 0, 0, 1),
        Alpha(0.5, 0.5, 0.5, 0.5, 45, 0.25, 1)
    ]

example = Gravlens(carargs, polarargs, modelargs)
example.run()

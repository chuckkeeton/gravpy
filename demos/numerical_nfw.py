from gravpy import NFW
from gravpy import Gravlens

carargs = [
        [-2.5, 2.5, 0.5],
        [-2.5, 2.5, 0.5]
    ]

polarargs = [[(0, 0), 1, 10, 42]]
modelargs = [NFW(1, 0.5, 0, 0, 0, 0, 0)]

example = Gravlens(carargs, polarargs, modelargs)
example.run()

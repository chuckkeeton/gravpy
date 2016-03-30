from gravpy import Gravpy
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
from models import SIE, Alpha, NFW

#############################
# Actual parameters needed for evaluation

# models' args, - look in models.py for parameter order
# (simliar to gravlens param order minus the two shear parameters)
pmodelargs = [
    # SIE(0.75, -0.45, -0.5, 0.1, 20, 0.1),
    # SIE(0.5, 1.5, 1.5, 0, 0.2, 0.0),
    # SIE(0.1, 0.25, 0.25, 0.1, 0, 0),
    # SIE(0.3, -0.4, 0.6, 0.5, 45, 0),
    NFW(2, 0, 0, 0, 0, 0.5, 1, 0.5)
]

# center position (coordinate pair), outer radius, number of divisions in radius,
# number of divisions in angle (for 360 degrees)
ppolargs = [
        [(0.0, 0.0), 0.9, 10, 42],
    ]
# lower bound, upper bound, initial spacing (two sets--for x and y axes)
pcarargs = [
        [-2.5, 2.5, 0.5],
        [-2.5, 2.5, 0.5]
    ]
pimage = [0.25, 0.25]  # image location -if- we want to specify
##############################

# If we want an image with runtime stats, set bool to true, otherwise runs the statement in the else branch.
callgraph = False
filepath = 'runs/class2.png'  # where we want to save the output image with the runtime breakdown
if callgraph:
    with PyCallGraph(output=GraphvizOutput(output_file=filepath)):
        example = Gravpy(pcarargs, ppolargs, pmodelargs, image=pimage, show_plot=False, include_caustics=False)
        example.run()
else:
    example = Gravpy(pcarargs, ppolargs, pmodelargs)
    example.run()
    # print pmodelargs[0].integrator.i.getcacheinfo()
    # print pmodelargs[0].integrator.j0.getcacheinfo()
    # print pmodelargs[0].integrator.j1.getcacheinfo()
    # print pmodelargs[0].integrator.k0.getcacheinfo()
    # print pmodelargs[0].integrator.k1.getcacheinfo()
    # print pmodelargs[0].integrator.k2.getcacheinfo()

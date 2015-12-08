#Welcome!

This project is a Python implementation of the general lens solver GravLens. 

*To run the program, execute the 'example_run.py' file.*

At this current version, this project does the following:
* Take the magnification function for given model(s) and compute the transformation of coordinates from the image plane to the source plane.
* Given a point in the source plane, find the exact position(s) in the image plane. 

# Supported Models
Alpha = 1 (Singular Isothermal Ellipsoid)

Models that are not fully functional in this version:
Alpha = -1 (Plummer Model)

# Dependencies:
Make sure to install these libraries before running!

* SciPy
* NumPy
* matplotlib
* pycallgraph (this library can be excluded if you are not interested in seeing runtimes, make sure to comment out the relevant lines in `example_run.py`)



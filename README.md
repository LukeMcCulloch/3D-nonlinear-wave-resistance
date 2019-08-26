# 3D nonlinear wave resistance


## This code is not finished!

But there is a lot of good stuff.  Most everything is there.

Here's what I remember it still needing:

 * update the hull mesh properly when the free surface changes
 * It's been 6 years since I worked on this.  I don't think it converged well before I had to stop work on it.  
 
 
 ## Two versions:
  * This one, I believe, implements the transom stern condition of Raven 1996.
  * Another version keeps the standard boundary conditions for a canoe body.



## Build

The included makefile builds the code on Linux with gfortran.


## Run

Run the an included offshore supply vessel hull example with

$ ./flowsolve ./hulls/OSV_Dec2015 AnyName.out .2 5
 
 where AnyName names an output file.  (Though most output goes to VTK paraview files.)
 
 .2 is the Froude number
 
 5 is the number of iterations allowed.  (I think this might not do anything right now.)

./flowsolve runs the executable
fifi.dat selects the included panel input file of the wigely hull
test1.out sets an output file
.2 sets the Froude number to 0.2

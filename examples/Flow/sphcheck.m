## check results from potential flow calculation

## set sphere radius (this is the value from the geometry initialization)
a = 1.5 ;
## read geometry
[x,n,w] = nbisurf("Geometries/sphere-7.nbi") ;
## read solution and extract computed surface potential
f = nbidata("solution.dat") ; pc = f(:,1) ;
## find analytical solution
pa = sphpot(a, x(:,1), x(:,2), x(:,3)) ;
## find error estimate
err = max(abs(pa - pc))/max(abs(pa))


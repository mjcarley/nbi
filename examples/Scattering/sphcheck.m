## check results from sphere scattering calculation

## set wavenumber (this is the value from test-all-helmholtz: change
## it if you change the value used for the test)
k = 10 ;
## read grid with computed scattered field
[x,n,w] = nbisurf("grid.nbi") ;
## find spherical polar coordinates
r = sqrt(sum(x.^2,2)) ;
th = acos(x(:,3)./r) ;
## read field data and convert to complex field amplitude
f = nbidata("grid.dat") ;
f = f(:,1:2)*[1; j] ;
## calculate scattered field for sphere of radius 0.5 and wavenumber k
p = sphscat(0.5, r, th, k) ;
## find points in field lying outside the scattering sphere
ii = find(r > 0.55) ;
## find error estimate for field: maximum difference scaled on maximum
## amplitude
err = max(abs(p(ii) - f(ii)))/max(abs(p(ii)))

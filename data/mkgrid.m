xmin = -2 ; xmax = 3 ; nx = 65 ;
ymin = -2 ; ymax = 3 ; ny = 65 ;
z = 0 ;

[x,y] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny)) ;

dat = [x(:) y(:) z*zeros(size(x(:)))] ;

fid = fopen("grid.dat", "w") ;

fprintf(fid, "%f %f %f\n", dat') ;

fclose(fid) ;


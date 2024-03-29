function [x,n,w]=nbisurf(file)

  ## -- [X, N, W] = nbisurf(FILE)
  ##
  ## Read the points, normals and quadrature weights of an NBI surface
  ##
  ## FILE is the name of an NBI surface, generated using nbi-surface,
  ## or similar
  ##
  ## X: matrix of surface nodes, one per row
  ## N: matrix of surface normals, one per row
  ## W: vector of quadrature weights
  
  nbi_width = 7 ;
  
  fid = fopen(file, "r") ;

  ## clear the ASCII header data
  fscanf(fid, "%c", 40) ;
  ## mesh data
  dat = fscanf(fid, "%d", 2) ;
  nnodes = dat(1) ; npatch = dat(2) ;

  dat = fscanf(fid, "%d", 2*npatch) ;
  p = reshape(dat, 2, npatch)' ;

  dat = fscanf(fid, "%f", nbi_width*nnodes) ;
  dat = reshape(dat, nbi_width, nnodes)' ;
  x = dat(:,1:3) ;
  n = dat(:,4:6) ;
  w = dat(:,  7) ;
  
  fclose(fid) ;

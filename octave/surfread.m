function [x,n]=surfread(file)

  nbi_width = 7 ;
  
  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 2) ;
  nnodes = dat(1) ; npatch = dat(2) ;

  dat = fscanf(fid, "%d", 2*npatch) ;
  p = reshape(dat, 2, npatch)' ;

  dat = fscanf(fid, "%f", nbi_width*nnodes) ;
  dat = reshape(dat, nbi_width, nnodes)' ;
  x = dat(:,1:3) ;
  n = dat(:,4:6) ;
  
  fclose(fid) ;

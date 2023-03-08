function f=datread(file)

  fid = fopen(file, "r") ;

  ## clear the ASCII header data
  fscanf(fid, "%c", 40) ;
  dat = fscanf(fid, "%d", 2) ;
  ## number of data entries per node, number of nodes
  ndat = dat(1) ; nnodes = dat(2) ;

  f = fscanf(fid, "%f", ndat*nnodes) ;

  f = reshape(f, ndat, nnodes)' ;

  fclose(fid) ;

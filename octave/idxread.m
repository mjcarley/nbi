function [idx,ip]=idxread(file, p)

  fid = fopen(file, "r") ;

  np = fscanf(fid, "%d", 1) ;

  ip = fscanf(fid, "%d", 2*(np+1)) ;
  ip = reshape(ip, 2, np+1)' ;

  idx = fscanf(fid, "%d", ip(end,2)+1) ;
  
  fclose(fid) ;

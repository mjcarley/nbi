function S=sptarea(th0,th1,ph0,ph1)

  gm = (th1*ph1 - th0*ph0)/(th1-th0) ;
  a = (ph1 - ph0)/(th1-th0) ;

  S = (sin(gm - a*th1) - sin(gm - a*th0))/a + (th1-th0)*cos(ph0) ;
  

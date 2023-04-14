function p=sphscat(a,r,th,k)

  ## P = SPHSCAT(a, r, th, k)
  ##
  ## Evaluate potential scattered from a sphere under plane wave excitation
  ##
  ## a:  sphere radius
  ## r:  field point radius from sphere centre
  ## th: field point angle from incident wave normal
  ## k:  wavenumber
  ##
  ## If the plane wave is taken to have the form exp(j*k*z), the boundary
  ## is -j*k*exp(j*k*z)*cos(th), with z=a*cos(th) on the sphere surface.
  
  M = 8 ;

  p = 0 ;
  ka = k*a ;
  kr = k*r ;
  for m=0:M
    hmm1 = besselh(m-1/2, 1, ka)*sqrt(pi/2/ka) ;
    hmp1 = besselh(m+3/2, 1, ka)*sqrt(pi/2/ka) ;
    A = (2*m+1)*j^m*real(m*hmm1 - (m+1)*hmp1)/(m*hmm1 - (m+1)*hmp1) ;
    P = legendre(m, cos(th)) ;
    P = P(1,:)' ;

    hm = besselh(m+1/2, 1, kr).*sqrt(pi/2./kr) ;
    p = p - A*P.*hm ;    
  endfor
 

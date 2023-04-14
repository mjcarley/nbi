function p=ring(a, n, k, x, y, z)

  nqp = 1024 ;
  p = zeros(size(x)) ;
  for i=0:nqp-1
    th1 = 2*pi*i/nqp ;
    xs = a*cos(th1) ; ys = a*sin(th1) ;
    R = sqrt((x-xs).^2 + (y-ys).^2 + z.^2) ;

    p += exp(j*(k*R + n*th1))/4/pi./R ;
    
  endfor

  p *= 2*pi/nqp ;

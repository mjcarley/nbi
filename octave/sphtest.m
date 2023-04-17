a = 0.7 ; k = 6.5 ;
z = linspace(-a, a, 65)' ;
th = acos(z/a) ;
x = a*sin(th) ;
r = sqrt(x.^2 + z.^2) ;

ee = 1e-6 ;

zz = (r + ee/2).*cos(th) ;
pp = sphscat(a, r+ee/2, th, k) + exp(j*k*zz) ;
zz = (r - ee/2).*cos(th) ;
pm = sphscat(a, r-ee/2, th, k) + exp(j*k*zz) ;

dp = (pp - pm)/ee ;

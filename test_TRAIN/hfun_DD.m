%
%
%
function res = hfun_DD( x, pars )
  res = 0 ;
  for j=1:2
    dz  = x-pars.zz(j) ;
    dz2 = (dz/pars.epsilon)^2 ;
    res = res + (pars.ss(j+1)-pars.ss(j))*dz/(1+dz2)^2 ;
  end
  res = -2*res / pi / pars.epsilon^3 ;
end

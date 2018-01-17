%
%
%
function res = hfun_D( x, pars )
  res = 0 ;
  for j=1:2
    res = res + (pars.ss(j+1)-pars.ss(j))/(1+((x-pars.zz(j))/pars.epsilon)^2) ;
  end
  res = res / pi / pars.epsilon ;
end

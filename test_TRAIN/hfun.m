%
%
%
function res = hfun( x, pars )
  res = 0 ;
  for j=1:2
    res = res + (pars.ss(j+1)-pars.ss(j))*atan((x-pars.zz(j))/pars.epsilon) ;
  end
  res = res / pi ;
end

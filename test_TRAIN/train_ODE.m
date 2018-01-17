% [ ODE, JAC ]
function ODE = train_ODE( tL, tR, XL, XR, UC, pars )
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  ua = UC(1) ;
  ub = UC(2) ;
  % ----------
  DT = tR - tL ;
  xM = (xR+xL)/2 ;
  vM = (vR+vL)/2 ;
  % ----------
  ODE    = zeros(2,1) ;
  acc    = hfun(xM, pars) - ( pars.alpha + pars.beta * vM + pars.gm * vM^2 ) ; 
  ODE(1) = (xR - xL)/DT - vM ;
  ODE(2) = (vR - vL)/DT - acc - ua + ub ;
end

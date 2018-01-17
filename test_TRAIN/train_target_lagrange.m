% [L,gradL,hessianL]
function L = train_target_lagrange( tL, tR, XL, XR, UC, pars )
  ua = UC(1) ;
  ub = UC(2) ;
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  v  = (vL+vR)/2 ;
  L  = ua * v ;
end
% [L,gradL,hessianL]
function hessL = train_target_lagrange_hessian( tL, tR, XL, XR, UC, pars )
  ua = UC(1) ;
  ub = UC(2) ;
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  v  = (vL+vR)/2 ;
  %gradL = [ 0, ua/2, 0, ua/2,  v, 0] ;
  hessL = [ 0,   0,   0,   0, 0,   0 ; ...
            0,   0,   0,   0, 0.5, 0 ; ...
            0,   0,   0,   0, 0,   0 ; ...
            0,   0,   0,   0, 0.5, 0 ; ...
            0, 0.5,   0, 0.5, 0,   0 ; ...
            0,   0,   0,   0, 0,   0 ] ;
end
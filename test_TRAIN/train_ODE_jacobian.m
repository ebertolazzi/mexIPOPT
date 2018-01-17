% [ ODE, JAC ]
function JAC = train_ODE_jacobian( tL, tR, XL, XR, UC, pars )
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
  tmp_x = -0.5 * hfun_D(xM, pars)  ; 
  tmp_v =  0.5 * pars.beta + pars.gm * vM  ; 
  JAC = [ -1/DT,        -0.5,  1/DT,       -0.5,  0, 0 ; ... 
           tmp_x, tmp_v-1/DT, tmp_x, tmp_v+1/DT, -1, 1 ] ;
end


% [ ODE, JAC ]
function H = train_ODE_hessian( tL, tR, XL, XR, UC, L, pars )
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
  nx = 2 ;
  nu = 2 ;
  tmp_xx = -0.25 * hfun_DD(xM, pars)  ; 
  tmp_vv =  0.5 * pars.gm ; 
  Ahess  = [ tmp_xx, 0 ; 0, tmp_vv ] ;
  H      = L(2) * [ Ahess, Ahess, zeros(nx,nu) ; ...
                    Ahess, Ahess, zeros(nx,nu) ; ...
                    zeros(nu,2*nx+nu) ] ;
end


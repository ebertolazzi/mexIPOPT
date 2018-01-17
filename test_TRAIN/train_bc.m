function bc = train_bc( tL, tR, XL, XR, pars )
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  bc = [ xL - pars.x_i ; ...
         xR - pars.x_f ; ...
         vL - pars.v_i ; ...
         vR - pars.v_f ] ;
end
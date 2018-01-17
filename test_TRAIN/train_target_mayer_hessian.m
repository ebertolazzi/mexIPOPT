% [ M, gradM, hessianM ]
function hessM = train_target_mayer_hessian( tL, tR, XL, XR, pars )
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  hessM = zeros(4,4) ;
end
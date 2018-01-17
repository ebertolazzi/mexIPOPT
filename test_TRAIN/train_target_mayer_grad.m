% [ M, gradM, hessianM ]
function gradM = train_target_mayer_grad( tL, tR, XL, XR, pars )
  xL = XL(1) ;
  vL = XL(2) ;
  xR = XR(1) ;
  vR = XR(2) ;
  gradM = zeros(1,4) ;
end
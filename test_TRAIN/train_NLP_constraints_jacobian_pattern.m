function Jac = train_NLP_constraints_jacobian_pattern( data )
  N    = data.N ;
  nx   = data.nx ;
  nu   = data.nu ;
  nbc  = data.nbc ;
  totx = N*nx ;
  totu = (N-1)*nu ;

  dim = (N-1)*nx+nbc ;
  Jac = sparse( dim, totx+totu  ) ;
  J   = ones(nx,2*nx+nu) ;
  for k=1:N-1
    imap = (k-1)*nx + (1:nx) ;
    jmap = [ (k-1)*nx + (1:2*nx), totx + (k-1)*nu + (1:nu) ] ;
    Jac(imap,jmap) = J ;
  end
  imap = totx - nx + (1:nbc) ;
  jmap = [ 1:nx, totx-nx+(1:nx) ] ;
  Jac(imap,jmap) = ones(nbc,2*nx)  ;
end
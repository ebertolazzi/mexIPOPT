function H = train_NLP_hessian_pattern( data )
  N    = data.N ;
  nx   = data.nx ;
  nu   = data.nu ;
  totx = N*nx ;
  totu = (N-1)*nu ;

  H = sparse( totx + totu, totx + totu ) ;
  imap         = [1:nx, (totx-nx)+(1:nx) ] ;
  H(imap,imap) = ones(2*nx,2*nx);
  for k=1:N-1
    imap = [ (k-1)*nx+(1:2*nx), totx+(k-1)*nu+(1:nu) ] ;
    H(imap,imap) = ones(2*nx+nu,2*nx+nu) ;
  end
  % variation for controls
  for k=2:N-1
    imap = [ totx+(k-2)*nu+(1:2*nu) ] ;
    H(imap,imap) = H(imap,imap) + [ eye(nu), eye(nu) ; eye(nu), eye(nu) ] ;
  end
  H = tril(H) ;
end

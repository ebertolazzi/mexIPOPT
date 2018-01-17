% [ M, gradM, hessianM ]
function H = train_NLP_hessian( Z, sigma, lambda, data )
  nodes = data.nodes ;
  N     = data.N ;
  nx    = data.nx ;
  nu    = data.nu ;
  nbc   = data.nbc ;
  totx  = N*nx ;
  totu  = (N-1)*nu ;
  % divido il vettore in celle di vettori
  X = mat2cell( Z(1:totx), nx*ones(1,N), 1 ) ;
  L = mat2cell( lambda, [nx*ones(1,N-1),nbc], 1 ) ;
  % divido il vettore in celle di vettori
  U = mat2cell( Z(totx+1:end), nu*ones(1,N-1), 1 ) ;
  H = sparse( totx+totu, totx+totu ) ;
  imap         = [1:nx, (totx-nx)+(1:nx) ] ;
  H(imap,imap) = sigma * train_target_mayer_hessian( nodes(1), nodes(end), X{1}, X{end}, data.pars ) + ...
                 train_bc_hessian( nodes(1), nodes(end), X{1}, X{end}, L{end}, data.pars ) ;
  for k=1:N-1
    tmp  = sigma * (nodes(k+1)-nodes(k)) ;
    imap = [ (k-1)*nx+(1:2*nx), totx+(k-1)*nu+(1:nu) ] ;
    H(imap,imap) = H(imap,imap) + ...
                   tmp * train_target_lagrange_hessian( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, data.pars ) + ...
                   train_ODE_hessian( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, L{k}, data.pars ) ;
  end
  % variation for controls
  for k=2:N-1
    imap = [ totx+(k-2)*nu+(1:2*nu) ] ;
    tmp  = 2 * sigma * (nodes(k+1)-nodes(k-1)) * data.pars.uepsi ;
    H(imap,imap) = H(imap,imap) + tmp * [ eye(nu), -eye(nu) ; -eye(nu), eye(nu) ] ;
  end
  H = tril(H) ;
end

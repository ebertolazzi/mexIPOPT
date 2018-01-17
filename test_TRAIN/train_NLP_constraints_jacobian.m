function Jac = train_NLP_constraints_jacobian( Z, data )
  nodes = data.nodes ;
  N    = data.N ;
  nx   = data.nx ;
  nu   = data.nu ;
  nbc  = data.nbc ;
  totx = N*nx ;
  totu = (N-1)*nu ;
  % divido il vettore in celle di vettori
  X = mat2cell( Z(1:totx), nx*ones(1,N), 1 ) ;
  % divido il vettore in celle di vettori
  U = mat2cell( Z(totx+1:end), nu*ones(1,N-1), 1 ) ;
  
  dim = (N-1)*nx+nbc ;
  Jac = sparse( dim, totx+totu ) ;
  for k=1:N-1
    J = train_ODE_jacobian( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, data.pars ) ;
    imap = (k-1)*nx + (1:nx) ;
    jmap = [ (k-1)*nx + (1:2*nx), totx + (k-1)*nu + (1:nu) ] ;
    Jac(imap,jmap) = J ;
  end
  J = train_bc_jacobian( nodes(1), nodes(end), X{1}, X{end}, data.pars ) ;
  imap = totx - nx + (1:nbc) ;
  jmap = [ 1:nx, totx-nx+(1:nx) ] ;
  Jac(imap,jmap) = Jac(imap,jmap) + J ;
end
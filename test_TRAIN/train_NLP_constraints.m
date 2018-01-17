function c = train_NLP_constraints( Z, data )
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
  C = mat2cell( zeros( (N-1)*nx, 1 ), nx*ones(1,N-1), 1 ) ;
  for k=1:N-1
    C{k} = train_ODE( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, data.pars ) ;
  end
  bc = train_bc( nodes(1), nodes(end), X{1}, X{end}, data.pars ) ;
  c  = [ cell2mat( C ) ; bc ] ;
end
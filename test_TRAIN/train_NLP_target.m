% [ M, gradM, hessianM ]
function res = train_NLP_target( Z, data )
  nodes = data.nodes ;
  N  = data.N ;
  nx = data.nx ;
  nu = data.nu ;

  totx = N*nx ;
  totu = (N-1)*nu ;
  % divido il vettore in celle di vettori
  X = mat2cell( Z(1:totx), nx*ones(1,N), 1 ) ;
  % divido il vettore in celle di vettori
  U = mat2cell( Z(totx+1:end), nu*ones(1,N-1), 1 ) ;

  res = train_target_mayer( nodes(1), nodes(end), X{1}, X{end}, data.pars ) ;
  for k=1:N-1
    res = res + (nodes(k+1)-nodes(k)) * train_target_lagrange( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, data.pars ) ;
  end
  % variation for controls
  for k=2:N-1
    res = res + (nodes(k+1)-nodes(k-1)) * data.pars.uepsi * sum( (U{k}-U{k-1}).^2 ) ;
  end
end

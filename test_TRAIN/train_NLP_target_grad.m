% [ M, gradM, hessianM ]
function g = train_NLP_target_grad( Z, data )
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

  gx = mat2cell( zeros( 1, nx*N ),     1, nx*ones(1,N) ) ;
  gu = mat2cell( zeros( 1, nu*(N-1) ), 1, nu*ones(1,N-1) ) ;

  gg      = train_target_mayer_grad( nodes(1), nodes(end), X{1}, X{end}, data.pars ) ;
  ggg     = mat2cell( gg, 1, [ nx, nx ] ) ;
  gx{1}   = ggg{1} ;
  gx{end} = ggg{2} ;
  for k=1:N-1
    gg      = (nodes(k+1)-nodes(k)) * train_target_lagrange_grad( nodes(k), nodes(k+1), X{k}, X{k+1}, U{k}, data.pars ) ;
    ggg     = mat2cell( gg, 1, [ nx, nx, nu ] ) ;
    gx{k}   = gx{k}   + ggg{1} ;
    gx{k+1} = gx{k+1} + ggg{2} ;    
    gu{k}   = ggg{3} ;    
  end
  % variation for controls
  for k=2:N-1
    tmp     = 2 * (nodes(k+1)-nodes(k-1)) * data.pars.uepsi * (U{k}-U{k-1}).' ;
    gu{k}   = gu{k}   + tmp ; 
    gu{k-1} = gu{k-1} - tmp ;
  end
  g = [ cell2mat(gx), cell2mat(gu) ] ;
  % variation for controls
end
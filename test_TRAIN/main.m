addpath('../lib');

x_i   = 0;
x_f   = 6;
v_i   = 0;
v_f   = 0;
alpha = 0.3;
beta  = 0.14;
gm    = 0.16;
uepsi = 1e-7 ;

nx    = 2 ;
nu    = 2 ;
nbc   = 4 ;
N     = 1000 ;
nodes = 4.8*(0:N-1)/N ;

auxdata.pars.x_i   = x_i ;
auxdata.pars.x_f   = x_f ;
auxdata.pars.v_i   = v_i ;
auxdata.pars.v_f   = v_f ;
auxdata.pars.alpha = alpha ;
auxdata.pars.beta  = beta ;
auxdata.pars.gm    = gm ;
auxdata.pars.uepsi = uepsi ;

auxdata.pars.epsilon = 0.05 ;
auxdata.pars.ss      = [ -2, 0, 2 ];
auxdata.pars.zz      = [ 2, 4 ];

auxdata.nodes      = nodes;

auxdata.nx         = nx ;
auxdata.nu         = nu ;
auxdata.nbc        = nbc ;
auxdata.N          = N ;

uaMax = 10 ;
ubMax = 2 ;

options.lb = [-ones(1,N*nx)*Inf, zeros( 1, (N-1)*nu ) ] ;  % Lower bound on the variables.
options.ub = [ ones(1,N*nx)*Inf, reshape( [uaMax;ubMax]*ones(1,N-1), 1, nu*(N-1)) ] ;  % Upper bound on the variables.

% The constraint functions are bounded to zero
options.cl = zeros(1,(N-1)*nx+nbc); %  constraints
options.cu = zeros(1,(N-1)*nx+nbc);

% Set up the auxiliary data.
options.auxdata = auxdata ;

% Set the IPOPT options.
options.ipopt.jac_d_constant   = 'no';
options.ipopt.hessian_constant = 'no';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = 400;
options.ipopt.tol              = 1e-10;%
%options.ipopt.linear_solver    = 'ma57';
%options.ipopt.linear_solver    = 'mumps';
%options.ipopt.linear_solver    = 'pardiso';
options.ipopt.pardiso_msglvl   = 4 ;

% The callback functions.
funcs.objective         = @train_NLP_target;
funcs.gradient          = @train_NLP_target_grad;

funcs.constraints       = @train_NLP_constraints;
funcs.jacobian          = @train_NLP_constraints_jacobian;
funcs.jacobianstructure = @train_NLP_constraints_jacobian_pattern;

if true
  %options.ipopt.derivative_test = 'second-order';
  funcs.hessian           = @train_NLP_hessian;
  funcs.hessianstructure  = @train_NLP_hessian_pattern;
else
  options.ipopt.derivative_test = 'first-order';
  options.ipopt.hessian_approximation      = 'limited-memory';
  %options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
  %options.ipopt.limited_memory_update_type = 'sr1' ;
  options.ipopt.limited_memory_update_type = 'bfgs' ; % {bfgs}, sr1 = 6; % {6}
end

% Run IPOPT.
xguess  = x_i+(x_f-x_i)*nodes./4.8 ;
vguess  = ones(1,N) ;
uaguess = zeros(1,N-1) ;
ubguess = zeros(1,N-1) ;

x0 = [ reshape( [ xguess,   vguess], 2*N ,1 ) ; ...
       reshape( [ uaguess, ubguess], 2*(N-1) ,1 ) ] ;

tic
[sol, info] = ipopt_auxdata(x0,funcs,options);
elapsed = toc ;

x  = sol(1:2:2*N) ;
v  = sol(2:2:2*N) ;
ua = sol(2*N+(1:2:2*N-2)) ;
ub = sol(2*N+(2:2:2*N-2)) ;

subplot( 3, 1, 1 );  
plot( nodes, x, 'Linewidth', 2 ) ;

subplot( 3, 1, 2 );  
plot( nodes, v, 'Linewidth', 2 ) ;

subplot( 3, 1, 3 );  
plot( nodes(1:end-1), ua, nodes(1:end-1), ub, 'Linewidth', 2 ) ;

info;

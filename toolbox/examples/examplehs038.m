% Test the "ipopt" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.
%
% Copyright (C) 2008 Peter Carbonetto. All Rights Reserved.
% This code is published under the Eclipse Public License.
%
% Author: Peter Carbonetto
%         Dept. of Computer Science
%         University of British Columbia
%         September 18, 2008
function [x, info] = examplehs038

  addpath('../lib');

  x0         = [-3  -1  -3  -1];   % The starting point.
  options.lb = [-10 -10 -10 -10];  % Lower bound on the variables.
  options.ub = [+10 +10 +10 +10];  % Upper bound on the variables.

  % The callback functions.
  funcs.objective = @objective;
  funcs.gradient  = @gradient;
  funcs.iterfunc  = @callback;

  if true
    funcs.hessian                 = @hessian;
    funcs.hessianstructure        = @hessianstructure;
    options.ipopt.derivative_test = 'second-order';
  else
    options.ipopt.hessian_approximation      = 'limited-memory';
    options.ipopt.limited_memory_update_type = 'bfgs';
    options.ipopt.derivative_test            = 'first-order';
  end

  % Set the IPOPT options.

  options.ipopt.linear_solver = 'mumps';
  % HSL solver family
  % to use this solvers see README_HSL.md
  %options.ipopt.linear_solver    = 'ma57';
  %options.ipopt.linear_solver    = 'ma77';
  %options.ipopt.linear_solver    = 'ma86';
  %options.ipopt.linear_solver    = 'ma97';

  % PARDISO solver
  % to use this solvers see README_HSL.md
  %options.ipopt.linear_solver    = 'pardiso';
  %options.ipopt.pardiso_msglvl   = 4;

  options.ipopt.mu_strategy   = 'adaptive';
  options.ipopt.print_level   = 0;
  options.ipopt.tol           = 1e-7;
  options.ipopt.max_iter      = 100;

  % Run IPOPT.
  [x info] = ipopt(x0,funcs,options);
end

% ----------------------------------------------------------------------
function f = objective(x)
  f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);
end

% ----------------------------------------------------------------------
function g = gradient(x)
  g(1) = -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
  g(2) = 200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
  g(3) = -360*x(3)*(x(4)-x(3)^2) -2*(1-x(3));
  g(4) = 180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1);
end

% ----------------------------------------------------------------------
function H = hessianstructure()
  H = sparse([ 1  0  0  0
               1  1  0  0
               0  0  1  0
               0  1  1  1 ]);
end

% ----------------------------------------------------------------------
function H = hessian(x, sigma, lambda)
  H = [ 1200*x(1)^2-400*x(2)+2  0       0                          0
        -400*x(1)               220.2   0                          0
         0                      0       1080*x(3)^2- 360*x(4) + 2  0
         0                      19.8   -360*x(3)                   200.2 ];
  H = sparse(sigma*H);
end

% ----------------------------------------------------------------------
function b = callback(varargin) %t, f, x)
  varargin{1}
  varargin{2}
  varargin{3}
  %fprintf('%3d  %0.3g \n',t,f);
  b = true;
end

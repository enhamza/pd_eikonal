%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Eikonal equation
% |Du(x)| = f(x) in Omega anf u = g on the boundary
% using the variational formulation proposed in 
% H.Ennaji, N.Igbida and V.T.Nguyen, Augmented Lagrangian methods for
% degenerate Hamilton-Jacobi equations, Calc. Var. Partial Differential
% Equations 2021
% Here  Chambolle-Pock's PD algorithm is used instead of ALG2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all; close all;
clc;

addpath('toolbox/');
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Omega = [0,L1]x[0,L2]
a = 0; b=1;c = a; d = b;
h=0.01;
nx = floor((b-a)/h+1);
ny = nx;
x = linspace(a,b,nx);
y = linspace(c,d,ny);
[X,Y] = meshgrid(x,y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%operators
K = @(u)  mygrad(u,h);
Ks = @(phi) -mydiv(phi(:,:,1),phi(:,:,2),h);
Prox_tau_F=@(v,tau)  v + tau;
%intialization
u0 = zeros(nx,ny);
phi = zeros(nx,ny,2);
u = u0;
uold =u0;
ubar = u;
g = zeros(nx,ny);%Boundary data
f = ones(nx,ny);%Rhs for Eikonal equation
%f =  2*pi*sqrt(cos(2*pi*X).^2*sin(2*pi*Y).^2 + cos(2*pi*Y).^2*sin(2*pi*X).^2);
%f = exp(-(X.^2+Y.^2)/0.3);
%uncomment to compare with the solution of |Du| = 1
%uexact = min(1-X,min(1-Y,min(X,Y)));%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of PD algorithm
tau = h/3;
sigma = h/3;
pd_itermax = 1e5;
tol = 1e-6;%tolerance used to stop iterations
E = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%index to subscript and vice versa
%L = [x0 x1 x2 x3 x4 x5 x6 x7];
%BC:    use set_bdcond(u,g) for Dirichelt BC, and full syntax
%       set_bdcond(u,g,L,0) for u(L) = 0,where L is a list of data points
%       (possibly a singleton)
bdidx = set_bdcond(u,g);



    %main loop
    tic;%time
    for k=1:pd_itermax
                 %

                %step1: Dual Step
                [d1,d2] = mygrad(ubar,h);
                phiaux  = phi + sigma*cat(3,d1,d2);
                phiaux = reshape(phiaux,nx*ny,2);
                p = pj2(phiaux./sigma,f(:));
                phi = phiaux - sigma.*p;


                 %step2: Primal step
                phi = reshape(phi,[nx,ny,2]);
                v = uold - tau*Ks(phi);
                u(bdidx) = Prox_tau_F( v(bdidx),tau);
                
                %error check
                E(k) = norm(u-uold)/norm(u);
                if(E(k)<tol)
                   break;
                end                 
                %update:
                ubar =2*u-uold;
                uold = u;
                %iterations disp (comment)
                  if mod(k,100) == 0 
                   fprintf('iter = %4d, \n', k);
                  end

                    
    end
    toc;

%solution display
figure;
mesh(X,Y,u); 
figure;
s = plot(E, '-r'); set(s, 'LineWidth', 2);
xlabel('Iterations');
ylabel('$\frac{\Vert u_{k} - u_{k-1}\Vert}{\Vert u_{k-1}\Vert}$','Interpreter','latex','FontSize',15);
set(gca,'yscale','log')
legend('PD');
figure;
contourf(u,40);

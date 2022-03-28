%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Eikonal equation
% F^*(x,Du(x)) = 1 in Omega anf u = g on the boundary
% where F is a crystalline norm
% using Chambolle-Pock's PD algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all; close all;
clc;

addpath('../toolbox/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Omega = [0,L1]x[0,L2]
a = -1; b=1;c = a; d = b;
h=0.02;
nx = floor((b-a)/h+1);
ny = nx;
x = linspace(a,b,nx);
y = linspace(c,d,ny);
[X,Y] = meshgrid(x,y);

%examples:
%s=[[1, -1];[1, 1];[-1, 1];[-1, -1];[1, -1];[1, 1]];
%vv=[[1, 0];[0, 1]; [-1,0]; [0, -1]]; 
%ss=[[0, 2];[-2, 0]; [0, -2]; [2, 0]; [0, 2]];
s= [[1 -1]; [1,-0.8]; [-0.8,1]; [-1,1]; [-1,-1]; [1,-1]; [1,-0.8]];
vv=[[1, 0]; [5, 5]; [0,1]; [-1, 0]; [0, -1]];
ss=[[0,0.2];[-1.8, 1.8]; [-0.2, 0]; [0, -2]; [2, 0]; [0,0.2]];
% s=[[1, 0]; 	[0.2, 0.2]; 	[-0.2, 0.2]; 	[-0.2, -0.2]; 	[0.2, -0.2]; 	[1, 0]; 	[0.2, 0.2]];
% ss=[[-0.8, 0.2];[-0.4, 0]; [0, -0.4]; [0.4, 0]; [0.8, 0.2]; [-0.8,0.2]];
% vv=[[1, 4]; [0, 5]; [-5,0]; [0, -5]; [1, -4]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%operators
K = @(u)  mygrad(u,h);
Ks = @(phi) -mydiv(phi(:,:,1),phi(:,:,2),h);
Prox_tau_F=@(v,tau)  v + tau;
%intialization
u0 = zeros(nx,ny);
phi = K(u0);%zeros(nx,ny,2);
u = u0;
uold =u0;
ubar = u;
g = zeros(nx,ny);%Boundary data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters of PD algorithm
tau = h/3;
sigma = h/3;
pd_itermax = 3000;
tol = 1e-5;%tolerance used to stop iterations
E = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%index to subscript and vice versa
%L = [x0 x1 x2 x3 x4 x5 x6 x7];
%BC:    use set_bdcond(u,g) for Dirichelt BC, and full syntax
%       set_bdcond(u,g,L,0) for u(L) = 0,where L is a list of data points
%       (possibly a singleton)
 mm = floor(size(u,1)/2);
x0 = [mm;mm];
bdidx = set_bdcond(u,g,x0,0);



    %main loop
    tic;%time
    for k=1:pd_itermax
                 %

                %step1: Dual Step
                [d1,d2] = mygrad(ubar,h);
                phiaux  = phi + sigma*cat(3,d1,d2);
                   for i=1:nx
                        for j=1:ny
                            V = [phiaux(i,j,1)/sigma;phiaux(i,j,2)/sigma];
                             P = Proj_Finsler_ball(V,vv,s,ss);
                             phi(i,j,1) = phiaux(i,j,1) - sigma * P(1);
                             phi(i,j,2) = phiaux(i,j,2) - sigma * P(2);
                 
                        
                       end
                    end
    


                 %step2: Primal step
                phi = reshape(phi,[nx,ny,2]);
                v = uold - tau*Ks(phi);
                u(bdidx) = Prox_tau_F( v(bdidx),tau);
                
%                 %error check
%                 E(k) = norm(u-uold)/norm(u);
%                 if(E(k)<tol)
%                    break;
%                 end                 
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
contourf(u,40);
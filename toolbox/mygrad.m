function [Dx,Dy] = mygrad(u,h)
  [ny,nx] = size(u);
   if nargin <2
       h = 1.;
   end
  Dx = (u(:,[2:nx,nx]) - u)/h;
  Dy = (u([2:ny,ny],:) - u)/h;
end
function d = mydiv(px,py,h) 
%// compute d = div(p) where p = (px,py) and size(px) == size(py)
  [ny,nx] = size(px);
   if nargin <3
       h = 1.;
   end
  div_x      =  (px - px(:,[1,1:(nx-1)]))/h;
  div_x(:,1) =  px(:,1)/h;
  div_x(:,nx) = -px(:,nx-1)/h;
  div_y      =  (py - py([1,1:(ny-1)],:))/h;
  div_y(1,:) =  py(1,:)/h;
  div_y(ny,:) = -py(ny-1,:)/h;
  d = div_x + div_y;
end
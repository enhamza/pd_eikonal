function [x] = pj2(mat,rad)
if nargin<2
    rad = 1;
end
nrm = sqrt(mat(:,1).^2 + mat(:,2).^2);
inside = nrm < rad(:);
x = rad .* mat./[nrm,nrm];
x(inside,:) = mat(inside,:);
        
end
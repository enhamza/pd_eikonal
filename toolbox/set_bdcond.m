%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function sets BC and gets the indices outside the boundary region:
% Inputs: u matrix (for instance n*n, general case later ?)
%         g boundary condition on ∂Ω (also n*n)
%         x0 point at which val is prescribed i.e.
%         u(x0) = val
% Output bdix: indices that are not affected by BC (to be applied for
%              prox_F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bdidx = set_bdcond(u,g,x0,val)

    m = size(u,1);n = size(u,2);
    if n~=m
        mm = sqrt(max(m,n));
        u = reshape(u,mm,mm);
    end
        sub2ind1 = @(x)(x(2)-1)*m+x(1);

    % BC on ∂Ω
    if nargin<3
        u(:,1) = g(:,1);
        u(:,end) = g(:,end);
        u(end,:) = g(end,:);
        u(1,:) = g(1,:); 
        isbd = true(size(u));
        isbd(:,1) = false;isbd(1,:) = false;
        isbd(end,:) = false;isbd(:,end) = false;
        %bdidx = find(isbd(:));
        bdidx = isbd;
   
    else %BC on x0 (or elements of x0)
        isbd = true(size(u));
        I = sub2ind(size(u), x0(1,:), x0(2,:));
        for k=1:size(I,2)
        u(I(k)) = val;
        isbd(I(k)) = false;
        
        %bdidx = find(isbd(:));
        bdidx = isbd;
        end
    end

end


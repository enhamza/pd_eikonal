function [x1,x2] = pj(y1,y2,rad)
nrm = sqrt(sum(sum(y1.^2 + y2.^2)));

    %if nrm <= rad, x1 = y1; x2 = y2;
    %else
        x1 = y1 .* (rad/max(nrm,1));
        x2 = y2 .* (rad/max(nrm,1));
    %end
        
end
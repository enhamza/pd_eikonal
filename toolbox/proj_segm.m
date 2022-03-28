function [p] = proj_segm(xx,v1,s1)
r = (scalar(xx,v1) - scalar(v1,s1))/scalar(v1,v1);

p = [xx(1) - r*v1(1), xx(2) - r*v1(2)];
end
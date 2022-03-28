function [rr1] = Proj_Finsler_ball(rr,v,s,ss)
rr1 = [ 0,0];
for j=1:size(v,1)

	if(scalar(rr, v(j,:)) >= scalar(s(j,:),v(j,:)))
            if (scalar(rr, ss(j,:))>= scalar(s(j,:),ss(j,:)))
                if (scalar(rr, ss(j,:))<=scalar(s(j+1,:), ss(j,:))) 
                rr1=proj_segm(rr,v(j,:),s(j,:));
                break;
                end
            end
      end
 
    

    if (scalar(rr, ss(j,:))>=scalar(s(j+1,:),ss(j,:))) 
            if (scalar(rr, ss(j+1,:))<=scalar(s(j+1,:),ss(j+1,:))) 
            rr1=s(j+1,:);
            break;
            end
    end
end

end



function out= E_MAX(E1,E2)
% function for calculating maximum envelope amplitude regardless of direction

out = zeros(length(E1),1);
for i =1:length(E1)
    if norm(E1(i,:))>norm(E2(i,:))
        e_max = E1(i,:);
        e_min = E2(i,:);
    else
        e_max = E2(i,:);
        e_min = E1(i,:);
    end
    if norm(e_max-e_min)>norm(e_max+e_min)
        e_min = -e_min;
    end
    
    if norm(e_min)^2 <= abs(e_max*e_min')
        out(i) = 2*norm(e_min);
    else
        out(i) = 2*sqrt(sum(cross(e_min,e_max-e_min).^2,2))/sqrt(sum((e_max-e_min).^2,2));
    end
end

end
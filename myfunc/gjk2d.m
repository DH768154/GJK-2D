function [dist,closestpt] = gjk2d(O,R)

simplex = NaN(2,3);

maxloop = 100;
indx = [1,3;2,3;1,3];
dist = nan;
closestpt = NaN(2,1);

%% Construct first 2 Vertex of the Simplex 

for n = 1:2
    if n == 1
        dirc = mean(O,2) - mean(R,2); % 1st dirction
    elseif n == 2
        dirc = -simplex(:,1); % 2nd dirction
    end
    simplex(:,n) = supportfunc(O,dirc) + supportfunc(-R,dirc);
end
dirc = axbxb2(simplex(:,1),simplex(:,2)-simplex(:,1));

%% 
for i = 1:maxloop

    %% 3rd Vertex of the Simplex
    simplex(:,3) = supportfunc(O,dirc) + supportfunc(-R,dirc);

    %% No more Better Points
    if any(all(simplex(:,3) == simplex(:,1:2)))
        closestpt = proj2line(simplex(:,1:2),[0;0],'bound');
        dist = norm(closestpt);
        return
    end

    %% Update Simplex

    [closest_pt,type,edge_ind] = proj2simplex2(simplex,[0;0]);
    dirc = -closest_pt;


%%

    % In the Triangle, Collision
    if type==1 
        closestpt = NaN(2,1);
        dist = 0; return;

    % Closest to Vertex
    elseif type==3 % 
        simplex(:,1:2) = simplex(:,indx(edge_ind(1),:));
    
    % Closest to Edge
    else 
        simplex(:,1:2) = simplex(:,edge_ind);
    end

end

end

function n = axbxb2(a,b)
% first 2 element of cross(cross([a;0],[b;0]),[b,0]) 
n = [-b(2,:);b(1,:)] .* (a(1,:).*b(2,:)-a(2,:).*b(1,:));
end
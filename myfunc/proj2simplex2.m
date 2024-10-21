function [closest_pt,type,edge_ind] = proj2simplex2(vertices,pt)
%
% type:
% 1: in simplex
% 2: project to edge
% 3: project to vertex
%%
closest_pt = NaN(size(pt));
type = deal(NaN);
edge_ind = NaN(1,2);
%%

[n,edgeV] = norm2edgev(vertices); % Normal Vector of each Edge

ABCP = pt - vertices; % AP, BP, CP
outind = sum(n.*ABCP);  % [dot(n(1),AP), dot(n(2),BP), dot(n(3),CP)]

AP_AB =  dot(ABCP(:,1),edgeV(:,1)); % dot(P-A,B-A) === A →
BP_BA = -dot(ABCP(:,2),edgeV(:,1)); % dot(P-B,B-A) === ← B
BP_BC =  dot(ABCP(:,2),edgeV(:,2)); % dot(P-B,C-B) === B →
CP_CB = -dot(ABCP(:,3),edgeV(:,2)); % dot(P-C,B-C) === ← C
CP_CA =  dot(ABCP(:,3),edgeV(:,3)); % dot(P-C,A-C) === C →
AP_AC = -dot(ABCP(:,1),edgeV(:,3)); % dot(P-A,C-A) === ← A

% inside Triangle
if all(outind<=0) 

    type = 1;
    edge_ind = [0,0];
    if size(pt,1) == 2
        closest_pt = pt;
    else % 3D case, project to the plane
        N = cross(edgeV(:,1),edgeV(:,2));
        D = -dot(N,vertices(:,1));
        dist = (dot(N, pt) + D) / norm(N);
        closest_pt = pt - dist * (N / norm(N));
    end

% BC Side, B → ← C === BC Edge
elseif all([outind(2),BP_BC,CP_CB]>0) 
    edge_ind = [2,3]; 
    type = 2;
    closest_pt = BP_BC/sum(edgeV(:,2).^2)*edgeV(:,2)+vertices(:,2);   

% AC Side, C → ← A === AC Edge    
elseif all([outind(3),CP_CA,AP_AC]>0) 
    edge_ind = [1,3];
    type = 2;
    closest_pt = CP_CA/sum(edgeV(:,3).^2)*edgeV(:,3)+vertices(:,3);

% B → ← B === Vertex B    
elseif all([BP_BA,BP_BC]<=0) 
    edge_ind = [2,2];
    type = 3;
    closest_pt = vertices(:,2);   

% C → ← C === Vertex C    
elseif all([CP_CB,CP_CA]<=0) 
    edge_ind = [3,3];
    type = 3;
    closest_pt = vertices(:,3);
   
% Following these will not be touch if use gjk, so check last

% A → ← B === AB Edge 
elseif all([outind(1),AP_AB,BP_BA]>0) 
    edge_ind = [1,2]; 
    type = 2;  
    closest_pt = AP_AB/sum(edgeV(:,1).^2)*edgeV(:,1)+vertices(:,1);   

% A → ← A === Vertex A    
elseif all([AP_AB,AP_AC]<=0) 
    edge_ind = [1,1];
    type = 3;
    closest_pt = vertices(:,1);

end

end
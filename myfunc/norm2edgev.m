function [n,edgeV] = norm2edgev(vertices,varargin)
% edgeV: edge vector of a shape, ccw
% n: Normal vector point out

shapeloop = [vertices,vertices(:,1)];
%edgeV = diff(shapeloop,1,2); % AB, BC, CA
edgeV = shapeloop(:,2:end) - shapeloop(:,1:end-1);

% cross(cross(a,b),b)
b = [edgeV(:,2:end),edgeV(:,1)];
if size(vertices,1)==3
    n = -cross(cross(edgeV,b),b);
else
    n = - [-b(2,:);b(1,:)] .* (edgeV(1,:).*b(2,:)-edgeV(2,:).*b(1,:));
end
n = [n(:,end),n(:,1:end-1)];

if ~isempty(varargin)
    if strcmpi(varargin{1},'norm')
        n = n ./ sqrt(sum(n.^2));
    end
end

end
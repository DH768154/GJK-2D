function s_pt = supportfunc(pts, dirc)
% support point on dirction
% max normal distance from the points to the normal vector of the dirction
% vector
d = pts'*dirc;
[~,ind] = max(d);
s_pt = pts(:,ind);
end
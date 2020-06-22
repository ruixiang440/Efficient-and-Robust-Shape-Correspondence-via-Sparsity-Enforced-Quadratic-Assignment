function [StiffM,MassM] = Stiff(pt,trg)

num_pt = size(pt,1);
N = num_pt;

%------------------------------
% Adjacency matrgx
% The incidence matrgx in the old code is the adjaency matrgx

Adj = sparse(trg(:,[1 2 3]),trg(:,[2 3 1]),1,num_pt,num_pt);
Adj = double(Adj|Adj'); % adjacent nodes
[i,j] = find(Adj);  %%% row and col index of nonzeros in Adj, [i(k) j(k)] is an edge
[ei,ej] = find(Adj(:,i) & Adj(:,j)); % common adjacent trgangle-edge
e1 = find([1; diff(ej)]); % 1:2:end  % ->ignors singularities
%e2 = e1+1;                % 2:2:end
e2 = find([diff(ej);1]);
ss = ones(size(e1))/2;
ss(find(e1 - e2)) = 1;
%------------------------------
% ptiantes of adjacent nodes and of common adjacent trgangle-edge
pi = pt(i,:);
pj = pt(j,:);
qi = pt(ei(e1),:); % qi = pt(ei(1:2:end),:);
qj = pt(ei(e2),:); % qj = pt(ei(2:2:end),:);

% distances
dii = pi-qi; dij = pi-qj;
dji = pj-qi; djj = pj-qj;
eij = pj-pi;


norm1 = @(x) sqrt(sum(x.^2,2)); % norm of a vector

%------------------------------
% % w = exp(2*u);   %%%%% conformal factor
% % Wqi = (w(ei(e1)) + w(i) + w(j))/3;
% % Wqj = (w(ei(e2)) + w(i) + w(j))/3;


% % Mass Matrix
sarea = @(vi,vj) norm1(cross(vi,vj))./2;
%sarea = @(vi,vj) abs(vi(:,1).*vj(:,2) - vi(:,2).*vj(:,1))./2;
MassM = sparse(i,j,ss.*(sarea(dii,dji) + sarea(dij,djj))./12,num_pt,num_pt);
MassM = MassM + diag(sum(MassM,2));



%------------------------------
% Stiffness matrgx
cotan = @(vi,vj) cot( acos(dot(vi,vj,2)./(norm1(vi).*norm1(vj))) );
StiffM = sparse(i,j,-ss.*(cotan(dii,dji) + cotan(dij,djj))./2,num_pt,num_pt);
StiffM = StiffM - diag(sum(StiffM,2));
end


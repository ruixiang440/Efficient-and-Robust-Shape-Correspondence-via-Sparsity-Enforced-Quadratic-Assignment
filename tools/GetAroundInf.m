function [aroundpt,aroundtrg, trgaroundtrg, EdgeTrg, EdgeLength] = GetAroundInf(Surf)
%% Rongjie Lai

if nargin < 1
    disp('Parameter:  aroundtrg = getaroundtrg(Surf)');
    return
end
pt = Surf.pt;
trg = Surf.trg;
num_pt = size(pt,1);

aroundtrg=cell(length(pt),1);
aroundpt=cell(length(pt),1);
EdgeTrg = sparse(num_pt,num_pt);
Edge = sparse(num_pt,num_pt);
EdgeLength = sparse(num_pt,num_pt);
for i=1:length(trg)
    p1=trg(i,1);
    p2=trg(i,2);
    p3=trg(i,3);
    v1 = pt(p1,:);
    v2 = pt(p2,:);
    v3 = pt(p3,:);
    aroundtrg{p1}=union(aroundtrg{p1},i);
    aroundtrg{p2}=union(aroundtrg{p2},i);
    aroundtrg{p3}=union(aroundtrg{p3},i);
    
    aroundpt{p1}=union(aroundpt{p1},[p1 p2 p3]);
    aroundpt{p2}=union(aroundpt{p2},[p1 p2 p3]);
    aroundpt{p3}=union(aroundpt{p3},[p1 p2 p3]);
    
    EdgeLength(p1,p2) = norm(v1 - v2);   EdgeLength(p2,p1) = EdgeLength(p1,p2);
    EdgeLength(p2,p3) = norm(v2 - v3);   EdgeLength(p3,p2) = EdgeLength(p2,p3);
    EdgeLength(p1,p3) = norm(v1 - v3);   EdgeLength(p3,p1) = EdgeLength(p1,p3);
    if Edge(p1,p2) ~= 1
        EdgeTrg(p1,p2) = i;
        Edge(p1,p2) = 1;
    else
        EdgeTrg(p2,p1) = i;
        Edge(p2,p1) = 1;
    end
    if Edge(p1,p2) ~= 1
        EdgeTrg(p1,p3) = i;
        Edge(p1,p3) = 1;
    else
        EdgeTrg(p3,p1) = i;
        Edge(p3,p1) = 1;
    end
    if Edge(p2,p3) ~= 1
        EdgeTrg(p2,p3) = i;
        Edge(p2,p3) = 1;
    else
        EdgeTrg(p3,p2) = i;
        Edge(p3,p2) = 1;
    end
end


% trgaroundtrg=zeros(size(trg));
% for i=1:length(trg)
%     p1=trg(i,1);
%     p2=trg(i,2);
%     p3=trg(i,3);
%     temptrg1=intersect(aroundtrg{p2},aroundtrg{p3});
%     temptrg2=intersect(aroundtrg{p3},aroundtrg{p1});
%     temptrg3=intersect(aroundtrg{p1},aroundtrg{p2});
%     temptrg1=setdiff(temptrg1,i);
%     temptrg2=setdiff(temptrg2,i);
%     temptrg3=setdiff(temptrg3,i);
%     if isempty(temptrg1)
%         temptrg1=0;
%     end
%     if isempty(temptrg2)
%         temptrg2=0;
%     end
%     if isempty(temptrg3)
%         temptrg3=0;
%     end
%     trgaroundtrg(i,:)=[temptrg1,temptrg2,temptrg3];
%     trgcenter(i,:)=(pt(p1,:)+pt(p2,:)+pt(p3,:))/3;
% end
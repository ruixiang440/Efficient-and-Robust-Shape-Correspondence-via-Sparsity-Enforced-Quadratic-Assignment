function D = FastMarchingSR(initialpt,Surf)
  if nargin<1
      disp('Parameter: D=FastMarch(p0,fname)')
      return;
  end
p0=initialpt;
pt = Surf.pt;
trg = Surf.trg;
aroundpt = Surf.aroundpt;
ptaroundtrg = Surf.aroundtrg;
% for i=1:num_pt
%     aroundpt{i}=setdiff(aroundpt{i},i);
% end
num_pt=size(aroundpt{initialpt},1);
% num_pt = n;
%########### Initializatoin ###############
X=[1:num_pt];

num=size(pt,1);
D=10^10*ones(num,1);
D(p0)=0;
Black=[p0];
Red=p0;
for i=1:length(p0)
    Red=union(Red,aroundpt{p0(i)});
end
% Green=setdiff(X,union(Black,Red));
% for i=1:size(Red)
%     D(Red(i))=norm(pt(p0,:)-pt(Red(i),:));
% end
%############  Iteration
visit = [];
while (length(Black)<num_pt)
    [~,ind]=min(D(Red));
    x1=Red(ind(1));
    trgaroundx1=ptaroundtrg{x1};

    for i=1:length(trgaroundx1)
        temptrg=trg(trgaroundx1(i),:);
        x=setdiff(temptrg,x1);
        x2=x(1); x3=x(2);
%         if ismember(x(2),Black) || ismember(x3,visit)
%             x2=x(2);
%             x3=x(1);
%         end
        if ~ismember(x3,Black) && ~ismember(x3,visit)
            %############## Update d3 on the triangle [x1 x2 x3]
            v1=pt(x1,:);
            v2=pt(x2,:);
            v3=pt(x3,:);
            V=[(v1-v3)',(v2-v3)'];
            Q=(V'*V)^(-1);
            d=[D(x1);D(x2)];
            I=[1;1];
            p=(I'*Q*d+sqrt((I'*Q*d)^2-I'*Q*I*(d'*Q*d-1)))/(I'*Q*I);
            tempvect=Q*(d-p*I);
            if (tempvect(1)<0 && tempvect(2)<0)
                D(x3)=min([D(x3),p]);
            else D(x3)=min([D(x3),D(x1)+ norm(v1-v3),D(x2)+ norm(v2-v3)]);
            end
            %##################################################
%             Red=union(Red,x3);
            visit = union(visit,x3);
        end
         x2=x(2);
         x3=x(1);
         if ~ismember(x3,Black) && ~ismember(x3,visit)
            %############## Update d3 on the triangle [x1 x2 x3]
            v1=pt(x1,:);
            v2=pt(x2,:);
            v3=pt(x3,:);
            V=[(v1-v3)',(v2-v3)'];
            Q=(V'*V)^(-1);
            d=[D(x1);D(x2)];
            I=[1;1];
            p=(I'*Q*d+sqrt((I'*Q*d)^2-I'*Q*I*(d'*Q*d-1)))/(I'*Q*I);
            tempvect=Q*(d-p*I);
            if (tempvect(1)<0 && tempvect(2)<0)
                D(x3)=min([D(x3),p]);
            else D(x3)=min([D(x3),D(x1)+ norm(v1-v3),D(x2)+ norm(v2-v3)]);
            end
            %##################################################
%             Red=union(Red,x3);
            visit = union(visit,x3);
         end
    end
    Red=setdiff(Red,x1);
    Black=union(Black,x1);
end
D(D==10^10)=0;
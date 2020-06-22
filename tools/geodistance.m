function [d] = geodistance(a,b,Surf)

if a == b
    d = 0;
else    
        p0=a;
        pt = Surf.pt;
        trg = Surf.trg;
        aroundpt = Surf.aroundpt;
        ptaroundtrg = Surf.aroundtrg;
        num=size(pt,1);
        D=10^10*ones(num,1);
        D(p0)=0;
        Black=[p0];
        Red=[];
        for i=1:length(p0)
            Red=union(Red,aroundpt{p0(i)});
        end

        while (D(b)==10^10)
            [~,ind]=min(D(Red));
            x1=Red(ind(1));
            trgaroundx1=ptaroundtrg{x1};

            for i=1:length(trgaroundx1)
                temptrg=trg(trgaroundx1(i),:);
                x=setdiff(temptrg,x1);
                x2=x(1); x3=x(2);
                if ismember(x(2),Black)
                    x2=x(2);
                    x3=x(1);
                end
                if ~ismember(x3,Black)
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
                    Red=union(Red,x3);
                end
            end
            Red=setdiff(Red,x1);
            Black=union(Black,x1);
        end
        d=D(b);     
end
end





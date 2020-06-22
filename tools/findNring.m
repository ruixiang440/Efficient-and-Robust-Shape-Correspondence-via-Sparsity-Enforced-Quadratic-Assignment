function [Nring] = findNring(idx, surf, N)

Nring = idx;
% dif = 0;
for i=1:N
%     t = length(Nring);
    for j=1:length(Nring)
        Nring = [Nring;surf.aroundpt{Nring(j)}];
    end
    Nring = unique(Nring);

end
% Nring = setdiff(Nring,idx);

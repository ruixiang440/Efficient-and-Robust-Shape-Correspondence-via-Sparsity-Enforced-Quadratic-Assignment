function [ M] = PerMatrix(P)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% n=length(P);
% M=zeros(n,n);
% for i=1:n;
%         M(i,P(i))=1;
% end
num = length(P);
M = sparse(1:num,P,1,num,num);

end


function [out]  =  Proj_DoubleStochasticM_reloc(Y)
global  bw 
% global sparsity
% Y = full(Y);
% ii = [];
% jj = [];
% ss = [];
% for i=1:length(sparsity)  
%     ii = [ii;ones(size(sparsity{i}))*i];
%     jj = [jj;sparsity{i}];
%     ss = [ss Y(i,sparsity{i}) - (sum(Y(i,sparsity{i}))-1)/length(sparsity{i})];
% end
% [n,m] = size(Y);
% out = sparse(ii, jj', ss', n,m);


Y=Y.*bw;
out=Y-bw.*(sum(Y,2)-1)./sum(bw,2);
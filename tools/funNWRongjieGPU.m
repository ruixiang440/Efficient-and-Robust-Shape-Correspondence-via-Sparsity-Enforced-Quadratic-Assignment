function [F,G] = funNWRongjieGPU(x)

global S1 S2 M1 M2 mu


% S1 = gpuArray(S1);
% S2 = gpuArray(S2);
% M1 = gpuArray(M1);
% M2 = gpuArray(M2);
% tic
% temp1 = x*S1-S2*x;
% temp2 = x*M1-M2*x;
% F = gather(0.5*norm(temp1,'fro')^2+0.5*mu*norm(temp2,'fro')^2);
% G = temp1*S1 - S2*temp1 + mu*(temp2*M1 - M2*temp2);

temp1 = S1*x-x*S2;
temp2 = M1*x-x*M2;
F = gather(0.5*norm(temp1,'fro')^2+0.5*mu*norm(temp2,'fro')^2);
G = S1*temp1 - temp1*S2 + mu*(M1*temp2 - temp2*M2);
% 
% temp1 = S1*x-x*S2;
% F = gather(0.5*norm(temp1,'fro')^2);
% G = S1*temp1 - temp1*S2;


% toc
% S1 = gather(S1);
% S2 = gather(S2);
% M1 = gather(M1);
% M2 = gather(M2);
% G = gather(G);

end
clc;
clear all;
addpath(genpath(pwd))
load('Cora.mat');% change dataset, "pol-book","Hamilton46.mat","Hamilton46.mat" is also allowed.
% here uG and ty is used, uG is Adjacency matrix, ty is truth-label
% if you have use your own dateset without uG and ty, they should be
% reassign or formatted .
%%
[LCC0, lcc_inds0, ci0, sizes0] = LargestConnectedComponent(sparse(uG)); 
% obtain Largest connected component as graph to be processed
A=sparse(LCC0);
ty_lcc=ty(lcc_inds0);
if find(ty_lcc==0)
    ty_lcc=ty_lcc+1;
end
lambda=0.8; % change between [0.0,1.0], larger value indicate more concentrate on Higher-order modularity
[COM,ending]=HoMod(A,1,lambda);
[a,ind]=max(COM.MOD);% choose community partition with max mocularity in hierarchical results.
com=COM.COM{ind};
F1 = F1Over(ty_lcc, com); 
nmi_score= NMI(ty_lcc,com);
pty = Purity(ty_lcc,com);
mod=COM.MOD(ind);
result = [nmi_score,F1,pty,mod];


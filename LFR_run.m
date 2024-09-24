clear all;
clc;
addpath(genpath(pwd))
edgelist = load('.\LFRdata\LFR_edgelist_30000_0.3.txt');
n = max(max(edgelist));
% constrcut sparse matrix
uG = sparse(edgelist(:,1),edgelist(:,2),edgelist(:,3),n,n);
uG = uG' + uG;
uty = load('.\LFRdata\LFR_label_30000_0.3.txt');
uty = (sortrows(uty,1));
ty = uty(:,2);
[LCC0, lcc_inds0, ci0, sizes0] = LargestConnectedComponent(sparse(uG)); 
A = sparse(LCC0);
num_nodes = length(A);
num_edges = sum(sum(A))/2;
ty_lcc=ty(lcc_inds0);%email-Eu-core has label 0
if find(ty_lcc==0)
    ty_lcc=ty_lcc+1;
end
num_com = length(unique(ty_lcc));
% motif instance
W = MotifAdjacency(sparse(A),'m4');
num_motif = sum(sum(W))/6;
cc_all = clustering_coefficients(A);
acc = sum(cc_all)/num_nodes;
description = [num_nodes,num_edges,num_com,num_motif,acc];
%%
lambda=0.2;
[COM,ending]=HoMod(A,1,lambda);
[a,ind]=max(COM.MOD);
com=COM.COM{ind};
F1 = F1Over(ty_lcc, com); 
nmi_score = NMI(ty_lcc, com);
pty = Purity(ty_lcc,com);
Q = COM.MOD(ind(1));
res1 = [nmi_score,F1,pty,Q];




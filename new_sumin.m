% this function is used inside Louvain_HO
function num=new_sumin(i,W,nodes)
% calculate # of triangle including node i and nodes that in nodes set.
% W is motifAdjacency matrix
num=0;
i_honb=find(W(i,:));
a=intersect(i_honb,nodes);
num=sum(W(i,a))/2;%triangle is counted repeatedly
end




     
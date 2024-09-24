function Q=modularity(A,C)
C_t=unique(C);
size_C=length(C_t);
m=sum(sum(A));
Q=0;
for i=1:size_C
    e_ii=0.0;
    a_i=0.0;
    nodes=find(C==C_t(i));
    n1=length(nodes);
    e_ii=sum(sum(A(nodes,nodes)));  
    a_i =sum(sum(A(nodes,:)));
    if a_i>0
        Q=Q+(e_ii/m-(a_i/m)^2);
    end
end
end


    
%this algorithm is completed by recursive methohd 
%INPUTs 
%A:adjacency matrix
%s if s=1 recursive mode;s=0 non-recursive mode
%OUPUT :
%COMMTY :struct 
%in each level i
%COMMTY.COM{i} :community assingment of each node in louvain iteration i
%COMMTY.SIZE{i}:community size
%COMMTY.MOD(i) community modularity
%COMMTY.ITER(i) iteration of process 1 in louvain.
function [COMMTY ending]=HoMod(A,s,lambda)

N=length(A);
ending=0;
debug=0;
ddebug=0;
%symmetrize
%process 1
A=A+A';%symmetrical matrix
M2=A;
M2((N+1).*[0:N-1]+1)=0;%make diagonal zero
m=sum(sum(A));
if (m==0||N==1)
    fprintf('No more possible decomposition\n');
    COMMTY=0;
    ending=1;
    return;
end
K=sum(A);%sum of weight incident to node i
SUM_IN=diag(A);
SUM_TOT=sum(A);
COM=1:N;
for i=1:N
    NB{i}=find(M2(i,:));
end

gain=1;% indicate the execution of while loop
iter=1;
    
%higher-order community structure
W=MotifAdjacency(A,'m4');

db=sum(W)/2;
HSUM_IN=zeros(1,length(W));
HSUM_TOT=db;
mbar=sum(db);

while (gain==1)
    gain=0;
    for i=1:N
        Ci=COM(i);
        nb=NB{i};
        Gain=zeros(1,N);
        most_increase=-1;
        C_new=Ci;
        COM(i) = -1;
        SUM_TOT(Ci)=SUM_TOT(Ci)-K(i);
        CNi=find(COM==Ci);
        SUM_IN(Ci)=SUM_IN(Ci)-2*sum(A(i,CNi))-A(i,i);
        if(iter>1)        
            %higher-order community structure
            HSUM_IN(Ci)=HSUM_IN(Ci)-new_sumin(i,W,CNi);%HSUM_IN>=0
            HSUM_TOT(Ci)=HSUM_TOT(Ci)-db(i);
        end
        for j=1:length(nb)
            Cj=COM(nb(j));
            if (Gain(Cj)==0)
                CNj=find(COM==Cj);
                ho_gain=0;
                lo_gain=(2*sum(A(i,CNj))/m-2*SUM_TOT(Cj)*K(i)/(m*m));
                if(iter>1)
                    d_i_in=new_sumin(i,W,CNj);
                    ho_gain=(d_i_in/(2*mbar)-3*((HSUM_TOT(Cj))^2*db(i)+ HSUM_TOT(Cj)*(db(i))^2)/(mbar)^3);
                    Gain(Cj)=(1-lambda)*lo_gain+lambda*ho_gain;
                else                   
                    Gain(Cj)=lo_gain;
                end
                %lambda_cluster=d_i_in;
                %fprintf(fid,"node %d-> community %d\nhigher-order gain is %g\nlower-order gain is %g\n",i,Cj,ho_gain,lo_gain);
                %Gain(Cj)=(1-lambda)*lo_gain+lambda*ho_gain;
                %fprintf("higher-order gain is %g\nlower-order gain is %g\n",full(ho_gain),full(lo_gain));
                if Gain(Cj)>most_increase
                    most_increase=Gain(Cj);
                    Cnew_t=Cj;
                end
            end
        end
        if most_increase > 0
            C_new = Cnew_t;
            if (debug)
                fprintf('Move %d => %d\n',i,C_new);
            end
        end
        Cn=find(COM==C_new);
        SUM_IN(C_new)=SUM_IN(C_new)+2*sum(A(i,Cn));
        SUM_TOT(C_new)=SUM_TOT(C_new)+K(i);
        if(iter>1)
            HSUM_IN(C_new)=HSUM_IN(C_new)+new_sumin(i,W,Cn);
            HSUM_TOT(C_new)=HSUM_TOT(C_new)+db(i);
        end
        COM(i)=C_new;
        if(C_new~=Ci)
            gain=1;
        end
    end
    iter=iter+1;
end
iter=iter-1;
[COM,SIZE_T]=reindex(COM);
MOD=modularity(A,COM);
COMMTY.COM{1}=COM;
COMMTY.SIZE{1}=SIZE_T;
COMMTY.MOD(1)=MOD;
COMMTY.ITER(1)=iter;
fprintf('current nodes num is %d\n',N);
%{
if (s~=0)
    fprintf(fid,'----------iter 1----------\n');
end
%}
%process 2
if(s==1)
    k=2;
    curCOM=COM;%size Node_cur, change during loop
    fullCOM=COM;%size N,don't change size in 'while' loop
    newM=A;
    oldM=newM;
    while 1
        oldM=newM;%oldM save information in last iteration
        Nnode=length(oldM);
        Ncom=length(unique(curCOM));
        ind_com=zeros(Ncom,Nnode);
        ind_com_full=zeros(Ncom,N);%memorize nodes set in last partition
        for p=1:Ncom
            set_com=find(curCOM==p);
            ind_com(p,1:length(set_com))=set_com;
        end
        for p=1:Ncom
            set_com=find(fullCOM==p);
            ind_com_full(p,1:length(set_com))=set_com;
        end
        newM=zeros(Ncom,Ncom);
        for m=1:Ncom
            for n=m:Ncom
                indm=ind_com(m,:);
                indn=ind_com(n,:);
                newM(m,n)=sum(sum(oldM(indm(indm>0),indn(indn>0))));
                newM(n,m)=sum(sum(oldM(indm(indm>0),indn(indn>0))));
            end
        end
        [COMt,e]=HoMod(newM,0,lambda);%recursive 
        if(e ~= 1)%not end
            fullCOM=zeros(1,N);
            curCOM=COMt.COM{1};%size of nodes in newM
            for p=1:Ncom
                indp=ind_com_full(p,:);
                fullCOM(indp(indp>0))=curCOM(p);
            end
            [fullCOM,size_t2]=reindex(fullCOM);
            COMMTY.COM{k}=fullCOM;
            COMMTY.SIZE{k}=size_t2;
            COMMTY.MOD(k)=modularity(A,fullCOM);
            COMMTY.ITER(k)=COMt.ITER(1);
            ind=(fullCOM==COMMTY.COM{k-1});
            %fprintf(fid,"-----------------iter %d-------------------\n",k);
            if(sum(ind)==length(ind))%no improvements
                return;
            end
            if (debug)
                fprintf('Identical segmentation => End\n');
            end
        else
            if (debug)
                fprintf('Empty matrix => End\n');
            end
            return;
        end
        k=k+1;
    end
end       
%fclose(fid);
end
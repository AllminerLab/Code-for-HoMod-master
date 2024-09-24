function [COM,sz]=reindex(COMold)
ucom=unique(COMold);
COM=zeros(1,length(COMold));
S=zeros(1,length(ucom));
for i=1:length(ucom)
    S(i)=length(COMold(COMold==ucom(i)));
end
[sz,index]=sort(S,'descend');
for i=1:length(ucom)
    COM(COMold==ucom(index(i)))=i;
end

   
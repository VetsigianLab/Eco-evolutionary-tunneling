function [ stoptime ] = identification_two( Traj,w )



slopelar1=zeros(16,size(Traj,2)-w+1);
intercept=zeros(16,size(Traj,2)-w+1);
ratiolar=zeros(16,size(Traj,2)-w+1);

x=[ones(w,1),[0:(w-1)]'];
ESCstrains=[3,7,9,10,11];

for xcv=1:5
    j=ESCstrains(xcv);
    
    for i=w:1:size(Traj,2)
        y=Traj(j,i-(w-1):i)';
        if sum(y<10^-4)==0
            a=x\y;
          
            slopelar1(j,i-(w-1))=abs(a(2));
            intercept(j,i-(w-1))=a(1);
            ratiolar(j,i-(w-1))=abs(a(2))/max(a(1),a(1)+w*a(2));
        else
            ratiolar(j,i-(w-1))=NaN;
        end
        
        
    end
end

ratiolar(:,end-w+1:end)=1;

q=find(sum(ratiolar([3,7,9,10,11],:)<0.2)==5,1,'first');
    if ~isempty(q)
           stoptime=q;
    else
        stoptime=Inf;
           
    end
    
end


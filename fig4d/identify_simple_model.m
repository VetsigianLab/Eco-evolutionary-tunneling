function [ R,x1,y1] = identify_simple_model( TrajX,N )
w=10;
trial=size(TrajX,2);
global maske

maske=(0:(w-1))'+((1:1:(size(TrajX,1)-w+1)));



for qaz=1:trial
    
     stoptime=NaN;
    
        
         Traj=squeeze(TrajX(:,qaz,:))';


         q=visual_function_basit_yeni_one_v3(Traj,w);
         ss=find(q,1,'first');
         
         if ~isempty(ss)
            stoptime=ss;
         else
             stoptime=Inf;
         end
           stoptimelar(qaz)= stoptime;
         
     qaz
   
end

stoptimelar=(stoptimelar-1)*10+1;
Y1=stoptimelar;
try
YY=Y1((Y1>0)&(Y1<300000));[f,x1]=ecdf(YY);f=f*length(YY)/(length(Y1));y1=log(1-f);

[a,b]=createFit(x1,y1);
title(['N=',num2str(N)])
R=a.p1;
catch
   R=0; x1=0;y1=0;
end
end


function [ A,B,C ] = evolution_simple( A2,B2,C2,M,N )
%log appx -C1-(1/2)*C1^2-(1/3)*C1^3

% change of variables
S=C2;
P=A2;
Deg=B2;


%growth rates in absance of any interactions
Rp=1.6;
Rd=1.81;
Rs=1.9;

% killing coeffs
kp=2*0.55;
kpp=1*1;

%eqsns
A1=Rp*(P+(kp)*1*P*S^2+(kpp)*1*S*P^2);

B1=Rd*(Deg);

C1=Rs*(S-1*2*P*S^2-1*P^2*S);

%normalization
T=A1+B1+C1;
A1=A1/T;
B1=B1/T;
C1=C1/T;



%appx mnrnd
if A1<B1
   if A1>C1
       
       C=binorandim(N,C1);
       A=binorandim(N-C,A1/(A1+B1));
       B=N-(A+C);
   else
       if B1<C1
           A=binorandim(N,A1);
           B=binorandim(N-A,B1/(C1+B1));
           C=N-(A+B);
       else
           A=binorandim(N,A1);
           C=binorandim(N-A,C1/(C1+B1));
           B=N-(A+C);
       end
   end
       
else
    if A1<C1
       B=binorandim(N,B1);
       A=binorandim(N-B,A1/(A1+C1));
       C=N-(A+B);
       
    else
        if B1>C1
            C=binorandim(N,C1);
            B=binorandim(N-C,B1/(A1+B1));
            A=N-(C+B);
        else 
            B=binorandim(N,B1);
            C=binorandim(N-B,C1/(A1+C1));
            A=N-(C+B);
        end
    end
    
end


%poisson rnd number generator
L=-M;
m=double(-1);
p=double(0);
while p>L
    u=log(rand('double'));
    p=p+u;
    m=m+1;
end


%calculation of mutants
for i=1:m
   x=rand('double');
   y=rand('double');
   
   if x<(A/N)
       A=A-1;
       if y<0.5
           B=B+1;
       else
           C=C+1;
       end
       
   elseif x<((A+B)/N)
       B=B-1;
       if y<0.5
           A=A+1;
       else
           C=C+1;
       end
       
   else
       C=C-1;
       if y<0.5
           A=A+1;
       else
           B=B+1;
       end
   end
    
end

A=round(A)/N;
B=round(B)/N;
C=round(C)/N;


end


function [x]=binorandim(N,p)
    
    if N*p>50
        x=round((randn('double')*sqrt(N*(1-p)*p)+N*p));
        
      
    else
    if p>0.01
        log_q=log(1-p);
    else
        log_q=-p-(1/2)*p^2-(1/3)*p^3;
    end
    
    x=double(0);
    summ=double(0);
    while summ>log_q
        summ=summ+log(rand('double'))/(N-x);
        x=x+1;
    end
    x=x-1;
    if x<0
       x=x+1; 
    end
  
    end


end

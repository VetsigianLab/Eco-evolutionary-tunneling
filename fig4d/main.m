clearvars -except r1

%initialization of variables
mler=gpuArray(double(0.008));
NN=gpuArray(double     (round(10.^[4:.2:8] )));
trial=gpuArray(double( 60000  ));
T=gpuArray(double(50000));

rng('shuffle')
totaltrial=trial*length(NN);

A=gpuArray(double(zeros(1,totaltrial)));
B=gpuArray(double(ones(1,totaltrial)));
C=1-(A+B);
N=[];
for i=1:length(NN)
    N=[N,NN(i)*gpuArray(double(ones(1,trial)))];
end

i=gpuArray(double(1));
ind=1;

m=mler(1)*gpuArray(double(ones(1,totaltrial)));
AA=(double(ones(floor(T/10)+1,totaltrial)));
BB=(double(ones(floor(T/10)+1,totaltrial)));

ts=1;
for i=1:T

%main evolution function
    [ A,B,C ] = arrayfun(@evolution_simple,A,B,C,m,N );

    
    if mod(ind,10)==1
        ts=ts+1;
        AA(ts,:)=gather(A);
        BB(ts,:)=gather(B);
    end

    i
    ind=ind+1;
end

%saving and loading takes long time

%save([num2str(zxc)],'AA','BB','-v7.3')

%Identifying the snapping times
trial=gather(trial);
for zxc=[1:length(NN)]


Trajlar=zeros([5001,trial,3]);
Trajlar(:,:,1)=AA(:,(zxc-1)*trial+1:(zxc)*trial);
Trajlar(:,:,2)=BB(:,(zxc-1)*trial+1:(zxc)*trial);
Trajlar(:,:,3)=1-(Trajlar(:,:,1)+Trajlar(:,:,2));




%function that identifies the community formation time.
[r,x1,y1]=identify_simple_model(Trajlar,NN(zxc));

x{zxc}=x1;
y{zxc}=y1;
r1(zxc)=r;



end

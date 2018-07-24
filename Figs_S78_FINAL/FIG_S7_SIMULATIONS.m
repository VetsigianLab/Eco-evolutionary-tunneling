
%These are the simulations underlying Fig_S7
%They will take a very long time to run, especially at large N.
%Notice the "parfor r0 = 1 : 28" loop. This is because the simulations were
%run on a machine with 28 cores.
%change parfor to for for non-parallel


exp_name = 'sym_a01_b03_v2';
alpha = 0.1;
beta = 0.3;
f2=1;
f3=1;
pwr = 1;

%the model is symmetric, so if a 3-strain community forms, the fixed point
%will be in the middle
xstar = [1;1;1]/3; 


N_ = [1000 3000 10000 30000 100000 300000 1000000 3000000 10000000];



%%

%Start at strain 1 and introduce one individual of strain 2. Simulate
%until strain 1 or 2 is lost. Record the probability to find the system in
%state x2, where x2 in (0,1) is the relative abundance of strain 2. 
%Px_{ni} contains N_(ni)-1 time points corresponding to number of times the
%system was at x2 = (0,1,2,...N-1)/N.
%
% TT_(ni) = sum(Px_{ni}) - total time spend in transit from 1-->2 for all
% replicate simulations

Px_ = cell(length(N_),1);



for ni = 1 : length(N_)
   
    N = N_(ni);
 
    Px = zeros(1, N-1);
    
    n0 = 1;
    
    parfor r0 = 1 : 28
        F = zeros(1,2);
       
        tmp = nan(1,2);
        WW = nan(1,2);
        Px0 = zeros(1, N-1);
        for r = 1 : 40000
            x = [N-n0; n0]/N;
            
            while 1 < 2
                
                %the fitness values of the two strains
                F(1) = 1 + alpha * x(1);
                F(2) = f2 + alpha * x(2) + beta * x(1);
                
                %replicator equation
                W = F.*x'./(F*x);
                
                
                [sW, si] = sort(W);
                
                %mean and standard deviation for binomial sampling
                sWN = N*sW;
                std_s = sqrt(sWN.*(1-sW));
                
                
                q = 1;
                %binomial sampling between generations is approximated by 
                %poisson (for low mean) or normal distributions to
                %speed up the simulations. Notice that these approximations
                %are super accurate for large N values used
                if sWN(q) < 100
                    WW(q) = poissrnd(sWN(q));
                else
                    WW(q) = round(sWN(q) + std_s(q)*randn);
                end
                
                
                WW(2) = N - WW(1);
                tmp(si) = WW;
                
                
                
                x = tmp'/N;
                if any(tmp==0)
                    if mod(r,1000)==0
                        disp([ni r0 r])
                    end
                    break
                end
                
                Px0(tmp(2)) = Px0(tmp(2)) + 1;
            end
            
        end
        Px = Px + Px0;
    end
    Px_{ni} = Px;
    TT_ = cellfun(@(x) sum(x), Px_);
    save([exp_name '_Px'], 'Px_', 'N_', 'TT_', 'alpha', 'beta', 'f2', 'f3', 'pwr')
end

%%
 
%From every possible x2 (for small N), and a uniform subset in (0,1) for large N,
%a single individual from strain 3 is introduced and the dynamics is
%followed to determine if a 3-strain community is reached. Many replicates
%are run (R_(ni)) from each x2 to determine the probability of community formation
%given a strain 3 mutatant appears when the system is at x2.
%The probability is SNAP_{ni}/R_(ni). At large N some positions are zero
%because no simulations were run at that x2.


SNAP_ = {}; 

R_ = round(10e6./N_);
for ni = 1 : length(N_)
   
    N = N_(ni);
    tmp = nan(1,3);
    WW = nan(1,3);
    
    R =  R_(ni);
    
    snap = zeros(1, N-1);
    
    parfor zi = 1 : 28
        snap0 = zeros(1, N-1);
        F = zeros(1,3);
        tmp = nan(1,3);
        WW = nan(1,3);
        for X2 = zi : 28: N-1
            
            x0 = [N-X2-1; X2; 1]/N;
            
            if mod(X2,1000) == 0
                disp([zi ni X2])
            end
            
            for r = 1: R
                
                x = x0;
                is_snap = 1;
                while 1 < 2
                    
                    F(1) = 1 + alpha * x(1) + beta * x(3);
                    F(2) = f2 + alpha * x(2) + beta * x(1);
                    F(3) = f3 + alpha * x(3) + beta * x(2);
                    
                    W = F.*x'./(F*x);
                    
                    [sW, si] = sort(W);
                    
                    sWN = N*sW;
                    std_s = sqrt(sWN.*(1-sW));
                    
                    
                    for q = 1 : 2
                        if sWN(q) < 100
                            WW(q) = poissrnd(sWN(q));
                        else
                            WW(q) = round(sWN(q) + std_s(q)*randn);
                        end
                    end
                    
                    WW(3) = N - WW(1)-WW(2);
                    tmp(si) = WW;
                    
                    x = tmp'/N;
                    
                    if any(tmp==0)
                        is_snap = 0;
                        break
                        %there is a danger of an infinite loop if there is no convergence to xstar
                    elseif max(abs(x - xstar))<0.02
                        break
                    end
                    
                end
                
                
                snap0(X2) = snap0(X2) + is_snap;
                
            end
            
        end
        snap = snap + snap0;
    end
    
  SNAP_{ni} = snap;
  save([exp_name '_SNAP'], 'N_', 'R_', 'SNAP_','alpha', 'beta', 'f2', 'f3', 'pwr')
end
    


%this is an example of how Fig S8 data is generated
%see the corresponding file for Fig S7 for comments on the code

alpha = 0;
beta = 0.3;
f2=1;
f3=1;
xstar = [1;1;1]/3;

pwr_ = [1 1.5 2 2.5];
List_ = {'sym_a0_b03_v3', 'sym_a0_b03_pwr_15_v3', 'sym_a0_b03_pwr_20_v3', 'sym_a0_b03_pwr_25_v3'};

for li = 1 : length(pwr_)
    
    pwr = pwr_(li);
    
    N_ = [1000 3000 10000 30000 100000 300000 1000000 3000000 10000000];
    
    
    exp_name =  List_{li};
   
    
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
                    
                    
                    F(1) = 1;
                    F(2) = f2  + beta * x(1).^pwr;
                    
                    
                    W = F.*x'./(F*x);
                    
                    %tmp = mnrnd(N,W);
                    
                    [sW, si] = sort(W);
                    
                    sWN = N*sW;
                    std_s = sqrt(sWN.*(1-sW));
                    
                    
                    q = 1;

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
    
    exp_name = List_{li};
    SNAP_ = {};
    R_ = round(40e6./N_);
    
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
                        
                        F(1) = 1  + beta * x(3).^pwr;
                        F(2) = f2  + beta * x(1).^pwr;
                        F(3) = f3  + beta * x(2).^pwr;
                        
                        W = F.*x'./(F*x);
                        
                        %tmp = mnrnd(N,W);
                        
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
    
    
end


function simulation = evolve_RS_MutOne_Space(fname, simulation, fname_OUT)

%Multiple localities/patches connected by migration are added to
%eco-evolutionary simulations.
%New parameters: simultion.PatchNum, simulation.MigrationMatrix
%if the MigrationMatrix is a scalar, assume that everyone is connected to
%everyone and construct the matrix.


global R_ 
global PROD_Gaussian RES_Gaussian DEG_Gaussian
global PROD_circ_pad
global RES_circ_pad

if nargin < 2 || isempty(simulation)
    load(fname, 'simulation')
end
if nargin < 3 || isempty(fname_OUT)
    fname_OUT = fname;
end

if isfield(simulation, 'rng_state')
   rng(simulation.rng_state)
elseif isfield(simulation, 'random_seed')
    rng(simulation.random_seed);
else
    simulation.random_seed = rng('shuffle');  
    rng(simulation.random_seed);
end


%mu - mutation rate per population per generation
%Pb is the probability to mutate production rather than attenuation
%Cp - cost of production
%Ca - cost of attenuating resistance


%default values
if ~isfield(simulation, 'max_change_factor')
    simulation.max_change_factor = 2;
end
if ~isfield(simulation, 'mutant_freq')
    simulation.mutant_freq = 1E-4;
end
if ~isfield(simulation, 'min_abundance')
   simulation.min_abundance = 1E-6;
end
if ~isfield(simulation, 'saveT')
   simulation.saveT = 10;
end


if ~isfield(simulation, 'Prob_loss_of_function_Mut')
    simulation.Prob_loss_of_function_Mut = 0.01; %loss of function towards S or R
end


if ~isfield(simulation, 'Prob_global_mut') 
  simulation.Prob_global_mut = 0.5; %with this prob a uniform phenotype is chosen. with 1-p probability phenotype changes by a log normal factor.

end


if ~isfield(simulation, 'TuningMutSize') % Three regions: production, attenuation, sensitivity
  simulation.TuningMutSize = 2; %the std of multiplactive factor is 2: most mutations change phenotype between 50% and 200%
end

if ~isfield(simulation, 'maxP')
  simulation.maxP = 1./simulation.Cp; 
end

if ~isfield(simulation, 'maxA')
  simulation.maxA = 1./simulation.Ca; 
end

if ~isfield(simulation, 'maxR')
  simulation.maxR = 1./simulation.Cr; 
end


if ~isfield(simulation, 'timesteps_completed')
  simulation.timesteps_completed = 0; 
end

if isscalar(simulation.MigrationMatrix)
   m =  simulation.MigrationMatrix;
   simulation.MigrationMatrix = m*ones(simulation.PatchNum)/(simulation.PatchNum-1);
   simulation.MigrationMatrix(1:(simulation.PatchNum+1):end) = 1-m;
end

if ~isfield(simulation, 'X')
  %simulation.X = cell(simulation.T+1,1); 
  simulation.X{1} = ones(1, simulation.PatchNum);
end

if ~isfield(simulation, 'Num_perm')
    simulation.Num_perm = 1;
end

if ~isfield(simulation, 'dose_response_factor')
    simulation.dose_response_factor = 20; %very steep
end

%other simulation fileds:
%X{ti} holds the relative abundance at every time point
%Phenotypes{ti} holds the phenotypes at every time point

%%


%the size of the PROD  and DEG filters needs to be the same in the current implementation
PROD_Gaussian = fspecial('gaussian', simulation.filter_size, simulation.PROD_sigma);
RES_Gaussian = fspecial('gaussian', simulation.filter_size, simulation.RES_sigma); %refers to the resources
DEG_Gaussian = fspecial('gaussian', simulation.filter_size, simulation.DEG_sigma);

[R_, simultion.XX_] = precalc_05(simulation, PROD_Gaussian);


dimNums = uint32(1:simulation.Nz);
padSize = (size(PROD_Gaussian,1)-1)/2;
PROD_circ_pad = dimNums(mod(-padSize:simulation.Nz+padSize-1, simulation.Nz) + 1);

padSize = (size(RES_Gaussian,1)-1)/2;
RES_circ_pad = dimNums(mod(-padSize:simulation.Nz+padSize-1, simulation.Nz) + 1);

%gaussian filters are separable and can be implemented as two 1D passes
%this removes significant overhead in using filter2 instead of conv2 in GPU calculations
PROD_Gaussian = decompose_filter(PROD_Gaussian);
RES_Gaussian = decompose_filter(RES_Gaussian);
DEG_Gaussian = decompose_filter(DEG_Gaussian);



disp('init done')
%%


%FIXED NUMBER OF ANTIBIOTICS
Na = size(simulation.Phenotypes{1},1);

%Continue an existing simulation from the last saved time-point
ts = 1 + simulation.timesteps_completed/simulation.saveT; 
X = simulation.X{ts};
p = double(simulation.Phenotypes{ts});

g = 1 - simulation.Cp0*sum(p>0,1) - simulation.Ca0*sum(p<-1,1) - simulation.Cp*sum((p>1).*(p-1),1) - simulation.Ca*sum((p<-1).*(abs(p)-1),1);
g = g';

Prob_GLOBAL = simulation.Prob_loss_of_function_Mut + (1-simulation.Prob_loss_of_function_Mut)*simulation.Prob_global_mut;


for t = (simulation.timesteps_completed+1) : simulation.T
    
    
    %introduce mutations between ecological cycles
    Nmut = poissrnd(simulation.mu);
 

    for mi = 1 : Nmut
        X = bsxfun(@times, X, 1./sum(X));
        
        %generate random locality
        mut_loc = randi(simulation.PatchNum);
        
        parent = find(mnrnd(1, X(:, mut_loc)')); %select a random individual to mutate
        

        
        %select random antibiotic
        z = randi(Na); 
        
        p1 = p(:,parent);        
        g1 = -1;


        
        if abs(p1(z)) == 1
            while g1 < 0
                rand1 = rand;
                if rand1 < simulation.Prob_loss_of_function_Mut
                     p1(z) = -p(z,parent); %swith between S and R
                elseif rand1 < Prob_GLOBAL
                    %global mutation
                    tmp = 2*rand - 1;
                    p1(z)=(tmp<0).*(-1+simulation.maxA*tmp) + (tmp>=0).*(1+simulation.maxP*tmp);
                else
                    break; % no mutation
                end
                g1 = 1 - simulation.Cp0*sum(p1>0,1) - simulation.Ca0*sum(p1<-1,1) - simulation.Cp*sum((p1>1).*(p1-1),1) - simulation.Ca*sum((p1<-1).*(abs(p1)-1),1);


            end
            
        else
            while g1 < 0
                rand1 = rand;
                if rand1 < simulation.Prob_loss_of_function_Mut
                    p1(z) = 2*(rand > 0.5) - 1;
                elseif rand1 < Prob_GLOBAL
                    tmp = 2*rand - 1;
                    p1(z)=(tmp<0).*(-1+simulation.maxA*tmp) + (tmp>=0).*(1+simulation.maxP*tmp);
                else
                    p0=p(z,parent);
                    p1(z) = sign(p0).*((abs(p0)-1).*exp(simulation.TuningMutSize*randn)+1); 
                end
                
                g1 = 1 - simulation.Cp0*sum(p1>0,1) - simulation.Ca0*sum(p1<-1,1) - simulation.Cp*sum((p1>1).*(p1-1),1) - simulation.Ca*sum((p1<-1).*(abs(p1)-1),1);

            end
        end
                
        if g1 > 0  
        p(:,end+1) = p1;
        g(end+1,1) = g1;
        X(end+1,mut_loc) = simulation.mutant_freq;        
        end
    end
    
    X = X * simulation.MigrationMatrix;
    X(X==0) = 1e-10;
    X = bsxfun(@times, X, 1./sum(X));
    

    a = (p<0).*(-p-1);
    b = (p>0).*(p-1);
    
    b(p==1) = 0.0001;
    
    %simulate one time-step (ecological cycle) for each patch
    for loc = 1 : simulation.PatchNum


        f = get_fitness(X(:,loc), g, b', a', simulation.alpha, simulation.Eff_op, simulation.Deg_op, simulation.Gamma, simulation.L.^2, simulation.Nz, ...
            simulation.dose_response_factor, simulation.Tp, simulation.Tr, simulation.Num_perm);
        
        X(:,loc) = f/sum(f);
    end

    

    tmp = all(X<simulation.min_abundance,2);
    X(tmp, :) = [];
    g(tmp) = [];
    p(:,tmp) = [];
    
    

   %record the dynamics every saveT time steps
   if mod(t,simulation.saveT)==0
        disp(t)
        
        %consolidate by combining independently arising identical strains
        [p,~,kq]=unique(p','rows');
        p = p';
        X1 = X;
        X = zeros(max(kq), simulation.PatchNum);
        for loc = 1 : simulation.PatchNum
             X(:,loc) = accumarray(kq,X1(:,loc), [size(p,2) 1]);
        end
        
        g = 1 - simulation.Cp0*sum(p>0,1) - simulation.Ca0*sum(p<-1,1) - simulation.Cp*sum((p>1).*(p-1),1) - simulation.Ca*sum((p<-1).*(abs(p)-1),1);
        g = g';
        
        
        ts = 1 + t/simulation.saveT; % the first indx stores the initial conditions
        simulation.timesteps_completed = t; 

        simulation.Phenotypes{ts} = single(p); %store in single precision to save space
    
        simulation.X{ts} = X;
        
        simulation.rng_state = rng;
        
        %while the data is recored every saveT time-steps, it is saved to
        %the disk only every saveMultiplier*saveT time-steps
        if mod(t,simulation.saveMultiplier*simulation.saveT)==0
                save(fname_OUT, 'simulation');
        end
     
        %if running purely ecological dynamics (mu=0) and there is only one
        %strain left, terminate the simulation
        if simulation.mu == 0 && size(X,1) == 1
            return;
        end

        X(X==0) = 1e-10;
        X = bsxfun(@times, X, 1./sum(X));

        %print the current diversity
        disp('Shannon Diversity')
        disp(exp(-sum(X.*log(X))))
    end
end


% a helper function that allows us to speed up the application for gaussian
% filter for diffusion
function [hcol, hrow] = decompose_filter(F)
stencil = rot90(F,2);
[u,s,v] = svd(stencil,'econ');
hcol = u(:,1) * sqrt(s(1));
%hrow = conj(v(:,1)) * sqrt(s(1));



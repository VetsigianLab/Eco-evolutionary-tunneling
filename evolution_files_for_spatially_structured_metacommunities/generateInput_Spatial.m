function [simulation]=generateInput_Spatial

%This file asks the user to key in the various variables to generate the Input file and variable 'sim'. Some of the paramters are 
%set to default number, and can be changed from the generalInput.m by overwriting.
%
%When run the file asks for:
%
%1) Enter the number of antibiotics that can be present in the sytem: (key in a number)
%2) Enter the number different phenotype you would like to key in: (key in a number)
%3) Enter the phenoypes you want to put in for species #1 : (key in an array of length that is equal to the answer of the first question, 
%where each entry is a number, [-Inf,-1) for degrader phenotype, [-1] for sensetive, [1] resistant, (1,Inf] for 
%antibiotic producer) (for example for two antibiotics, '[100,100]' is a strain that produces both antibiotics, '[1,-100]' is a strain 
%that is resistant to first antibiotic and degrader of the second antibiotic). 
%4) Enter the abundance of the phenotype #1:
%
%(repeat 3 and 4 for the number which was the answer of the second question)

%abundance of introduced mutant
simulation.mutant_freq=4e-06;
%min abundance limit
simulation.min_abundance=4e-7;

%maximum level of production
simulation.maxP=1000;
%maximum level of degradation
simulation.maxA=1000;


%number of spores per patch dimension, 2D patches have L^2 spores to start
%the feast and famine cycles
simulation.L=60;

%physical dimensions of patch dimension, Nz=100, L=60 => Nz^2/L^2=2.8, on
%average the size of a colony is 2.8 'pixels'
simulation.Nz=100;

%Number of Patches
simulation.NN=10; %since this simulations take more time, the NN is lowered to 10

%
simulation.alpha=0.5;

%paramter that determines the speed of diffusion of antibiotic in one step. Higher the
%number flatter the filter will be
simulation.PROD_sigma=1;

%paramter that determines the speed of diffusion of degrader molecules in one step (only executed ones). Higher the
%number flatter the filter will be
simulation.DEG_sigma=1;

%paramter that determines the speed of nutriens in one step. Higher the
%number flatter the filter will be
simulation.RES_sigma=1;


%size of the diffusion filters
simulation.filter_size=5;

%number of antibiotic diffusion/degradation cycles
simulation.Tp=40;

%number of resource diffusion/comsumption cycles for the antibiotics
simulation.Tr=10;

%The steepness of the dose-response curve
simulation.dose_response_factor=20;

%Number of permutaitions on using precomputed spaces
simulation.Num_perm=1;

%factor controlling the rate of consumption versus the rate of
% diffusion for resources.
simulation.Gamma=0.3;

% Probability that the mutation will be a loss of function mutation
simulation.Prob_loss_of_function_Mut=0.33;

%simulation.Prob_global_mut=Probability that the mutation will be a global mutation+simulation.Prob_loss_of_function_Mut
simulation.Prob_global_mut=0.5;

%is a constant specifying the weight of degradation opperational costs
simulation.Deg_op=0;

%is a constant specifying the weight of efflux opperational costs
simulation.Eff_op=0;

% Initial fitness cost of being a degrader phenotype
simulation.Ca0=0.06;

% Fitness cost assciated with degree of degradation
simulation.Ca=4e-5;

% Initial fitness cost of being a producer/resistant phenotype
simulation.Cp0=0.3;

% Fitness cost assciated with degree of production
simulation.Cp=8e-4;


% time intervals of saving into simulation.Phenotype and simulation.X
simulation.saveT=1;

%the std of multiplactive factor is 2: most mutations change phenotype between 50% and 200%
simulation.TuningMutSize=2;

% number of timesteps completed.
simulation.timesteps_completed=0;

% time intervals to save the data into hardisk, simulation.saveT*simulation.saveMultiplier
simulation.saveMultiplier=100;




success=0;
while success==0
    prompt='1) Enter the number feast and famine cycles you want to simulate: (key in a number) > (example=10000)';
    temp=input(prompt);
    
    if (rem(temp,1)==0)&&temp>0
       success=1;
       simulation.T=temp;
    else
         fprintf('error\n');
    end
end
    
success=0;
while success==0
    prompt='2) Enter average number of mutants (simulated stochastically) to be introduced at start of each cycle: > (default=[1])';
    temp=input(prompt);

    if (temp>=0)
        simulation.mu=temp;
        success=1;
    else
        fprintf('error\n');
    end
end



simulation.m=1e-5; %migration rate


success=0;
while success==0
    prompt='3) Enter the number sub-communities to simulate: (key in a number) > (example=10)';
    temp=input(prompt);
    
    if (rem(temp,1)==0)&&temp>0
       success=1;
       simulation.PatchNum=temp;
    else
         fprintf('error\n');
    end
end



simulation.timesteps_completed = 0;
simulation.saveMultiplier = 100;
simulation.saveT=1;

simulation.min_abundance=simulation.min_abundance*4; %since the NN is 4 times less than our regular simulations (40 vs 10)
simulation.mutant_freq=simulation.mutant_freq*4;






success=0;
while success==0
    prompt='4) Enter the number of antibiotics that can be present in the sytem: (key in a number, 1 or 2) > (example=1)';
    temp=input(prompt);
    if (temp==[1])||(temp==[2])
        simulation.NumAnti=temp;
        success=1;
    else
        fprintf('error\n')
    end
end







if simulation.NumAnti==1
simulation.Phenotypes = {[-1]};%every sub-community starts from the wild type

 % Initial fitness cost of being a degrader phenotype
simulation.Ca0=0.06;

% Fitness cost assciated with degree of degradation
simulation.Ca=4e-5;

% Initial fitness cost of being a producer/resistant phenotype
simulation.Cp0=0.3;

% Fitness cost assciated with degree of production
simulation.Cp=8e-4;
   
elseif simulation.NumAnti==2
    simulation.Phenotypes = {[-1;-1]};%every sub-community starts from the wild type
    simulation.Ca0=0.06;
    simulation.Ca=8e-5;
    simulation.Cp0=0.35;
    simulation.Cp=2e-4;
else
    fprintf('error!');
    return
end
 



simulation.Traj=zeros(4^size(simulation.Phenotypes{1},1),10,simulation.PatchNum); %initialization
simulation.X = {([ones(1,simulation.PatchNum)])};




save('sim','simulation')


end


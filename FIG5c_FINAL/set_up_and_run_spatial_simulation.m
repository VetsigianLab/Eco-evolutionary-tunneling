%load a non-spatial simulation that reached an ESC
load('ESS_3_153_mu4_r10.mat'); 

%extract the ESC from its last time point
[x, ind] = sort(simulation.X{end}, 'descend');
Ph = simulation.Phenotypes{end}(:,ind);
x = x(1:5); x = x/sum(x); %relative abundance of strains of ESC community
Ph = Ph(:,1:5);

%remove previous simulation data
simulation = rmfield(simulation, 'X');
simulation = rmfield(simulation, 'Phenotypes');
simulation = rmfield(simulation, 'rng_state');
simulation = rmfield(simulation, 'random_seed');
simulation.timesteps_completed = 0;

%put ESC phenotypes at the start of the spatial sim 
simulation.Phenotypes{1} = Ph;

simulation.PatchNum = 20;
X = zeros(length(x),simulation.PatchNum);
X(3,:) = 1; %strain 3 of ESC is SS. Put SS in all patched
X(:,10) = x; % Put the ESC community in patch #10 (The patches form a ring)
simulation.X{1} = X;


simulation.T = 40000; %run for maximum of 40000 steps. Can safely kill/restart as needed (see comment at bottom)
simulation.saveT = 1; %record every time step
simulation.saveMultiplier = 10; %save to file every 10 time-steps
simulation.Nz = 315;
simulation.NN = 4; % 142884 population size = simulation.L^2*simulation.NN
simulation.L = 3.15 * 60;
%simulation.mutant_freq = 1 /(simulation.L^2*simulation.NN);



%set up a migration matrix with ring topology
m = 0.001; %migration rate between neighboring communities
m_ = m*ones(simulation.PatchNum,1);
M=diag(1-m_) + diag(m_(1:end-1)/2,1) + + diag(m_(1:end-1)/2,-1);
M(1,end) = m/2;
M(end,1) = m/2;
simulation.MigrationMatrix = M;

%set up muN = simulation.mu
simulation.mu = 0.1 * simulation.PatchNum;

%run a simulation and save it as RING_01.mat
evolve_RS_MutOne_Space('RING_01', simulation)

% %to continue an interrupted simulation
% load RING_01
% evolve_RS_MutOne_Space('RING_01_continued', simulation)

%to extend an existing simulation that has completed:
%load RING_01
%simulation.T =  simulation.T + 10000 %(extend by 10000 steps)
%evolve_RS_MutOne_Space('RING_01', simulation)



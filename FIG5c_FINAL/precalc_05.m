function [R_, XX_] = precalc_05(simulation, PROD_Gaussian)
%fixed periodic boundary condition
%precompute colony regions assuming equal colonization rates
%reduce peak memory usage by adding one more loop (would allow us to work with
%larger N and Nz on the cluster where there is 2GB limitation)

Nz = simulation.Nz;  %use a NzxNz descritization to figure fitness from resource competition
           %increasing Nz (while keeping other parameters fixed) should improve accuracy 
NN = simulation.NN; %number of spaces.
%We support multiple areas because the fitness function chokes for many
%points and large grids. N = 100, Nz = 100 works great

L = simulation.L; %"physical size of the grid"
N = L^2; %population size. desniity of cells is one by default
%N/L is the density of spores in the beginning of a cycle


%z = linspace(0, L, Nz); %coordinates of the grid points
z = L/(2*Nz):(L/Nz):L; %coordinates of the grid points


XX_={}; R_ = {}; 
for NNi = 1 : NN
    if ~isfield(simulation, 'XX_')
        XX_{NNi} = L*rand(N, 2); % Randomly scattered spores
    end
    XX = XX_{NNi};
    

    %Compute the distance between all the grid points and spores assuming
    %periodic boundary conditions
    R2 = zeros(Nz,Nz);
    X = XX(:,1); Y = XX(:,2);
  
    [zz, Xz] = meshgrid(z, X);
    tmp = zz - Xz;
    dX1 = (tmp).^2;
    dX2 = (tmp+L).^2;
    dX3 = (tmp-L).^2;
    dX = min(cat(3, dX1, dX2, dX3), [], 3);
    
    [zz, Yz] = meshgrid(z, Y);
    tmp = zz - Yz;
    dY1 = (tmp).^2;
    dY2 = (tmp+L).^2;
    dY3 = (tmp-L).^2;
    dY = min(cat(3, dY1, dY2, dY3), [], 3);
    
 
    for ix = 1 : Nz
        dX_ = dX(:, ix);
        [~, q] = min(bsxfun(@plus, dY, dX_));
        R2(ix, :) = q;
    end


    R_{NNi} = padarray(R2, (size(PROD_Gaussian)-1)/2, 'circular');
     
 
end %done with pre-computing



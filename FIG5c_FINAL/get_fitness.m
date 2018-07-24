function final_yield = get_fitness(x, g,  PROD, DEG, alpha, Eff_op, Deg_op, Gamma, N, Nz, dose_response_factor, Tp, Tr, Num_perm)

%computed the yield of all strains after one ecological cycle

global R_ %pre-computed colony clash matrices; from precalc_04.m
global PROD_Gaussian RES_Gaussian DEG_Gaussian
global PROD_circ_pad
global RES_circ_pad

% x is a column vector specifying the relative abundance of each species
%g is a column vector specifying the mycellium growth rate of species
%PROD is a (species, antibiotic) matrix speicfying the production level of an antibiotic by a species
%EFFLUX is a (species, antibiotic) matrix speicfying the Factor by which
%the internal concentration is reduced
%DEG is a (species, antibiotic) matrix speicfying the antibiotic degradation level
%alpha is a number controlling the relative weight of availabel resources
%if not inhibited and those released from killed cells.



%N = size(R_{1},1); % number of spores in a periodic boundary condition space
%Nz = sqrt(size(R_{1},2));

% dose_response_factor: The steepness of the dose-response curve
%Tp: number of antibiotic diffusion/degradation cycles
%Tr: number of resource diffusion/comsumption cycles for the antibiotics

%Eff_op : is a constant specifying the weight of efflux opperational costs
%Deg_op: is a constant specifying the weight of degradation opperational costs
%These are both zero in this manuscript

% Gamma: factor controlling the rate of consumption versus the rate of
% diffusion for resources. There is no such factor for the antibiotics
% because they are assumed to be produced in one step.
% The investement costs in efflux and degradation are already incorporated in g

final_yield = zeros(size(PROD,1),1);
DRP = dose_response_factor/Tp;


for ni = 1 : length(R_)
    
    for perm_i = 1 : Num_perm % reuse every field several times to increase samling size
        
        %assign spore spots to species in proportion to the species abundances

        [S,~] = find(mnrnd(1, x/sum(x), N)');

        S = reshape(S, N, 1);
        
        
        %Each point in the grid is assigned a species number
        q_ext = S(R_{ni});
        f_size = (length(PROD_Gaussian)-1)/2;
        q = q_ext(f_size+1:end-f_size, f_size+1:end-f_size);
        

        G = g(q); 

        pS = zeros(Nz); %log probability of survival
        for a = 1 : size(PROD,2)
            PRODa = PROD(:,a);
            DEGa = DEG(:,a);
            
            
            A_ext = DEGa (q_ext)/(Tp+1);
            %A = filter2(DEG_Gaussian, A_ext,'valid');
            %two 1D convolutions (for each dimension) are used instead of
            %one 2D Gausssian filter
            A = conv2(DEG_Gaussian, DEG_Gaussian, A_ext,'valid');
            
            exp_A = exp(-A);
            A = Deg_op.*A;
            
            E = double(PRODa(q)>0); % E = 1 means having adaptive efflux, E = 0 means no adaptive efflux
            
            P_ext = PRODa(q_ext);
            
            %P = filter2(PROD_Gaussian, P_ext,'valid'); % can be speed up by calling conv2 directly
            P = conv2(PROD_Gaussian, PROD_Gaussian, P_ext,'valid');
            P = P.* exp_A;
            
            
            for ti = 1 : Tp % number of degradation attenuation cycles
                
                %P_ext = padarray(P, (size(PROD_Gaussian)-1)/2, 'circular');
                P_ext = P(PROD_circ_pad, PROD_circ_pad);
                
                %P = filter2(PROD_Gaussian, P_ext,'valid');
                P = conv2(PROD_Gaussian, PROD_Gaussian, P_ext,'valid');
                
                P = P.*exp_A;
                
                %do some killing of cells
                
                %compute the internal cellular concentration from the drug
                %efflux
                
                Pint = P.*(1-E) + E; %if adaptive immunity Pint = 1, Pint = P otherwise
                
                W = Eff_op.*E.*(P-1);%power consumed for efflux
                
                G = G - max(W,0) - A.*P;
                

                
 %               Pint(Pint < 1) = 1; % 1 is the internal level at which an antibiotic starts having an effect
 %               pS = pS + (1-Pint);
 
                 %pS is the log probability to survive. We are adding the
                 %contributions during antibiotic diffusion, which is
                 %equivalent to multiplying by a survivial factor each time
                 %step
                 pS = pS + min(1-Pint,0);
               
                
            end
            
        end
        
  
        pS = exp(DRP.*pS);

        G(G<0) = 0;
       
        
        %Rate of consumption of resources after antibiotic interactions
        R = (1-alpha + alpha.*(1-pS)); %alpha = 0 : no resources are released from dead guys.
        % alpha < 1/2 ensured that being partially killed is never a good thing
        % (despite resources released)
         
        
        spore_yield = R.*pS;
        
        
        %Diffusion and consumption of leftover resources after antibiotic
        %interactions
        for ti = 1 : Tr
            
            %consumption
            R = R .* exp(-Gamma *pS);
            
            %diffusion of the resource
            %R_ext = padarray(R, (size(RES_Gaussian)-1)/2, 'circular');
            R_ext = R(RES_circ_pad, RES_circ_pad);
            %R  = filter2(RES_Gaussian, R_ext,'valid');
            R  = conv2(RES_Gaussian, RES_Gaussian, R_ext,'valid');
            
            %consumed resources increase spore yield
            spore_yield = spore_yield + R.*pS;
            
            
        end
        
        %incorporate the cost of antibiotic production, degradation and
        %efflux
        spore_yield = spore_yield .* G;
        
        final_yield = final_yield + accumarray(q(:), spore_yield(:)', [length(x) 1]);
        
    end
    
end


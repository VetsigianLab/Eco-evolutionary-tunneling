% 1) load simulation_output_Fig5c.mat 
% 2) Run this file to produce Fig. 5c

%colors for the 16 possible phenotypic classes (see Methods)
col16 = [0 0 0; 0 102 102; 153 153 255; 204 0 204; 102 102 0; 160 160 160; 204 255 153; 0 102 0; 102 255 255; 0 255 102; 255 255 0; 255 128 0; ...
    127 0 255; 0 0 153; 255 51 51; 102 0 0]/255;

figure('Position', [1 41 1840 963], 'Color', 'w');
set(gcf, 'PaperUnits', 'normalized', 'PaperOrientation', 'landscape', 'PaperPosition', [0.05 0.05 0.95 0.95])

%find the number of completed time-steps (ecological cycles)
T1 = find(~cellfun(@isempty,simulation.Phenotypes),1,'last');


%string containing simulation parameters is put on top of the figure
params = ['  Cp=' num2str(simulation.Cp) '; eps=' num2str(simulation.Eff_op) '; Cd=' num2str(simulation.Ca) '; Cd0=' num2str(simulation.Ca0) ...
    '; Cr=' num2str(simulation.Cp0)];
set(gcf, 'Name', params)


T_ = 1 :2: T1; %display every second cycle


for si = 1 : simulation.PatchNum %PatchNum is the number of metacommunities (H)
    k = 1; TR = []; XX = [];
    for t1 = T_
        
        x= simulation.X{t1}(:,si);
        x = x/sum(x);
        
        tmp = double(simulation.Phenotypes{t1});
        tmp(tmp<-1) = (tmp(tmp<-1)+1)/4-1; %The degradation levels are rescaled by a factor of 4 for visualization.
                                           %This is allowed since degrdation and production levels have different units
        
        
        TR = [TR [tmp; t1*ones(1,size(tmp,2))]];
        XX = [XX; x];
        k = k + 1;
        
    end
    
    %do not plot strains that are below 10^-4 relative abundance
    tmp = XX>0.0001;
    XX1 = XX (tmp);
    P1 = TR(1,tmp);
    P2 = TR(2,tmp);
    P_ = TR(1:2,tmp);
    
    %switch to a log scale for P and D levels
    P_ = sign(P_) .* (1+log10(abs(P_)));
    

    %compute the phenotypic class and corresponding color for every strain
    %SARP = 0123
    type_ = nan(size(P_));
    type_(P_ <= 0) = 0; %S
    type_(P_ < -1.001) = 1; %A
    type_(P_ > 0) = 2; %R
    type_(P_ > (1.00+log10(2))) = 3; %P
    
    Col_ind = 1 + type_(1,:) + 4*type_(2,:);
    Col = col16(Col_ind,:);
    
    %plot linages over time for one patch in a custom subplot
    kalin_subplot(1,simulation.PatchNum, si);
    scatter3(P_(1,:), P_(2,:), simulation.saveT*TR(3,tmp), max(round(20 + 6*log10(XX1)),1), Col, 'filled'); grid on;

    %project 2D phenotypes to 1D (using a line at 60 deg)
    view(60,0)
    
    zlim([simulation.saveT*T_(1) simulation.saveT*T_(end)]);

 
 set(gca, 'XTick', []);
 set(gca, 'YTick', []);
 if si > 1
      set(gca, 'ZTickLabel', []);
 end
    
end

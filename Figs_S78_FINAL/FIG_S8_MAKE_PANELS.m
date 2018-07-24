

%pwr is the exponent alpha from Methods, "Alternative models without spatial structure"
pwr_ = [0.5 1 1.5 2 2.5];
List_ = {'sym_a0_b03_pwr_05_v2_Px', 'sym_a0_b03_v4_Px', 'sym_a0_b03_pwr_15_v2_Px', 'sym_a0_b03_pwr_20_v2_Px' 'sym_a0_b03_pwr_25_v2_Px'};
List2_ = {'sym_a0_b03_pwr_05_v2_SNAP', 'sym_a0_b03_v3_SNAP', 'sym_a0_b03_pwr_15_v3_SNAP', 'sym_a0_b03_pwr_20_v3_SNAP', 'sym_a0_b03_pwr_25_v3_SNAP'};




%%

SPtot_ = {};
SPTXtot_ = {};
SPXtot_ = {};
T_ = {};
gamma_ = [];
gamma1_ = [];
NN_ = {};
L = {};
for li = 1 : length(List_)
    
    
    
    L{li} = ['\alpha = ' num2str(pwr_(li))];

    load(List_{li});
    load(List2_{li});
    
    NN_{li} = N_;
    
    SPx = {}; SPxT = {};
    Sx={};
    SPtot = [];
    SPTtot = [];
    SPTXtot = [];
    SPXtot = [];
    for ni = 1 : length(N_)
        
        
        Sx{ni} = SNAP_{ni}/R_(ni); %(1 x N-1)
        %this is the probability of escape and snap of an existing mutant at
        %abundance 1.
        

        SPx{ni} =  Sx{ni}.*Px_{ni}/sum(Px_{ni});
        SPxT{ni} = Sx{ni}.*Px_{ni}; %this incorportes the fact that transitions take longer at large N (the log(N) factor)
        
      
        SPtot(ni) = sum(SPx{ni});
        SPTtot(ni) = sum(SPxT{ni});
        SPTXtot(ni) = sum((1-(1:(N_(ni)-1))/N_(ni)).*SPxT{ni});
        SPXtot(ni) = sum((1-(1:(N_(ni)-1))/N_(ni)).*SPx{ni});
    end
    
    SPtot_{li} = SPtot; 
    T_{li} = TT_/28/40000; %normalize by number of replicas
    
    %extract the the power exponent gamma from the last two data points,
    %which are closest to the power law asymptote for high N
    indx = (length(N_)-1):length(N_);
    p = polyfit(log(N_(indx)), log(T_{li}(indx))', 1);
    gamma_(li) = p(1); %this is for panel S8c
    p = polyfit(log(N_(indx)), log(SPTXtot(indx)), 1);
    gamma1_(li) = p(1); %this is for panel S8f
    SPTXtot_{li} = SPTXtot;
    SPXtot_{li} = SPXtot;
end
%% Panel S8e
figure('Color', 'w', 'Position', [1000  938    560   400]);
semilogx(NN_{1}(1:7), SPTXtot_{1}(1:7), 'o-'); hold on
semilogx(NN_{2}(1:7), SPTXtot_{2}(1:7), 'o-'); hold on

legend(L(1:2), 'Location', 'southeast');
grid on
ylabel('r (µN)^-^2 for µ_2_\rightarrow_3 = 0, a.u.');
xlabel('population size, N')

set(gca,'FontSize', 16)
saveas(gcf, 'Panel_S8e.pdf')

%% Panel S8b + extra lone for alpha = 2.5
%powerlaw (N) dependence of 1->2 transition time
figure('Color', 'w', 'Position', [1000  938    560   400]);
for li = 1 : length(List_)
    
    loglog(NN_{li}(1:9), T_{li}(1:9), 'o-'); hold on
    
end
legend(L, 'Location', 'southeast')
xlabel('population size, N');
ylabel('T_1_2(N)');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S8b.pdf')
%% Panel S8d (P_CF)
figure('Color', 'w', 'Position', [1000  938    560   400]);
for li = 1 : length(List_)-1
    semilogx(NN_{li}(1:9), SPtot_{li}(1:9), 'o-'); hold on
end
legend(L(1:end-1), 'Location', 'southeast')
xlabel('population size, N');
ylabel('P_C_F(N)');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S8d.pdf')


%% Panel S8c + extra data point.
figure('Color', 'w', 'Position', [1000  938    560   400]);
plot(1./pwr_, 1./gamma_, 'ko', 'MarkerSize', 10); hold on
xl = [0 3];
plot(xl, xl+1, 'k-');
for li = 1 : length(List_)
plot(1./pwr_(li), 1./gamma_(li), 'o', 'MarkerSize', 10); hold on
end



xlabel('\alpha^-^1');
ylabel('\gamma^-^1');
legend('simulations', '\gamma^-^1 = \alpha^-^1 + 1', 'Location', 'northwest');

set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S8c.pdf')

%% Inset of Panel S8f
figure('Color', 'w', 'Position', [ 707   470   197   134]);
plot((pwr_(2:end)-1)./(pwr_(2:end)+1), gamma1_(2:end), 'ko', 'MarkerSize', 10); hold on
xlim([0.1 0.5]);
ylim([0.1 0.5]);
xl = xlim;
plot(xl, xl, 'k-');
for li = 2 : length(List_)
plot((pwr_(li)-1)./(pwr_(li)+1), gamma1_(li), 'o', 'MarkerSize', 10); hold on
end

set(gca, 'XTick', [0.1 0.3 0.5], 'YTick', [0.1 0.3 0.5])

axis square

xlabel('(\alpha-1)/(\alpha+1)');
ylabel('\gamma');
%L=legend('simulations', '\gamma = (\alpha-1)/(\alpha+1)', 'Location', 'northwest');

set(gca, 'FontSize', 12)
saveas(gcf, 'Panel_S8f_inset.pdf')

%% Panel S8f
figure('Color', 'w', 'Position', [1000  938    560   400]);;
for li = 3 : length(List_)
loglog(NN_{li}, SPTXtot_{li}, 'o-'); hold on
end
legend(L(3:end), 'Location', 'southeast');
grid on
ylabel('r (µN)^-^2 for µ_2_\rightarrow_3 = 0, a.u.');
xlabel('population size, N')
ylim([3e6 1e9]);
set(gca,'FontSize', 16)
saveas(gcf, 'Panel_S8f.pdf')







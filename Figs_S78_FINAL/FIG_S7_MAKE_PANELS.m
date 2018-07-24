load('sym_a01_b03_v2_Px.mat')
load('sym_a01_b03_v2_SNAP.mat')


numTries = cellfun(@(x)numel(x), SNAP_); %roughly the same



SPx = {}; SPxT = {}; SPtot = []; SPTXtot = [];
Sx={};
for ni = 1 : length(N_)
    
    Sx{ni} = SNAP_{ni}/R_(ni); %(1 x N-1)
    %this is the probability of escape and snap of an existing mutant at
    %abundance 1.
    
    SPx{ni} =  Sx{ni}.*Px_{ni}/sum(Px_{ni});
    SPxT{ni} = Sx{ni}.*Px_{ni}; %this incorportes the fact that transitions take longer at large N (the log(N) factor)
%     
% 
%     
    SPtot(ni) = sum(SPx{ni});
    SPTtot(ni) = sum(SPxT{ni});
    SPTXtot(ni) = sum((1-(1:(N_(ni)-1))/N_(ni)).*SPxT{ni});
end

%% log(N) dependence of 1->2 transition time
figure('Color', 'w', 'Position', [1000  938    560   400]);
semilogx(N_, TT_/28/40000, 'o-');
xlabel('population size, N');
ylabel('T_1_2(N)');
set(gca, 'FontSize', 16)
%% Panel S7a
figure('Color', 'w', 'Position', [1000  938    560   400]);
L = plotyy(N_, TT_/28/40000, N_, 1.1e-6*SPTtot);
set(L(1), 'XScale', 'log');
set(L(2), 'XScale', 'log')
set(L(1).Children, 'Marker', 'o');
set(L(2).Children, 'Marker', '^');
L(2).FontSize = 16;
YL = get(L(2), 'Ylabel');
YL.String = 'r (µN)^-^2, a.u.';
YL.FontSize = 18;
xlabel('population size, N');
ylabel('T_1_2(N)');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S7a.pdf')
%% This should be constant if if r ~ (muN)^2 log(N)
% We see a small residual dependence, which decreases at larger N
figure('Color', 'w');
semilogx(N_, SPtot, 'o-');
xlabel('population size, N');
ylabel('P_C_F(N)');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S7b.pdf')

%% panel S7e
figure('Color', 'w', 'Position', [1000  938    560   400]);
semilogx(N_, SPTXtot, 'o-');
xlabel('population size, N');
ylabel('r (µN)^-^2 for µ_2_\rightarrow_3 = 0, a.u.');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S7e.pdf')

%% Expect these lines to converge to same integral, if r ~ (muN)^2 log(N)
L = {};
figure
for i = 1 : 8
    %the factor of N_(i) accounts for more mutations in a large population
    semilogy((1:(N_(i)-1))/N_(i), N_(i)*smooth(SPxT{i}, N_(i)/200)); hold on
    L{end+1} = num2str(N_(i));
end
legend(L)


%% PDF version of panel S7c in which all curves are apparently on top of each other
figure('Color', 'w', 'Position', [1000  938    560   400]);
for i = 1:7
    semilogy((1:(N_(i)-1))/N_(i), N_(i)*smooth(Px_{i},N_(i)/200)/40000/28); hold on
end
xlabel('x_2');
ylabel('Relative time spent at location, a.u. ');
set(gca, 'FontSize', 16)
%% Panel S7c
figure('Color', 'w', 'Position', [1000  938    560   400]);
for i = 1:2:9
    Fx = cumsum(Px_{i}/sum(Px_{i}));
    H=plot((1:(N_(i)-1))/N_(i), Fx); hold on
    COI = get(gca, 'ColorOrderIndex');
    set(gca, 'ColorOrderIndex', COI+1);
end
H.LineWidth = 2;
ylim([0 1])
xlabel('x_2');
ylabel('CDF(x_2) for strain 3 appearance');
set(gca, 'FontSize', 16)
legend('N = 10^3',...
       'N = 10^4',...
       'N = 10^5',...
       'N = 10^6',...
       'N = 10^7', 'Location', 'northwest');
saveas(gcf, 'Panel_S7c.pdf')

%% Panel_S7d. Note that it takes a LONG TIME to generate (plotting at an excessive resolution)
%The snap probability converges to the same thing, but the region of convergence is extended towards 1 at high N
figure('Color', 'w', 'Position', [1000  938    560   400]);
Xmax = [];
for i = 1:9%length(N_)
    s = smooth(Sx{i},N_(i)/200);

    x = (1:(N_(i)-1))/N_(i);
    
    x1 = 1-logspace(log10(1/N_(i)), 0, 10000);
    s1 = interp1(x,s, x1);

   plot(x1, s1); hold on
   
   s1p = polyval(polyfit(x1, s1, 4), x1);
   [~,indx] = max(s1p); 
   Xmax(i) = x1(indx);
   
end

legend('N = 1x10^3',...
       'N = 3x10^3',...
       'N = 1x10^4',...
       'N = 3x10^4',...
       'N = 1x10^5',...
       'N = 3x10^5',...
       'N = 1x10^6',...
       'N = 3x10^6',...
       'N = 1x10^7',...
       'Location', 'northwest');
xlabel('x_2');
ylabel('Prob. of comm. formation after mutation');
set(gca, 'FontSize', 16)
saveas(gcf, 'Panel_S7d.pdf')
%% Inset in Panel S7d
figure('Color', 'w', 'Position', [1000  1140   303   198]);
loglog(N_, 1-Xmax, 'o'); hold on;
loglog(N_, (1-Xmax(end)).*N_(end)./N_)
xlabel('N');
ylabel('x_1 = 1 - x_2 for max prob.');
legend('data', 'x_1 ~ 1/N fit')
axis square
xlim([1e2 1e8]);
set(gca, 'FontSize', 12, 'XTick', [1e2 1e4 1e6 1e8]); 
saveas(gcf, 'Panel_S7d_inset.pdf')

%SARP = 0123
close all
col16 = [0 0 0; 0 102 102; 153 153 255; 204 0 204; 102 102 0; 160 160 160; 204 255 153; 0 102 0; 102 255 255; 0 255 102; 255 255 0; 255 128 0; ...
    127 0 255; 0 0 153; 255 51 51; 102 0 0]/255;


%determine number of simulations
A=dir('result.mat');
A={A.name};
N_sim = length(A);







SIM_params = struct;


    %figure('Position', [ 394          96        1698        1242]);
figure('Position', [252          67        1840        1236]);
set(gcf, 'PaperUnits', 'normalized')
 set(gcf, 'PaperPosition', [0.05 0.05 0.95 0.95])
missing = [];
maxP_hit = [];

for i = 1:N_sim  %i=scan_around_QSoutput/415 Error using matlab.graphics.axis.Axes/set
%Value must be a 1x2 vector of numeric type in which the second element is larger than the first and may be Inf
    

    
   
    
    
    TTT = simulation.T;
    
    
    T1 = find(~cellfun(@isempty,simulation.X),1,'last');
    T111 = T1;
    
  
    
    
    params = ['  Cp=' num2str(simulation.Cp) '; Eff_o_p=' num2str(simulation.Eff_op) '; Deg_o_p=' num2str(simulation.Deg_op)...
        '; Ca=' num2str(simulation.Ca) '; Ca0=' num2str(simulation.Ca0) '; Cp0=' num2str(simulation.Cp0) '; maxP=' num2str(simulation.maxP) '; maxA=' num2str(simulation.maxA)];
  
    
    SIM_params(i).Cp = simulation.Cp;
    SIM_params(i).Ca = simulation.Ca;
    SIM_params(i).Eff_op = simulation.Eff_op;
    
    clf;
    
    set(gcf, 'Name', params)
    
    
    T1 = find(~cellfun(@isempty,simulation.X),1,'last');
    disp([i, T1]) 

     if T1 < 100
         continue
     end
    
   %T_ = 100 :4: T1;
    T_ = 1 :1: T1;
    b1=nan(size(T_));
    b2 = b1;
    a1 = b1;
    a2 = b1;
    Div = b1;
    
    x= simulation.X{T1}; x = x/sum(x);
    Div_end = exp(-x'*log(x))
    
    
    k = 1; TR = []; XX = [];
    for t1 = T_

        x= simulation.X{t1}; x = x/sum(x);
        
        tmp = double(simulation.Phenotypes{t1});
        tmp(tmp<-1) = (tmp(tmp<-1)+1)/4-1; %Degradation is divided by 4
       
        
        TR = [TR [tmp; t1*ones(1,size(tmp,2))]];
        XX = [XX; x];
        k = k + 1;
        
    end
    
    tmp = XX>0.001;
    
    XX1 = XX (tmp);
    P1 = TR(1,tmp);
    P2 = TR(2,tmp);
    
    P_ = TR(1:2,tmp);


    P_ = sign(P_) .* (1+log10(abs(P_)));
    
    
        

%SARP = 0123
type_ = nan(size(P_));
type_(P_ <= 0) = 0; %S
type_(P_ < -1.001) = 1; %A
type_(P_ > 0) = 2; %R
type_(P_ > (1.00+log10(2))) = 3; %P

Col_ind = 1 + type_(1,:) + 4*type_(2,:);

Col = col16(Col_ind,:);

subplot(3,3,[1 2 4 5]);
scatter3(P_(1,:), P_(2,:), TR(3,tmp), round(20 + 6*log10(XX1)), Col, 'filled'); grid on;
% try
% set(gca, 'xlim', [min([0 1.1*min(TR(1,tmp))]) max([0 1.1*max(TR(1,tmp))])]);
% set(gca, 'ylim', [min([0 1.1*min(TR(2,tmp))]) max([0 1.1*max(TR(2,tmp))])]);
% end
view(60,0)

xlabel('Phenotypes')
zlabel('time')


title(params)    
zlim([0 T_(end)]);
subplot(3,3,[3 6]);
scatter3(P_(1,:), P_(2,:), TR(3,tmp), round(20 + 6*log10(XX1)), Col, 'filled'); grid on;
% try
% set(gca, 'xlim', [min([0 1.1*min(TR(1,tmp))]) max([0 1.1*max(TR(1,tmp))])]);
% set(gca, 'ylim', [min([0 1.1*min(TR(2,tmp))]) max([0 1.1*max(TR(2,tmp))])]);
% end
view(0,90)
axis square 
try
    set(gca, 'xlim', 1.05*[min(P_(:)) max(P_(:))]);
    set(gca, 'ylim', 1.05*[min(P_(:)) max(P_(:))]);
end

xlabel('Phenotypes')
ylabel('Phenotypes')

q2=(1+log10(1+simulation.maxP));
q1 = -(1+log10(1+simulation.maxA/4));
rectangle('Position', [q1 q1 q2-q1 q2-q1]);

%    plot3(TR(1,tmp), TR(2,tmp), TR(3,tmp), '.'); grid on;

%     cmap = hsv; cmap(end,:) = [0 0 0];
%     Ct = round(TR(3,tmp)/T_(end)*64); Ct(Ct==0) = 1;
%     scatter(TR(1,tmp), TR(2,tmp), sqrt(XX(tmp)')*100, cmap(Ct,:)); axis square; grid on;
%     xlim([-simulation.maxA/4 simulation.maxP]);
%     ylim([-simulation.maxA/4 simulation.maxP]);
%     colormap(cmap); h=colorbar;

subplot(3,3,[7 8 9]);

P = TR(1:2,:);

%SARP = 0123
type_ = nan(size(P));
type_(P < 0) = 0; %S
type_(P < -1.0001) = 1; %A
type_(P > 0) = 2; %R
type_(P > 2.001) = 3; %P

Col_ind = 1 + type_(1,:) + 4*type_(2,:);

W = TR(1:end-1,:)';
[uPh, ia] = unique(W,'rows','first','legacy');
T1 = TR(end,ia);
[~, ia] = unique(W,'rows','last','legacy');
T2 = TR(end,ia);
uT=[T1;T2];
uPh = uPh';
Col_ind_ = Col_ind(ia);
dT = diff(uT);


tmp=dT > -1; %gives the option to plot only long lived lineg
%tmp=dT > 10;
uT = uT(:,tmp);
uPh = uPh(:,tmp);
Col_ind_ = Col_ind_(tmp);



for i0 = 1 : size(uPh,2);
   tmp = find(uPh(1,i0)==TR(1,:) & uPh(2,i0)==TR(2,:));
   
   semilogy(TR(3,tmp),XX(tmp), '-', 'Color', col16(Col_ind_(i0),:));
   hold on;
end


xlabel('time')
ylabel('abundances')










% %Col = col16(Col_ind,:);
% for z = 1 : 16
%     
%     tmp = find(Col_ind ==z);
%     %x=accumarray(TR(3,tmp)',XX(tmp), [T_(end) 1]);
%     if isempty(x), continue, end
%     %x(x==0) = 1e-6;
%    % semilogy(1:T_(end), x,'o-', 'Color', col16(z,:));
%     semilogy(TR(3,tmp), XX(tmp), '.', 'Color', col16(z,:));
%     hold on
% end

ylim([1e-5 1]); xlim([0 T_(end)]);
%set(gca, 'YScale', 'log')
    

end
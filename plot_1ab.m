close all
col16 = [0 0 0; 0 102 102; 153 153 255; 204 0 204; 102 102 0; 160 160 160; 204 255 153; 0 102 0; 102 255 255; 0 255 102; 255 255 0; 255 128 0; ...
    127 0 255; 0 0 153; 255 51 51; 102 0 0]/255;



figure('Position', [ 1          41        1920         963]);

set(gcf, 'PaperUnits', 'normalized')
 set(gcf, 'PaperPosition', [0.05 0.05 0.95 0.95])


    

  
    
    
    T1 = find(~cellfun(@isempty,simulation.Phenotypes),1,'last');
    T111 = T1;
    
    
    
    params = ['Cp=' num2str(simulation.Cp) '; Eff_o_p=' num2str(simulation.Eff_op) '; Deg_o_p=' num2str(simulation.Deg_op)...
        '; Ca=' num2str(simulation.Ca) '; Ca0=' num2str(simulation.Ca0) '; Cp0=' num2str(simulation.Cp0) '; maxP=' num2str(simulation.maxP) '; maxA=' num2str(simulation.maxA)];
  
   
    
    clf;
    
    set(gcf, 'Name', params)
    
    
  
    disp([i, T1]) 


    T_ = 1 :5: T1;
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
        tmp(tmp<-1) = (tmp(tmp<-1)+1)/4-1;
       
        
        TR = [TR [tmp; t1*ones(1,size(tmp,2))]];
        XX = [XX; x];
        k = k + 1;
        
    end
    
    tmp = XX>0.001;

    
    XX1 = XX (tmp);
    P1 = TR(1,tmp);
    
    P_ = TR(1,tmp);


    P_ = sign(P_) .* (1+log10(abs(P_)));
    
    
        

type_ = nan(size(P_));
type_(P_ <= 0) = 0; %S
type_(P_ < -1.001) = 1; %A
type_(P_ > 0) = 2; %R
type_(P_ > (1.00+log10(2))) = 3; %P

Col_ind = 1 + type_(1,:);

Col = col16(Col_ind,:);

subplot(3,3,[1 2 4 5]);
scatter(P_(1,:), TR(2,tmp), round(20 + 6*log10(XX1)), Col, 'filled'); grid on;

subplot(3,3,[7 8 9]);

P = TR(1,:);

%SARP = 0123
type_ = nan(size(P));
type_(P < 0) = 0; %S
type_(P < -1.0001) = 1; %A
type_(P > 0) = 2; %R
type_(P > 2.001) = 3; %P

Col_ind = 1 + type_(1,:); %+ 4*type_(2,:);

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

uT = uT(:,tmp);
uPh = uPh(:,tmp);
Col_ind_ = Col_ind_(tmp);



for i0 = 1 : size(uPh,2);
   tmp = find(uPh(1,i0)==TR(1,:));
   
   semilogy(TR(2,tmp),XX(tmp), '-', 'Color', col16(Col_ind_(i0),:));
   hold on;
end

    

ylim([1e-5 1]); xlim([0 T_(end)]);
    set(gcf, 'PaperOrientation', 'landscape')

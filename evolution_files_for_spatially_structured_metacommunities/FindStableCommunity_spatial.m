function [time]=FindStableCommunity_spatial( simulation )

    
if size(simulation.Phenotypes{1},1)==2    

    Traj=NaN;

    antibn=size(simulation.Phenotypes{1},1)+1;
    Typelarinfo=cell(1,length(simulation.Phenotypes));
    slopelar=[];
    Prod=[];zz=0;

    combinasyonlar=[[1,1];[1,2];[1,3];[1,4];[2,1];[2,2];[2,3];[2,4];[3,1];[3,2];[3,3];[3,4];[4,1];[4,2];[4,3];[4,4];];
    WQ=unique(combinasyonlar,'rows');
    Traj=zeros(size(WQ,1),length(Typelarinfo));  
    
    
    fprintf('Mergind same phenotypes...\n')
    for i=1:1:length(simulation.X)
        %%%
        zz=zz+1;
        X1=simulation.X{i};
        p=simulation.Phenotypes{i};

        Typelarinfo=zeros(size(p));
        Typelarinfo(p==-1)=3;
        Typelarinfo(p<-1)=2;
        Typelarinfo(p>2)=1;
        Typelarinfo(Typelarinfo==0)=4;

        for zxa=1:simulation.PatchNum
            X=X1(:,zxa);
            Typelarinfo(antibn,:)=X;

            WW=Typelarinfo(1:(antibn-1),:); 

            [~,a]=ismember(WW',WQ,'rows');

            TotalTypelar=WQ';

            [a,~,kq]=unique(a');
            a = a';
            Traj(a,zz,zxa) = accumarray(kq,Typelarinfo(antibn,:));
        end
    end


   QAZ=simulation.X;
  
  Tplot=Traj;
    for zxa=1:simulation.PatchNum
        time(zxa)=NaN;

      
        Traj=squeeze(Tplot(:,1:10:end,zxa));
        %identifying the community assembly.
        try
             w=3; %number of time points to be considered for identification
             time(zxa)=identification_two(Traj,w);
             
             fprintf(['The time of community assembly is predicted as: ',num2str(time),' in sub_community #' num2str(zxa) '\n']);
        catch
             fprintf('No ESC was found \n');
        end
        % visualization
        time=(time-1)*10+1;
        fprintf('The plots for visualizing the trajectory are being computed...\n')
        figure(zxa)
        
        
        for mm=1:size(Tplot,2)
        simulation.X{mm}=QAZ{mm}(:,zxa);
        end
        plot_2ab_spatial
        
    end
    
    time=min(time);



end

if size(simulation.Phenotypes{1},1)==1   

    Traj=NaN;

    antibn=size(simulation.Phenotypes{1},1)+1;
    Typelarinfo=cell(1,length(simulation.Phenotypes));
    slopelar=[];
    Prod=[];zz=0;

    combinasyonlar=[[1;2;3;4]];
    WQ=unique(combinasyonlar,'rows');
    Traj=zeros(size(WQ,1),length(Typelarinfo));  
    
    
    fprintf('Mergind same phenotypes...\n')
    for i=1:1:length(simulation.X)
        %%%
        zz=zz+1;
        X1=simulation.X{i};
        p=simulation.Phenotypes{i};

        Typelarinfo=zeros(size(p));
        Typelarinfo(p==-1)=3;
        Typelarinfo(p<-1)=2;
        Typelarinfo(p>2)=1;
        Typelarinfo(Typelarinfo==0)=4;

        for zxa=1:simulation.PatchNum
            X=X1(:,zxa);
            Typelarinfo(antibn,:)=X;

            WW=Typelarinfo(1:(antibn-1),:); 

            [~,a]=ismember(WW',WQ,'rows');

            TotalTypelar=WQ';

            [a,~,kq]=unique(a');
            a = a';
            Traj(a,zz,zxa) = accumarray(kq,Typelarinfo(antibn,:));
        end
    end


   QAZ=simulation.X;
  
  Tplot=Traj;
    for zxa=1:simulation.PatchNum
        time(zxa)=NaN;

      
        Traj=squeeze(Tplot(:,1:10:end,zxa));
        %identifying the community assembly.
        try
             w=3; %number of time points to be considered for identification
             time(zxa)=identification_two(Traj,w);
             
             fprintf(['The time of community assembly is predicted as: ',num2str(time),' in sub_community #' num2str(zxa) '\n']);
        catch
             fprintf('No ESC was found \n');
        end
        % visualization
        time=(time-1)*10+1;
        fprintf('The plots for visualizing the trajectory are being computed...\n')
        figure(zxa)
        
        
        for mm=1:size(Tplot,2)
        simulation.X{mm}=QAZ{mm}(:,zxa);
        end
        plot_1ab_spatial
        
    end
    
    time=min(time);



end
 
end
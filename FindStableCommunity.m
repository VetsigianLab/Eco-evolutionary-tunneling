function [time]=FindStableCommunity( simulation )

    
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
        zz=zz+1;
        Typelarinfo{zz}=zeros(size(simulation.Phenotypes{i}));
        Typelarinfo{zz}(simulation.Phenotypes{i}==-1)=3;
        Typelarinfo{zz}(simulation.Phenotypes{i}<-1)=2;
        Typelarinfo{zz}(simulation.Phenotypes{i}>2)=1;
        Typelarinfo{zz}( Typelarinfo{zz}==0)=4;

        Typelarinfo{zz}(antibn,:)=simulation.X{i};

        WW=Typelarinfo{zz}(1:(antibn-1),:); 

        [~,a]=ismember(WW',WQ,'rows');

        [a,~,kq]=unique(a');
        a = a';
        Traj(a,zz) = accumarray(kq,Typelarinfo{zz}(antibn,:));



    end
    time=NaN;
     
    Tplot=Traj;
    Traj=Tplot(:,1:10:end);
    %identifying the community assembly.
    try
        w=3; %number of time points to be considered for identification
         time=identification_two(Traj,w);
         time=(time-1)*10+1;
         fprintf(['The time of community assembly is predicted as: ',num2str(time),'\n']);
    catch
         fprintf('No ESC was found \n');
    end
    % visualization
    fprintf('The plots for visualizing the trajectory are being computed...\n')
    plot_2ab
    



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
        zz=zz+1;
        Typelarinfo{zz}=zeros(size(simulation.Phenotypes{i}));
        Typelarinfo{zz}(simulation.Phenotypes{i}==-1)=3;
        Typelarinfo{zz}(simulation.Phenotypes{i}<-1)=2;
        Typelarinfo{zz}(simulation.Phenotypes{i}>2)=1;
        Typelarinfo{zz}( Typelarinfo{zz}==0)=4;

        Typelarinfo{zz}(antibn,:)=simulation.X{i};

        WW=Typelarinfo{zz}(1:(antibn-1),:); 

        [~,a]=ismember(WW',WQ,'rows');

        [a,~,kq]=unique(a');
        a = a';
        Traj(a,zz) = accumarray(kq,Typelarinfo{zz}(antibn,:));



    end
    time=NaN;
     
    Tplot=Traj;
    Traj=Traj(:,1:10:end);
    %identifying the community assembly.
    try
        w=3; %number of time points to be considered for identification
         time=identification_one(Traj,w);
         time=(time-1)*10+1;
         fprintf(['The time of community assembly is predicted as: ',num2str(time),'\n']);
    catch
         fprintf('No ESC was found \n');
    end
    % visualization
    fprintf('The plots for visualizing the trajectory are being computed...\n')
    plot_1ab
 



end
 
end
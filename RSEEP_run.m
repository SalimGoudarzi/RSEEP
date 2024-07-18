function RSEEP_run
warning('off','all')
PATH                   = [pwd '\DATA\'];
%--------------------------------------------------------------------------
nSims                  = 10000;
species                = 'conifer';
% species                = 'beech';
% species                = 'grassland';
%--------------------------------------------------------------------------
%soil moisture observation to which the model is calibrated (i.e., surface soil moisture)
sensorDepth            = 0.1; 
%soil thickness [m]
zmax                   = 1.2; 
%extract soil moisture at these depths [m]
z                      = [sensorDepth 0.2 0.4 0.6 1]';
%--------------------------------------------------------------------------
%load observed data (soil moisture available at 10, 20, 40, 60 a,d 100cm depth)
load([pwd '\DATA\' species '_' 'botanic' '_observations.mat' ],'DT','rain','PET','T','TSM','obsS','SM20','SM30','SM60','SM100','pClay')
%calculate cell boundaries
zB                    = [0;(z(1:end-1)+z(2:end))/2;zmax];
%change in depth
dz                    = diff(zB);
%weights to calculate observed bulk soil moisture
wz                    = dz./zmax;
obsS                  = (wz(1)*TSM + wz(2)*SM20 + wz(3)*SM30 + wz(4)*SM60 + wz(5)*SM100)*zmax*1000;
%--------------------------------------------------------------------------
%load parameter-sets (10,000 in this case)
load([pwd '\DATA\' 'PARAMset' '_'  num2str(nSims) '.mat'],'PARAMset');
%load behavioural model IDs
load([PATH species '_' 'botanic' '_' num2str(nSims) '.mat'],'idx')
%number of behavioural models
Nb                    = length(idx);
%initialise variables
NO                    = length(TSM);
S                     = zeros(NO,1);
VSM1                  = TSM*0;
VSM2                  = TSM*0;
VSM3                  = TSM*0;
VSM4                  = TSM*0;
VSM5                  = TSM*0;
%--------------------------------------------------------------------------
%loop through behavioural models
for kk=1:Nb
    %parameter-set
    params            = PARAMset(idx(kk),:);
    %run model for it
    OUT               = RSEEP_ode(rain,PET,species,sensorDepth,params,zmax,TSM,1,1,24);
    %disaggregate
    VSM1              = VSM1 + OUT(:,1);
    VSM2              = VSM2 + OUT(:,11);
    VSM3              = VSM3 + OUT(:,12);
    VSM4              = VSM4 + OUT(:,13);
    VSM5              = VSM5 + OUT(:,14);
    S                 = S    + OUT(:,9);
    %kling-Gupta 
    [KGE(kk)]         = objective_fun_calc(TSM,OUT(:,1),0,0,0);
    %store fluxes
    OLF(kk)           = sum(OUT(:,2));
    EC(kk)            = sum(OUT(:,3));
    ES(kk)            = sum(OUT(:,4));
    TR(kk)            = sum(OUT(:,5));
    EcET(kk)          = EC(kk)./(EC(kk)+ES(kk)+TR(kk));
    EsET(kk)          = ES(kk)./(EC(kk)+ES(kk)+TR(kk));
    TrET(kk)          = TR(kk)./(EC(kk)+ES(kk)+TR(kk));
    TrPPT(kk)         = TR(kk)./(sum(rain));
    QV(kk)            = sum(OUT(:,7));
    Qb(kk)            = sum(OUT(:,6));
    ERR(kk)           = mean(OUT(:,end));
    %----------------------------------------------------------------------
    % plot
    %----------------------------------------------------------------------
    figure(44)
    %--------
    if kk==1
        clf
        tiledlayout(6,1,'TileSpacing','none');
        %--------------
        nexttile(1)
        %--------------
        yyaxis right
        hold on
        gg0=area(DT,rain,'facecolor','b','edgecolor','b','facea',0.5,'edgea',0.5);
        gg1=area(DT,PET,'facecolor','g','edgecolor','g','facea',0.90,'edgea',0.80);
        set(gca,'ydir','reverse','ycolor','k','ytick',[0 10 20 30],'xticklabels',[])
        set(gca,'color',[0.9 0.9 0.9]);
        ylabel('PPT/PET [mm/day]')
        set(gca,'ycolor','b')
        xlim([DT(1) DT(end)])
        ylim([0 50])
    end
    %--------------
    ax=nexttile(1);
    %--------------
    yyaxis left
    hold on;box on;
    area(DT,OUT(:,1),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none','linewidth',2)
    set(gca,'xticklabels',[])
    ax=nexttile(2);
    hold on;box on;
    area(DT,OUT(:,11),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none','linewidth',2)
    set(gca,'xticklabels',[])
    %--------------
    ax=nexttile(3);
    %--------------
    hold on;box on;
    area(DT,OUT(:,12),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none','linewidth',2)
    set(gca,'xticklabels',[])
    %--------------
    ax=nexttile(4);
    %--------------
    hold on;box on;
    area(DT,OUT(:,13),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none','linewidth',2)
    set(gca,'xticklabels',[])
    %--------------
    ax=nexttile(5);
    %--------------
    hold on;box on;
    area(DT,OUT(:,14),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none','linewidth',2)
    set(gca,'xticklabels',[])
    %--------------
    ax=nexttile(6);
    %--------------
    hold on; box on;
    area(DT,OUT(:,9),'edgecolor','#4DBEEE','edgea',0.5,'facecolor','none')
    if kk==Nb
        %--------------
        figure(44)
        %--------------
        ax=nexttile(1);
        %--------------
        yyaxis left
        hold on;box on;
        hh0=plot(DT,TSM,'color','#EDB120','linewidth',2,'marker','none','linestyle','-');
        hh1=plot(DT,TSM*nan,'-','color','#4DBEEE','marker','none');
        hh2=plot(DT,VSM1/Nb,':','color','#7E2F8E','marker','none','linewidth',2);
        legend([gg0 gg1 hh0 hh1 hh2],'PPT','PET','Measurement','Individual models','Mean of models','location','best','orientation','horizontal')
        title([species  '  |  ' location '  |  ' '# models = '  num2str(Nb) ' (top 1%)' '  |  ' 'KGE = ' num2str(round(min(KGE),2)) ' : ' num2str(round(max(KGE),2))])
        set(gca,'ycolor','k')
        set(gca,'ytick',[0.2 0.4 0.6 0.8 1])
        xlim([DT(1) DT(end)])
        ylim([0 1])
        ylabel('VMC [-]')
        text(DT(floor(NO/2)-7),0.6,['@ - ' num2str(z(1)*100) 'cm'])
        text(DT(floor(NO-30)),0.8,'(a)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
        %--------------
        ax=nexttile(2);
        %--------------
        hold on;box on;
        hh2=plot(DT,SM20,'-','color','#EDB120','marker','none','linewidth',2);
        plot(DT,VSM2/Nb,':','color','#7E2F8E','marker','none','linewidth',2);
        set(gca,'ytick',[0.2 0.4 0.6 0.8])
        set(gca,'xticklabels',[])
        xlim([DT(1) DT(end)])
        ylim([0 1])
        ylabel('VMC [-]')
        text(DT(floor(NO/2)-7),0.6,['@ -' num2str(round(z(2)*100,0)) 'cm'])
        text(DT(floor(NO-30)),0.8,'(b)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
        %--------------
        ax=nexttile(3);
        %--------------
        hold on;box on;
        hh4=plot(DT,SM30,'-','color','#EDB120','marker','none','linewidth',2);
        plot(DT,VSM3/Nb,':','color','#7E2F8E','marker','none','linewidth',2);
        set(gca,'ytick',[0.2 0.4 0.6 0.8])
        set(gca,'xticklabels',[])
        xlim([DT(1) DT(end)])
        ylim([0 1])
        ylabel('VMC [-]')
        text(DT(floor(NO/2)-7),0.6,['@ -' num2str(round(z(3)*100,0)) 'cm'])
        text(DT(floor(NO-30)),0.8,'(c)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
        %--------------
        ax=nexttile(4);
        %--------------
        hold on;box on;
        hh6=plot(DT,SM60,'-','color','#EDB120','marker','none','linewidth',2);
        plot(DT,VSM4/Nb,':','color','#7E2F8E','marker','none','linewidth',2);
        set(gca,'ytick',[0.2 0.4 0.6 0.8])
        set(gca,'xticklabels',[])
        xlim([DT(1) DT(end)])
        ylim([0 1])
        ylabel('VMC [-]')
        text(DT(floor(NO/2)-7),0.6,['@ -' num2str(round(z(4)*100,0)) 'cm'])
        text(DT(floor(NO-30)),0.8,'(d)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
        %--------------
        ax=nexttile(5);
        %--------------
        hold on;box on;
        hh8=plot(DT,SM100,'-','color','#EDB120','marker','none','linewidth',2);
        plot(DT,VSM5/Nb,':','color','#7E2F8E','marker','none','linewidth',2);
        set(gca,'ytick',[ 0.2 0.4 0.6 0.8])
        set(gca,'xticklabels',[])
        xlim([DT(1) DT(end)])
        ylim([0 1])
        ylabel('VMC [-]')
        text(DT(floor(NO/2)-7),0.6,['@ -' num2str(round(z(5)*100,0)) 'cm'])
        text(DT(floor(NO-30)),0.8,'(e)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
        %--------------
        ax=nexttile(6);
        %--------------
        hold on; box on;
        ylabel('Storage [mm]')
        hhh0=plot(DT,obsS,'color','#EDB120','linewidth',2,'marker','none','linestyle','-'); 
        hhh1=plot(DT,S/Nb*nan,'-','marker','none','color','#4DBEEE');
        hhh2=plot(DT,S/Nb,':','marker','none','color','#7E2F8E','linewidth',2);
        legend([hhh0 hhh1,hhh2],'Measurement-based estimate','Individual models','Mean of models','location','southoutside','orientation','horizontal')
        set(gca,'ycolor','k','ytick',[0 200 400 600 800 1000])
        ylim([0 zmax*1000])
        xlim([DT(1) DT(end)])
        text(DT(floor(NO/2)-20),0.8*zmax*1000,['Bulk storage (soil thickness=' num2str(round(zmax,2)) 'm)'])
        text(DT(floor(NO-30)),0.8*zmax*1000,'(f)')
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        ax.GridLineStyle = ':';
        ax.GridColor='k';
        set(gca, 'FontWeight', 'Bold')
    end
end
function parallel_MASTER
% clc
warning('off','all')
PATH                          = pwd;
FOLDER                        = '\DATA\';
%--------------------------------------------------------------------------
p = gcp('nocreate');
if isempty(p)
    %check if a  pool is open
    defaultProfile = parallel.defaultClusterProfile;
    %get cluster profile
    myCluster      = parcluster(defaultProfile);
    %close it
    delete(myCluster.Jobs)
    %get cluster profile
    c              = parcluster();
    %set path to matlab default
    cd(c.JobStorageLocation);
    %start parpool
    parpool(7)
    %redefine storage location
    c.JobStorageLocation = [PATH FOLDER];
    %change current directory
    cd(PATH)
end
%--------------------------------------------------------------------------
nSims                       = 10000;
location                    = {'botanic'};
species                     = {'grassland','beech'};
% species                     = {'conifer'};
% %--------------------------------------------------------------------------
% nSims                       = 10001;
% location                    = {'middle'};
% species                     = {'larch','sycamore','pine','grass'};
% species                     = {'grass'};
% year                        = 2022;
%--------------------------------------------------------------------------
%start of the record
for il=1:length(location)
    %-------------------------------
    if strcmp(location{il},'botanic')
        sensorDepth         = 0.1;
        zmax                = 1.2;
        age                 = 100;
    else
        sensorDepth         = 0.06;
        zmax                = 0.5;
        age                 = 34;
    end
    %-------------------------------
    for is=1:length(species)
        load([pwd '\DATA\' species{is} '_' location{il} '_LOA.mat' ],'TSM','rain','PET','DT')
        %---------------------------------
        N                   = length(rain);
        %---------------------------------
        %filename
        filename            = ['EcoHydroPlot_' species{is} '_' location{il} '_' num2str(nSims) ];
        %--------------------------------------------------------------------------
        %                             parameter-set
        %--------------------------------------------------------------------------
        %parameter-sets
        load([pwd '\DATA\' 'PARAMset' '_'  num2str(nSims) '.mat'],'PARAMset')
        %--------------------------------------------------------------------------
        %                                run model
        %--------------------------------------------------------------------------
        %initialise arrays
        KGE                 = nan(N,1);
        LOA                 = nan(N,1);
        %------------------------------------------------------------------
        %run parallel
        parfor iSim         = 1:nSims
            [LAI,LAIg]      = seasonal_sinewave_LAI_FFT(PET,species{is},age);
            %GW recharge limit mod factor relative to grass
            qvMOD           = 1-sum(LAI-LAIg)./sum(LAI+LAIg);
            %allocate parameter
            param           = PARAMset(iSim,:);
            %run
            out             = EcoHydroPlot_ode(rain,PET,species{is},sensorDepth,param,zmax,TSM,1,1,age);
            vsm             = out(:,1);
            qvv             = sum(out(:,7));
            %GW recharge min and max
            minQV           = qvMOD*0.45*sum(rain);
            maxQV           = qvMOD*0.65*sum(rain);
                        %performance metrics
            [KGE(iSim),LOA(iSim)]= objective_fun_calc(TSM,vsm,qvv,minQV,maxQV);
        end
        %lielihood weights
        [W,Nb]              = LikelihoodWeightCal(KGE,LOA);
        %weights
        idx                 = find(W>sqrt(eps));
        idm                 = find(W==max(W));
        idm                 = idm(1);
        %         %save objective functions
        %         save([PATH FOLDER filename  '_ofs.mat'],'KGE','LOA','W','idx','idm');
        %------------------------------------------------------------------
        VSM                 = nan(Nb,N);
        Ec                  = nan(Nb,N);
        OLF                 = nan(Nb,N);
        Es                  = nan(Nb,N);
        Tr1                 = nan(Nb,N);
        Tr2                 = nan(Nb,N);
        Qv                  = nan(Nb,N);
        SMD                 = nan(Nb,N);
        %         %loop again to save the variables this time
        %         load([PATH FOLDER filename  '_ofs.mat'],'W');
        parfor iSim         = 1:Nb
            %run mode
            %             if ismember(iSim,idb)
            %allocate parameter
            param       = PARAMset(idx(iSim),:);
            OUT         = EcoHydroPlot_ode(rain,PET,species{is},sensorDepth,param,zmax,TSM,1,1,age);
            VSM(iSim,:) = OUT(:,1);
            OLF(iSim,:) = OUT(:,2);
            Ec(iSim,:)  = OUT(:,3);
            Es(iSim,:)  = OUT(:,4);
            Tr1(iSim,:) = OUT(:,5);
            Tr2(iSim,:) = OUT(:,6);
            Qv(iSim,:)  = OUT(:,7);
            SMD(iSim,:) = OUT(:,end);
            %             end
        end
        mtemp0             = median(SMD,1,'omitnan');
        mtemp0             = median(mtemp0);
        mtemp              = median(SMD,2,'omitnan')';
        ofs                = abs(mtemp-mtemp0)';
        id_med             = find(ofs==min(ofs),1);
        
        mtemp0             = min(SMD,[],1,'omitnan');
        mtemp0             = min(mtemp0);
        mtemp              = min(SMD,[],2,'omitnan')';
        ofs                = abs(mtemp-mtemp0)';
        id_min             = find(ofs==min(ofs),1);
        
        mtemp0             = max(SMD,[],1,'omitnan');
        mtemp0             = max(mtemp0);
        mtemp              = max(SMD,[],2,'omitnan')';
        ofs                = abs(mtemp-mtemp0)';
        id_max             = find(ofs==min(ofs),1);
        %save objective functions
        save([PATH FOLDER filename  '_ofs.mat'],'KGE','LOA','W','idx','idm','id_med','id_min','id_max');
        %save discharge
        save([PATH FOLDER filename  '.mat'    ],'VSM','Ec','OLF','Es','Tr1','Tr2','Qv');
        %------------------------------------------------------------------
        disp(filename)
    end
end

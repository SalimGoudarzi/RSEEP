function RothC
clc

PATH                          = pwd;
FOLDER                        = '\DATA\';
% %--------------------------------------------------------------------------
% %check if a  pool is open
% defaultProfile = parallel.defaultClusterProfile;
% %get cluster profile
% myCluster      = parcluster(defaultProfile);
% %close it
% delete(myCluster.Jobs)
% %get cluster profile
% c              = parcluster();
% %set path to matlab default
% cd(c.JobStorageLocation);
% %start parpool
% parpool(7)
% %redefine storage location
% c.JobStorageLocation = [PATH FOLDER];
% %change current directory
% cd(PATH)
%--------------------------------------------------------------------------
location    = 'middle';
species     = 'pine';


% load([pwd '\DATA\' 'PARAMset' '_'  num2str(10001) '.mat'],'PARAMset')
% load([pwd '\DATA\' species '_' location '_' num2str(10001) '.mat'],'idx')
% load([pwd '\DATA\' species '_' location '_LOA.mat' ],'maxTSM')
% OUT      = EcoHydroPlot_ode(rain,PET,species,sensorDepth,param,zmax,maxTSM);
% SMD      = OUT(:,end);
% %(convert to mm)      1     2     3    4     5    6    7     8    9   10     11   12   13   14   15             16   
% %OUT               = [VSM1  LPOF   Ec   Es    Tr   Qb   Qv    Sc   S   err   VSM2 VSM3 VSM4 VSM5 zmax*ones(N,1) SMD]*1000


pClay       = 4.48;
if strcmp(species,'grass')
    TOC1997 = 72.3;
    TOC2012 = 68;
    %max height per year [cm]
    dmax        = 20;
elseif strcmp(species,'larch')
    TOC1997 = 78;
    TOC2012 = 79;
    %max tree diamter [cm]
    dmax        = 70;
elseif strcmp(species,'pine')
    TOC1997 = 72.9;
    TOC2012 = 75;
    %max tree diamter [cm]
    dmax        = 70;
elseif strcmp(species,'sycamore')
    TOC1997 = 76.2;
    TOC2012 = 70;
    %max tree diamter [cm]
    dmax        = 70;
end
%soil thickness [m]
zmax        = 0.5;

[DT1,SMD1,DPM1,RPM1,BIO1,HUM1,TOC1,CO21,INPUT1,DT2,SMD2,DPM2,RPM2,BIO2,HUM2,TOC2,CO22,INPUT2...
    ,DT3,SMD3,DPM3,RPM3,BIO3,HUM3,TOC3,CO23,INPUT3,DT4,SMD4,DPM4,RPM4,BIO4,HUM4,TOC4,CO24,INPUT4,AB,BL]=RothC_ode(species,TOC1997,TOC2012,zmax,pClay,dmax,location);

%--------------------------------------------------------------------------
%                             plot results
%--------------------------------------------------------------------------
LW          = 1.5;
if strcmp(species,'pine')
    figure(35)
elseif strcmp(species,'larch')
    figure(36)
elseif strcmp(species,'sycamore')
    figure(37)
elseif strcmp(species,'grass')
    figure(38)
end
clf
tiledlayout(3,1,'TileSpacing','tight');
%------------------
ax(1)=nexttile(1);
hold on; box on;
hh(1)=patchPlot(DT3,SMD3,SMD4,'b',0.25);
% hh(2)=plot(DT2,SMD2,'color','#0072BD','linestyle','--');
hh(3)=plot(DT1,SMD1,'color','k','linestyle','-');

ylabel('[mm]')
ylim([-10 70])
if strcmp(species,'grass')
    title([species ' | growth/yr = ' num2str(dmax) ' cm' ])
else
    title([species ' | DBH @yr200 = ' num2str(dmax) ' cm' ])
end
legend([hh(1) hh(3)],'SDM:future (rcp6)','SMD:past data looped')

ax(1).YGrid = 'on';
ax(1).XGrid = 'on';
ax(1).GridLineStyle = ':';
ax(1).GridColor='k';
set(gca, 'FontWeight', 'Bold')

%----------------
ax(2)=nexttile(2);
hold on; box on;

hh(1)=patchPlot(DT3,DPM3,DPM4,'r',0.25);
% hh(2)=plot(DT2,DPM2,'-','marker','none','color','#A2142F','linewidth',LW,'linestyle','--');
hh(3)=plot(DT1,DPM1,'-','marker','none','color','#A2142F','linewidth',LW,'linestyle','-');


hh(4)=patchPlot(DT3,BIO3,BIO4,'g',0.25);
% hh(5)=plot(DT2,BIO2,'-','marker','none','color','#77AC30','linewidth',LW,'linestyle','--');
hh(6)=plot(DT1,BIO1,'-','marker','none','color','#77AC30','linewidth',LW,'linestyle','-');


set(gca,'ycolor','k','xticklabel',[])
ylabel('[tC/ha]')
ylim([0 Inf])

legend([hh(1) hh(3) hh(4) hh(6)],'DPM:future (rcp6)','DPM:past data looped','BIO:future (rcp6)','BIO:past data looped','location','east')

ax(2).YGrid = 'on';
ax(2).XGrid = 'on';
ax(2).GridLineStyle = ':';
ax(2).GridColor='k';
set(gca, 'FontWeight', 'Bold')
%------------------
ax(3)=nexttile(3);
hold on; box on;

hh(1)=patchPlot(DT3,RPM3,RPM4,'c',0.25);
% hh(2)=plot(DT2,RPM2,'-','marker','none','color','#4DBEEE','linewidth',LW,'linestyle','--');
hh(3)=plot(DT1,RPM1,'-','marker','none','color','#4DBEEE','linewidth',LW,'linestyle','-');

hh(4)=patchPlot(DT3,HUM3,HUM4,[0.4941    0.1843    0.5569],0.25);
% hh(5)=plot(DT2,HUM2,'-','marker','none','color','#7E2F8E','linewidth',LW,'linestyle','--');
hh(6)=plot(DT1,HUM1,'-','marker','none','color','#7E2F8E','linewidth',LW,'linestyle','-');


set(gca,'ycolor','k','xticklabel',[])
ylabel('[tC/ha]')
ylim([0 Inf])
legend([hh(1) hh(3) hh(4) hh(6) ],'RPM:future (rcp6)','RPM:past data looped','HUM:future (rcp6)','HUM:past data looped','location','southeast')


ax(3).YGrid = 'on';
ax(3).XGrid = 'on';
ax(3).GridLineStyle = ':';
ax(3).GridColor='k';
set(gca, 'FontWeight', 'Bold')
%------------------
linkaxes(ax(1:end),'x')


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
perYr       = 360;
LW          = 1.5;
if strcmp(species,'pine')
    figure(45)
elseif strcmp(species,'larch')
    figure(46)
elseif strcmp(species,'sycamore')
    figure(47)
elseif strcmp(species,'grass')
    figure(48)
end
clf
tiledlayout(2,1,'TileSpacing','tight');
%------------------
ax(1)=nexttile(1);
hold on; box on;


hh(1)=patchPlot(DT3,TOC3,TOC4,[ 0.9294    0.6941    0.1255],0.25);
% hh(2)=plot(DT2,TOC2,'-','marker','none','color','#EDB120','linewidth',LW,'linestyle','--');
hh(3)=plot(DT1,TOC1,'-','marker','none','color','#EDB120','linewidth',LW,'linestyle','-');


hh(4)=patchPlot(DT3,cumsum(CO23*perYr),cumsum(CO24*perYr),'r',0.25);
% hh(5)=plot(DT2,cumsum(CO22*perYr),'-','marker','none','color','#D95319','linewidth',LW,'linestyle','--');
hh(6)=plot(DT1,cumsum(CO21*perYr),'-','marker','none','color','#D95319','linewidth',LW,'linestyle','-');

hh(7)=plot(DT2,AB,'-','marker','none','color','g','linewidth',LW,'linestyle','-');
hh(8)=plot(DT2,BL,'-','marker','none','color','b','linewidth',LW,'linestyle','-');

ylabel('[tC/ha]')
ylim([0 Inf])

legend([hh(1) hh(3) hh(4) hh(6) hh(7:8)],'Soil C:future (rcp6)','Soil C:past data looped'...
    ,'Emitted CO2:future (rcp6)','Emitted CO2:past data looped','Above-ground biomass','Below-ground biomass','location','northwest')

if strcmp(species,'grass')
    title([species ' | growth/yr = ' num2str(dmax) ' cm' ])
else
    title([species ' | DBH @yr200 = ' num2str(dmax) ' cm' ])
end
ax(1).YGrid = 'on';
ax(1).XGrid = 'on';
ax(1).GridLineStyle = ':';
ax(1).GridColor='k';
set(gca, 'FontWeight', 'Bold')


%----------------
ax(2)=nexttile(2);
hold on; box on;


hh(1)=patchPlot(DT3,TOC3+AB+BL-cumsum(CO23*perYr),TOC4+AB+BL-cumsum(CO24*perYr),'k',0.25);
% hh(2)=plot(DT2,TOC2+AB+BL-cumsum(CO22*perYr),'-','marker','none','color','k','linewidth',LW,'linestyle','--');
hh(3)=plot(DT1,TOC1+AB+BL-cumsum(CO21*perYr),'-','marker','none','color','k','linewidth',LW,'linestyle','-');



legend([hh(1) hh(3)],'Net C:future (rcp6)','Net C:past data looped','location','northwest')


ylabel('[tC/ha]')
% ylim([0 Inf])


ax(2).YGrid = 'on';
ax(2).XGrid = 'on';
ax(2).GridLineStyle = ':';
ax(2).GridColor='k';
set(gca, 'FontWeight', 'Bold')
%------------------
linkaxes(ax(1:end),'x')
%**************************************************************************
function [DT1,SMDo1,DPMo1,RPMo1,BIOo1,HUMo1,TOCo1,CO2o1,plantCo1,DT2,SMDo2,DPMo2,RPMo2,BIOo2,HUMo2,TOCo2,CO2o2,plantCo2...
    ,DT3,SMDo3,DPMo3,RPMo3,BIOo3,HUMo3,TOCo3,CO2o3,plantCo3,DT4,SMDo4,DPMo4,RPMo4,BIOo4,HUMo4,TOCo4,CO2o4,plantCo4,AB,BL]=RothC_ode(species,TOC1997,TOC2012,zmax,pClay,dmax,location)
%--------------------------------------------------------------------------
%                               parameters
%--------------------------------------------------------------------------
%R:DPM/RPM ratio [-] for:
%   -agri crop/improved grassland: 1.44
%   -un-improved grassland:        0.67
%   -dedicious/tropical woodland:  0.25
if strcmp(species,'grass')
    R            = 0.67;
else
    R            = 0.25;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                     Equilibrium run (for 1997)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%equilibrium run for xxx number of years
nYrsEQ           = 10000;

[plantRe,raine,PETe,Te,agee,LAIe]=read_wheather_data('equilibrium',species,[]);
%replicate to make 10000 yrs
Nrep             = ceil(nYrsEQ/(length(raine)/360));
LAIe             = repmat(LAIe,Nrep,1);
plantRe          = repmat(plantRe,Nrep,1);
raine            = repmat(raine,Nrep,1);
PETe             = repmat(PETe,Nrep,1);
Te               = repmat(Te,Nrep,1);
%soil cover
cover            = LAIe>min(LAIe)*1.01;
manure           = cover*0;
%equilibrium run
[SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,EQinput,plantCE]=initialise(TOC1997,R,zmax,pClay,raine,Te,PETe,cover,manure,plantRe);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                Forward run until the observation period (2012)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[plantRe,raine,PETe,Te,agee]=read_wheather_data('obervation',species,[]);

if strcmp(species,'grass')
    maxInput     = Inf;
    inputs       = plant_input(EQinput,maxInput,agee,species,dmax);
    plantC       = plantRe.*inputs;
    %initialise variables
    SMD          = zeros(size(PETe)); SMD(1)  = SMDE(end);
    DPM          = zeros(size(PETe)); DPM(1)  = DPME(end);
    RPM          = zeros(size(PETe)); RPM(1)  = RPME(end);
    BIO          = zeros(size(PETe)); BIO(1)  = BIOE(end);
    HUM          = zeros(size(PETe)); HUM(1)  = HUME(end);
    %run model
    [SMD,DPM,RPM,BIO,HUM,TOC,CO2]=carbon_turnover_cal(SMD,DPM,RPM,BIO,HUM,plantC,manure,raine,PETe,Te,cover,R,zmax,pClay);
else
    simTOC       = Inf;
    IOM          = 0.049*TOC2012.^(1.139);
    maxInput     = EQinput;
    while abs(simTOC-(TOC2012-IOM))./(TOC2012-IOM)*100>1%to within the standard error of the Beckert et al. values
        inputs   = plant_input(EQinput,maxInput,agee,species,dmax);
        plantC   = plantRe.*inputs;
        %initialise variables
        SMD      = zeros(size(PETe)); SMD(1)  = SMDE(end);
        DPM      = zeros(size(PETe)); DPM(1)  = DPME(end);
        RPM      = zeros(size(PETe)); RPM(1)  = RPME(end);
        BIO      = zeros(size(PETe)); BIO(1)  = BIOE(end);
        HUM      = zeros(size(PETe)); HUM(1)  = HUME(end);
        %run model
        [SMD,DPM,RPM,BIO,HUM,TOC,CO2]=carbon_turnover_cal(SMD,DPM,RPM,BIO,HUM,plantC,manure,raine,PETe,Te,cover,R,zmax,pClay);
        simTOC   = TOC(end);
        %         [simTOC (TOC2012-IOM)]
        maxInput = maxInput*1.005;
    end
end
[EQinput maxInput]
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                    Forward run with past data looped 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[DT1,SMDo1,DPMo1,RPMo1,BIOo1,HUMo1,TOCo1,CO2o1,plantCo1]=run_future('futureFpast',species...
    ,[], EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                    Forward run with future data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[DT2,SMDoh,DPMoh,RPMoh,BIOoh,HUMoh,TOCoh,CO2oh,plantCoh]=run_future('futureFchess',species...
    ,'01', EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax);

[DT3,SMDoi,DPMoi,RPMoi,BIOoi,HUMoi,TOCoi,CO2oi,plantCoi]=run_future('futureFchess',species...
    ,'04', EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax);

[DT4,SMDoj,DPMoj,RPMoj,BIOoj,HUMoj,TOCoj,CO2oj,plantCoj]=run_future('futureFchess',species...
    ,'06', EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax);

[DT5,SMDok,DPMok,RPMok,BIOok,HUMok,TOCok,CO2ok,plantCok]=run_future('futureFchess',species...
    ,'15', EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax);


SMDo2    = median([SMDoh SMDoi SMDoj SMDok],2);
DPMo2    = median([DPMoh DPMoi DPMoj DPMok],2);
RPMo2    = median([RPMoh RPMoi RPMoj RPMok],2);
BIOo2    = median([BIOoh BIOoi BIOoj BIOok],2);
HUMo2    = median([HUMoh HUMoi HUMoj HUMok],2);
TOCo2    = median([TOCoh TOCoi TOCoj TOCok],2);
CO2o2    = median([CO2oh CO2oi CO2oj CO2ok],2);
plantCo2 = median([plantCoh plantCoi plantCoj plantCok],2);


SMDo3    = min([SMDoh SMDoi SMDoj SMDok],[],2);
DPMo3    = min([DPMoh DPMoi DPMoj DPMok],[],2);
RPMo3    = min([RPMoh RPMoi RPMoj RPMok],[],2);
BIOo3    = min([BIOoh BIOoi BIOoj BIOok],[],2);
HUMo3    = min([HUMoh HUMoi HUMoj HUMok],[],2);
TOCo3    = min([TOCoh TOCoi TOCoj TOCok],[],2);
CO2o3    = min([CO2oh CO2oi CO2oj CO2ok],[],2);
plantCo3 = min([plantCoh plantCoi plantCoj plantCok],2);


SMDo4    = max([SMDoh SMDoi SMDoj SMDok],[],2);
DPMo4    = max([DPMoh DPMoi DPMoj DPMok],[],2);
RPMo4    = max([RPMoh RPMoi RPMoj RPMok],[],2);
BIOo4    = max([BIOoh BIOoi BIOoj BIOok],[],2);
HUMo4    = max([HUMoh HUMoi HUMoj HUMok],[],2);
TOCo4    = max([TOCoh TOCoi TOCoj TOCok],[],2);
CO2o4    = max([CO2oh CO2oi CO2oj CO2ok],[],2);
plantCo4 = max([plantCoh plantCoi plantCoj plantCok],2);


%above and below ground biomass
AGE0     = datetime(1998,1,1,0,0,0);
AGE      = DT2>AGE0;
AGE      = cumsum(AGE);


[AB,BL]  = biomass(species,AGE,dmax);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                             sub-funcitons
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [SMDi,DPMi,RPMi,BIOi,HUMi,TOCi,CO2i,annualInput,plantC]=initialise(targetTOC,R,H,pClay,rain,T,PET,cover,manure,plantR)
%number of timesteps
N                    = length(PET);
%IOM estimate from total organic carbon (TOC)
IOM                  = 0.049*targetTOC.^(1.139);
%deduct IOM from target for equilibrium run
targetTOC            = targetTOC - IOM;
%plant input carbon initial guess [tC/ha/month]
annualInput          = eps;
%initial error [%]
error                = 100;
%continue while error is above threshold
simTOC               = targetTOC;
while error>1%
    annualInput      = annualInput*(targetTOC)./(simTOC +eps);
    %annual average plantC input (1st guess) [tC/ha/month]
    plantC           =  annualInput.*plantR;
    %top-soil moisture deficit
    SMD0             = zeros(N,1);
    %decomposable plant material
    DPM0             = zeros(N,1);
    %resistant plant material
    RPM0             = zeros(N,1);
    %biomass
    BIO0             = zeros(N,1);
    %humus
    HUM0             = zeros(N,1);
    %load inputs
    [SMDi,DPMi,RPMi,BIOi,HUMi,TOCi,CO2i]=carbon_turnover_cal(SMD0,DPM0,RPM0,BIO0,HUM0,plantC,manure,rain,PET,T,cover,R,H,pClay);
    %total organic carbon at the end of simulation [tC/ha]
    simTOC           = TOCi(end);
    error            = abs(simTOC-targetTOC)./(targetTOC)*100;
end
%**************************************************************************
function [SMD,DPM,RPM,BIO,HUM,TOC,CO2,INPUT]=carbon_turnover_cal(SMD,DPM,RPM,BIO,HUM,plantC,manureC,rain,PET,T,cover,R,H,pClay)
%initialse CO2 release vector
CO2                  = zeros(size(DPM));
INPUT                = zeros(size(DPM));
INPUT(1)             = plantC(1);
day                  = 1;
for ii = 2:length(PET)
    %update month number
    day              = day + 1;
    %----------------------------------------------------------------------
    %               rate modifying factor for temperature
    %----------------------------------------------------------------------
    a                = 47.91./(1 + exp(106.06./(T(day) + 18.27)));
    %----------------------------------------------------------------------
    %                 top-soil moisture deficit (TSMD)
    %----------------------------------------------------------------------
    SMD(ii)          = SMD(ii-1) + 0.75*PET(day) - rain(day);
    %deficit cannot go negative
    if SMD(ii)<0
        SMD(ii)      = 0;
    end
    %max top-soil moisture deficit (TSMD) for the 0-23cm of a soil
    TSMD_max         = 20 + 1.3*pClay - 0.01*pClay.^2;
    %adjust if top soil thickness is not 23cm
    TSMD_max         = TSMD_max./0.23*H;
    %soil cover rate modifier for vegetated cover
    c                = 0.6;
    if cover(day)==0
        TSMD_max     = max(TSMD_max/1.8,SMD(ii-1));
        c            = 1;
    end
    %soil moisture deficit cannot exceed max
    if SMD(ii)>TSMD_max
        SMD(ii)     = TSMD_max;
    end
    %rate modifying factor for top-soil moisture deficit
    if SMD(ii)<0.444*TSMD_max
        b            = 1;
    else
        b            =  0.2 + (1-0.2)*(TSMD_max-SMD(ii))./(TSMD_max-0.444*TSMD_max);
    end
    %----------------------------------------------------------------------
    %                              decay rates
    %----------------------------------------------------------------------
    t                = 1./365; %daily
    dDPM             = DPM(ii-1).*( 1 - exp(- a.*b.*c*10.0*t ) );
    dRPM             = RPM(ii-1).*( 1 - exp(- a.*b.*c*0.30*t ) );
    dBIO             = BIO(ii-1).*( 1 - exp(- a.*b.*c*0.66*t ) );
    dHUM             = HUM(ii-1).*( 1 - exp(- a.*b.*c*0.02*t ) );
    %total decayed material
    TOTAL            = dDPM + dRPM + dBIO + dHUM;
    %----------------------------------------------------------------------
    %                             partitioning
    %----------------------------------------------------------------------
    %partitioning
    part             = 1.67*(1.85 + 1.60*exp( - 0.0786*pClay ));
    %fraction into CO2
    X_CO2            = part./(1 + part);
    %fraction of BIO+HUM that is BIO
    X_BIO            = 0.46*(1-X_CO2);
    %fraction of BIO+HUM that is HUM
    X_HUM            = 0.54*(1-X_CO2);
    %----------------------------------------------------------------------
    %update values
    DPM(ii)          = DPM(ii-1)  - dDPM + manureC(day)*0.49 + plantC(day).*R./(R+1);
    RPM(ii)          = RPM(ii-1)  - dRPM + manureC(day)*0.49 + plantC(day)./(R+1);
    BIO(ii)          = BIO(ii-1)  - dBIO + TOTAL.*X_BIO;
    HUM(ii)          = HUM(ii-1)  - dHUM + TOTAL.*X_HUM + manureC(day)*0.02;
    %     CO2(ii)          = CO2(ii-1)  + X_CO2*TOTAL;
    CO2(ii)          = X_CO2*TOTAL;
    INPUT(ii)        = plantC(day);
end
%total organic carbon
TOC                  = DPM + RPM + BIO + HUM;
%**************************************************************************
function Input       = plant_input(EQ,MAX,t,species,dmax)
%X: time in yrs
%Input: is also per year
alpha                = 1;
t0                   = 5;
if strcmp(species,'grass')
    Input            = ones(size(t))*EQ;
else
    Y                = biomass(species,t,dmax);
    YN               = Y./max(Y);
    Input            = YN*(MAX-EQ)+EQ;
    Input(t<=t0)     = -alpha*EQ+alpha*EQ/t0.*t(t<=t0)+EQ;
end
%**************************************************************************
function [AB,BL]     = biomass(species,t,dmax)

switch species
    case 'sycamore'
        a            = -2.455;
        b            = 2.354;
    case 'larch'
        a            = -2.26;
        b            = 2.298;
    case 'pine'
        a            = -2.029;
        b            = 2.289;
    case 'grass'
        a0            = 0.0174;
        b0            = 0.5242;       
        a1            = 0.0185;
        b1            = 0.9377;
end

if strcmp(species,'grass')
    %from 2018, Bathelomee et al.,
    %total biomass in tC/ha per year
    TOT               = exp(dmax*a1+b1)* t;
    AB                = exp(dmax*a0+b0)* t;
    BL                = TOT-AB;
else
    
    beta                 = 0.02;
    t0                   = 5;
    DBH                  = dmax*(1-exp(-beta.*[t-t0]));
    DBH(DBH<0)           = 0;
    %abive ground biomass (kg)
    AB                   = exp(a+b*log(DBH));
    %below ground biomass (roots)
    BL                   = exp(-1.0587 + 0.8836*log(AB+0.2840));
    %convert to tC/ha (400 tree/ha)
    AB                   = AB/1000*400;
    BL                   = BL/1000*400;
end
%**************************************************************************
function [DT1,SMDo1,DPMo1,RPMo1,BIOo1,HUMo1,TOCo1,CO2o1,plantCo1]=run_future(scenario,species,ens...
    , EQinput,maxInput,SMD,DPM,RPM,BIO,HUM,TOC,CO2,plantC,SMDE,DPME,RPME,BIOE,HUME,TOCE,CO2E,plantCE,R,zmax,pClay,dmax)

[plantRe,raine,PETe,Te,agee,LAIe,manure,cover]=read_wheather_data(scenario,species,ens);


inputs            = plant_input(EQinput,maxInput,agee,species,dmax);
plantCf1          = plantRe.*inputs;
%initialise variables
SMDf1             = zeros(size(PETe)); SMDf1(1)  = SMD(end);
DPMf1             = zeros(size(PETe)); DPMf1(1)  = DPM(end);
RPMf1             = zeros(size(PETe)); RPMf1(1)  = RPM(end);
BIOf1             = zeros(size(PETe)); BIOf1(1)  = BIO(end);
HUMf1             = zeros(size(PETe)); HUMf1(1)  = HUM(end);
%run model
[SMDf1,DPMf1,RPMf1,BIOf1,HUMf1,TOCf1,CO2f1]=carbon_turnover_cal(SMDf1,DPMf1,RPMf1,BIOf1,HUMf1,plantCf1,manure,raine,PETe,Te,cover,R,zmax,pClay);
%aggregate output
SMDo1             = [SMDE(end-5*365:end);SMD;SMDf1];
DPMo1             = [DPME(end-5*365:end);DPM;DPMf1];
RPMo1             = [RPME(end-5*365:end);RPM;RPMf1];
BIOo1             = [BIOE(end-5*365:end);BIO;BIOf1];
HUMo1             = [HUME(end-5*365:end);HUM;HUMf1];
TOCo1             = [TOCE(end-5*365:end);TOC;TOCf1];
CO2o1             = [CO2E(end-5*365:end);CO2;CO2f1];
plantCo1          = [plantCE(end-5*365:end);plantC;plantCf1];
%add IOM
IOMo1             = 0.049*TOCo1.^(1.139);
TOCo1             = TOCo1 + IOMo1;
%datetime vector
DT                = datetime(1992,1,1,0,0,0)+minutes(linspace(0,length(TOCo1)*24*60,length(TOCo1)))';
%yearly averages
win               = 360;
SMDo1             = windowed_average(SMDo1,win);
DPMo1             = windowed_average(DPMo1,win);
RPMo1             = windowed_average(RPMo1,win);
BIOo1             = windowed_average(BIOo1,win);
HUMo1             = windowed_average(HUMo1,win);
TOCo1             = windowed_average(TOCo1,win);
CO2o1             = windowed_average(CO2o1,win);
plantCo1          = windowed_average(plantCo1,win);
DT1               = windowed_average(DT,win);
%**************************************************************************
function [plantRe,raine,PETe,Te,agee,LAIe,manure,cover]=read_wheather_data(scenario,species,ens)

switch scenario
    case 'equilibrium'
        %equilibrium run for xxx number of years
        raine            = [];
        Te               = [];
        PETe             = [];
        LAIe             = [];
        plantRe          = [];
        years            = 1992:1997;
        for iy           = 1:length(years)
            load([pwd '\DATA\' 'climate_historical_daily_' num2str(years(iy)) '.mat'],'DT','rain','PET','T')
            LAI   = seasonal_sinewave_LAI_FFT(PET,species);
            %daily plant input distribution
            dLAI         = [-diff(LAI);0];
            dLAI(dLAI<0) = 0;
            plantR       = dLAI./sum(dLAI);
            plantRe      = [plantRe;plantR];
            LAIe         = [LAIe;LAI];
            raine        = [raine;rain];
            PETe         = [PETe;PET];
            Te           = [Te;T];
        end
        %soil cover
        cover            = LAIe>min(LAIe)*1.01;
        manure           = cover*0;
        agee             = [];
%--------------------------------------------------------------------------        
    case 'obervation'
        raine            = [];
        Te               = [];
        PETe             = [];
        LAIe             = [];
        plantRe          = [];
        agee             = [];
        years            = 1997;
        while years<2012
            years        = years+1;
            load([pwd '\DATA\' 'climate_historical_daily_' num2str(years) '.mat'],'DT','rain','PET','T')
            LAI          = seasonal_sinewave_LAI_FFT(PET,species);
            %daily plant input distribution
            dLAI         = [-diff(LAI);0];
            dLAI(dLAI<0) = 0;
            plantR       = dLAI./sum(dLAI);
            plantRe      = [plantRe;plantR];
            LAIe         = [LAIe;LAI];
            raine        = [raine;rain];
            PETe         = [PETe;PET];
            Te           = [Te;T];
            age          = 0*T+years-1998;
            agee         = [agee;age];
        end
        %soil cover
        cover            = LAIe>min(LAIe)*1.01;
        manure           = cover*0;
        
%--------------------------------------------------------------------------  
    case 'futureFpast'
        raine             = [];
        Te                = [];
        PETe              = [];
        LAIe              = [];
        plantRe           = [];
        agee              = [];
        years             = 2012;
        fyr               = 1991;
        while years<2078
            years         = years+1;
            if years<2023
                load([pwd '\DATA\' 'climate_historical_daily_' num2str(years) '.mat'],'DT','rain','PET','T')
            else
                fyr       = fyr+1;
                if fyr>2022
                    fyr   = 1992;
                end
                load([pwd '\DATA\' 'climate_historical_daily_' num2str(fyr) '.mat'],'DT','rain','PET','T')
            end
            LAI           = seasonal_sinewave_LAI_FFT(PET,species);
            %daily plant input distribution
            dLAI          = [-diff(LAI);0];
            dLAI(dLAI<0)  = 0;
            plantR        = dLAI./sum(dLAI);
            plantRe       = [plantRe;plantR];
            LAIe          = [LAIe;LAI];
            raine         = [raine;rain];
            PETe          = [PETe;PET];
            Te            = [Te;T];
            age           = 0*T+years-1998;
            agee          = [agee;age];
        end
        %soil cover
        cover            = LAIe>min(LAIe)*1.01;
        manure           = cover*0;
        
%--------------------------------------------------------------------------        
    case 'futureFchess'
        raine             = [];
        Te                = [];
        PETe              = [];
        LAIe              = [];
        plantRe           = [];
        agee              = [];
        years             = 2012;
        rcp               = '60';
        bias              = '_bias-corrected';
        while years<2079
            years         = years+1;
            if years<2023
                load([pwd '\DATA\' 'climate_historical_daily_' num2str(years) '.mat'],'DT','rain','PET','T')
            else
                load([pwd '\DATA\combo\' 'climate_future_daily_' num2str(years) '_rcp' rcp  '_ens' ens bias '.mat'],'DT','rain','PET','T')
                
            end
            LAI           = seasonal_sinewave_LAI_FFT(PET,species);
            %daily plant input distribution
            dLAI          = [-diff(LAI);0];
            dLAI(dLAI<0)  = 0;
            plantR        = dLAI./sum(dLAI);
            plantRe       = [plantRe;plantR];
            LAIe          = [LAIe;LAI];
            raine         = [raine;rain];
            PETe          = [PETe;PET];
            Te            = [Te;T];
            age           = 0*T+years-1998;
            agee          = [agee;age];
        end 
        %soil cover
        cover            = LAIe>min(LAIe)*1.01;
        manure           = cover*0;  
end


% sensorDepth              = 0.06;
% load([pwd '\DATA\' 'PARAMset' '_'  num2str(10001) '.mat'],'PARAMset')
% load([pwd '\DATA\' species '_' location '_' num2str(10001) '.mat'],'idx')
% load([pwd '\DATA\' species '_' location '_LOA.mat' ],'maxTSM')
% 
% Nb                       = length(idx);
% SMDmod                   = zeros(length(raine),Nb);
% tic
% parfor kk=1:Nb
%     %parameterset
%     param                = PARAMset(idx(kk),:);
%     OUT                  = EcoHydroPlot_ode(raine,PETe,species,sensorDepth,param,zmax,maxTSM);
%     SMDmod(:,kk)          = OUT(:,end);
% end
% toc
SMDmod=[];
%**************************************************************************
function b           = windowed_average(a,n)
% b                    = arrayfun(@(i) mean(a(i:i+n-1)),1:n:length(a)-n+1)';
b                    = arrayfun(@(i) median(a(i:i+n-1)),1:n:length(a)-n+1)';
%**************************************************************************
function hh          = patchPlot(X,Y0,Y1,col,alpha)
%shade area between two lines
hh=patch([X' fliplr(X')],[Y0' fliplr(Y1')],col,'Edgecolor','none','facea',alpha);
%set x-axis limit
% xlim([datetime(2022,01,01) X(end)])
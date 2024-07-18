function OUT      = RSEEP_ode(rain,PET,species,sensorDepth,param,zmax,TSM,Oold,Onew,age)
%--------------------------------------------------------------------------
%                       calibration parameters
%--------------------------------------------------------------------------
%decay exponent of porosity with depth, [-]
n                 = param(1);
%bulk-soil saturated hydraulic conductivity, [m/s]
k0                = param(2);
%exponent of decay soil water potential with saturation
p                 = param(3);
%--------------------------------------------------------------------------
%                       fixed parameters
%--------------------------------------------------------------------------
%porosity at observed depth
phid              = max(TSM);
%max porosity at the surface
phim              = phid^(1/n);
%chnage in max porosity due to OM change
Dphi              = 0.1224*(Onew-Oold)./(Onew+Oold);
%update
phi0              = phim + Dphi;
%conductivity modifier based on OM content change
M_K0              = 0.66*(phi0-phim)./(phi0+phim);
%--------------------------------------------------------------------------
%                            unit conversions
%--------------------------------------------------------------------------
%convert rain from mm/day to m/day
rain              = rain/1000;
%convert PET from mm/day to m/day
PET               = PET/1000;
%timestep [s]
dt                = 24*60*60;
%convert conductivity from m/s to m/day
k0                = k0*dt+M_K0;
%--------------------------------------------------------------------------
%                              pre-processing
%--------------------------------------------------------------------------
%leaf area index 
[LAI,LAIg,rE]     = seasonal_sinewave_LAI_FFT(PET,species,age);
%depths of equal storage control volumes
[SM,z,sm,phi,zB]  = equi_storage_depths(zmax,phi0,n);
%elevation
z0                = zmax-z;
%location of the soil moisture sensors
ids1              = find(abs(z-sensorDepth)   == min(abs(z-sensorDepth)),1);
ids2              = find(abs(z-0.2)   == min(abs(z-0.2)),1);
ids3              = find(abs(z-0.4)   == min(abs(z-0.4)),1);
ids4              = find(abs(z-0.6)   == min(abs(z-0.6)),1);
ids5              = find(abs(z-1)     == min(abs(z-1)),1);
%surface cover fraction (SCF)
SCF               = 1-exp(-rE*LAI);
%--------------------------------------------------------------------------
%                                 initialise
%--------------------------------------------------------------------------
%number of timesteps
N                 = length(rain);
%canopy storage
Sc                = zeros(N,1);
%canopy evap
Ec                = zeros(N,1);
%bulk soil soil evap
Es                = zeros(N,1);
%bulk soil transpiration
Tr                = zeros(N,1);
%bulk soil vertical flux
Qv                = zeros(N,1);
%bulk soil lateral flux
Qb                = zeros(N,1);
%bulk lateral/preferential/overland flow
LPOF              = zeros(N,1);
%layer-wise storage
s                 = zeros(N,length(z));
%layer-wise actual soil evaporation
ase               = zeros(N,length(z));
%layer-wise actual transpiration
astr              = zeros(N,length(z));
%bulk soil storage
S                 = zeros(N,1);
%bulk soil storage at time=0
S(1)              = SM;
%canopy storage at time=0
Sc(1)             = rain(1);
%eliminate first rain for mass balance
rain(1)           = 0;
%--------------------------------------------------------------------------
%                               march in time
%--------------------------------------------------------------------------
%machine precision
e                 = 1e-16;
%time loop
for ii            = 2:N
    %update bulk storage
    S(ii)         = S(ii-1);
    %update canopy storage
    Sc(ii)        = Sc(ii-1) + rain(ii);
    %----------------------------------------------------------------------
    %potential transpiration
    Mod           = 1 + (LAI(ii)-LAIg(ii))./(LAI(ii)+LAIg(ii));
    %adjust potential transpiration of the 'species' relative to grass
    PET0          = PET(ii)*Mod;
    %----------------------------------------------------------------------
    %max interception storage (converted from mm to m)
    Smax          = 0.2*LAI(ii)/1000;
    %find max canopy storage
    ScM           = min(Sc(ii),Smax);
    %throughfall
    tf            = Sc(ii) - ScM;
    %update canopy
    Sc(ii)        = Sc(ii) - tf;
    %update bulk
    S(ii)         = S(ii)  + tf;
    %canopy evap
    Ec(ii)        = min(PET0,Sc(ii));
    %update canopy storage
    Sc(ii)        = Sc(ii) - Ec(ii);
    %update remaining PET
    PET0          = PET0-Ec(ii);
    %potential transpiration
    Tr0           = SCF(ii)*PET0;
    %potential evaporation
    Ep0           = (1-SCF(ii))*PET0;
    %----------------------------------------------------------------------
    %storage cannot exceed maximum in the bulk soil
    %if it does, remove as lateral/preferential/overland flow (LPOF)
    if S(ii)>SM
        LPOF(ii)  = S(ii)-SM;
        %reset to max
        S(ii)     = SM;
    end
    %----------------------------------------------------------------------
    %bulk unsaturate hydraulic conductivity
    Ku            = k0*(S(ii)./SM).^p;
    %vertical flux across the soil profile
    qv            = Ku.*S(ii)./zmax;
    %it cannot exceed totalstorage
    Qv(ii)        = min(qv,S(ii));
    %update bulk storage
    S(ii)         = S(ii) - Qv(ii);
    %----------------------------------------------------------------------
    %reconstruct profile soil storage from bulk storage
    [s(ii,:),SAT] = reconstruct(S(ii),SM,z,sm,zmax);
    %saturation of each layer
    sat           = s(ii,:)'./(sm+e);
    %weight for transpiration
    if strcmp(species,'grassland') || strcmp(species,'grass')
        %short-rooted plants can't drink from deeper depth
        w2        = z0./zmax ;
    else
        %longer-rooted plants can alter their water source according to saturation
        w2        = SAT.*(z0./zmax) + (1-SAT).*sat;
    end
    %weights should add up to one
    w2            = w2./(sum(w2)+e);
    %layer-wise potential transpiration
    pstr          = w2.*Tr0;
    %actual transpiration
    astr(ii,:)    = min(pstr.*sat.^p,s(ii,:)');
    %total transpiration
    Tr(ii)        = min(sum(astr(ii,:)),S(ii));
    %update bulk storage
    S(ii)         = S(ii) - Tr(ii);
    %----------------------------------------------------------------------
    %reconstruc profile soil storage from bulk storage
    s(ii,:)       = reconstruct(S(ii),SM,z,sm,zmax);
    %saturation of each layer
    sat           = s(ii,:)'./(sm+e);
    %weights should add up to one
    pse=integrate_soil_evap(zB,Ep0);
    
    %actual soil evap for each layer
    ase(ii,:)     = min(pse.*sat.^p,s(ii,:)');
    %total soil evap
    Es(ii)        = min(sum(ase(ii,:)),S(ii));
    %update bulk storage
    S(ii)         = S(ii) - Es(ii);
end
SMD               = SM-S;
%--------------------------------------------------------------------------
%calculate mass error
err               = mass_balance(rain,S,Sc,Ec,Es,Tr,LPOF,Qv,Qb);
err               = err/1000 + rain*0;
%soil moisture content at target layers (observed layers)
VSM1              = s(:,ids1)./sm(ids1)*phi(ids1)./1000;
VSM2              = s(:,ids2)./sm(ids2)*phi(ids2)./1000;
VSM3              = s(:,ids3)./sm(ids3)*phi(ids3)./1000;
VSM4              = s(:,ids4)./sm(ids4)*phi(ids4)./1000;
VSM5              = s(:,ids5)./sm(ids5)*phi(ids5)./1000;
%replace the first value
VSM1(1)           = VSM1(2);
LPOF(1)           = LPOF(2);
Ec(1)             = Ec(2);   
Es(1)             = Es(2);
Tr(1)             = Tr(2);
Qb(1)             = Qb(2);
Qv(1)             = Qv(2);
Sc(1)             = Sc(2);
S(1)              = S(2);
err(1)            = err(2);
VSM2(1)           = VSM2(2);
VSM3(1)           = VSM3(2);
VSM4(1)           = VSM4(2);
VSM5(1)           = VSM5(2);
SMD(1)            = SMD(2);
%(convert to mm)      1     2     3    4     5    6    7     8    9   10     11   12   13   14   15             16   
OUT               = [VSM1  LPOF   Ec   Es    Tr   Qb   Qv    Sc   S   err   VSM2 VSM3 VSM4 VSM5 zmax*ones(N,1) SMD]*1000;
%eliminate spinup period
OUT               = real(OUT);
OUT(isnan(OUT))   = 0;
%**************************************************************************
function Ep = integrate_soil_evap(zB,Epmax)
%integrates f=2./(1+exp(alpha*z))
z1          = zB(1:end-1);
z2          = zB(2:end);
dz          = z2-z1;
alpha       = 500;
I1          = 2*Epmax*(z1-1/alpha*log(exp(alpha*z1)+1));
I2          = 2*Epmax*(z2-1/alpha*log(exp(alpha*z2)+1));
Ep          = (I2-I1)./dz;
%**************************************************************************
function err      = mass_balance(rain,S,Sc,Ec,Es,Tr,OLF,Qv,Qb)
%change in total storage rom t=0 to t=end
dS                = S(end)-S(1) + Sc(end)-Sc(1);
%total influx to the system
Qin               = sum(rain);
%total outflux
Qout              = sum(Ec+Es+Tr+OLF+Qv +Qb);
%percentage error relative initial storage
err               = abs(dS+Qout-Qin)./(S(1)+Sc(1))*100;
%**************************************************************************
function [s,SAT]  = reconstruct(S,SM,z,sm,zm)
e                 = eps;
%bulk saturation
SAT               = S./(SM+eps);
%mix them
W                 = SAT.*sm./(max(sm)+e) + (1-SAT).*z./(zm+e);
%weights should add up to one for mass conservation
W                 = W./(sum(W)+e);
%storage of each layer
s                 = W.*S;
%**************************************************************************
function [SM,z,sm,phi,zB] = equi_storage_depths(zm,phi0,n)
%--------------------------------------------
%calculates equi-storage layers
%--------------------------------------------
%thickness of the first layer
z1                     = 0.02;
%average porosity of the first layer
phi1                   = integrate_fun(n,phi0,0,z1);
%max storage of the first layer
sm1                    = z1*phi1;
%average porosity of the entire soil profile
PHI                    = integrate_fun(n,phi0,0,zm);
%max storage in the soil
SM                     = zm*PHI;
%total number of layers
nL                     = ceil(SM/sm1);
sB                     = repmat(SM/nL,nL,1);
sB                     = cumsum(sB);
RHS                    = (log(sB*(1-n)./n./phi0+1) + (1-n)*log(n))/(1-n+eps);
zB                     = [0;exp(RHS)-n];
%calculate mid-point depth
z                      = ((zB(1:end-1)+zB(2:end))/2);
%average porosity of all layers
phi                    = integrate_fun(n,phi0,zB(1:end-1),zB(2:end));
sm                     = phi.*diff(zB);
SM                     = sum(sm);
%**************************************************************************
function phi     = integrate_fun(n,phi0,z1,z2)
%--------------------------------------------
%integrates phi=phi0./(1+zm/n)^n
%--------------------------------------------
%machine precision
e                = 1e-16;
%range
dz               = z2-z1;
%lower limt of the integral
X1               = -(n+z1).*((n+z1)./n).^(-n)./(n-1)*phi0;
%upper limit
X2               = -(n+z2).*((n+z2)./n).^(-n)./(n-1)*phi0;
%sum of porosity
phi              = (X2-X1);
%average porosity
phi              = phi./(dz+e);
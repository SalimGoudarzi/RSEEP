%**************************************************************************
function [LAI,LAIg,rE,R]     = seasonal_sinewave_LAI_FFT(PET,species,age)
%diamter at breast height (DBH) at current age
[DBH]                        = biomass_fun(age);
%DBH at maturation (after 200 years)
[DBHmax]                     = biomass_fun(200);
%fractions of the two
frac                         = (DBH./DBHmax)+0.01;
%--------------------------------------------------------------------------
%get the shape of the PET wave for seasonality
WAVE                         = Signal_in_noise_processing(PET);
%leaf area index of species
LAI                          = sine_wave(species,WAVE);
%leaf area index of grass (reference crop)
LAIg                         = sine_wave('grass',WAVE);

if strcmp(species,'grass') || strcmp(species,'grassland')
    frac=0*age;
end
LAI                          = frac.*(LAI) +(1-frac).*LAIg;
rE                           = frac.*0.7+(1-frac).*0.3;
R                            = frac.*0.25+(1-frac).*0.67;
%**************************************************************************
function [LAI,minLAI,maxLAI] = sine_wave(species,WAVE)
%minimum and maximum LAI
switch species
    %----------------------------------------------------------------------
    case 'conifer'
        maxLAI               = 7;
        minLAI               = 6;
    %----------------------------------------------------------------------
    case 'beech'
        maxLAI               = 9;
        minLAI               = 1.5;
    %----------------------------------------------------------------------
    case 'grassland'
        maxLAI               = 3;
        minLAI               = 1.5;
    %----------------------------------------------------------------------
    case 'grass'
        maxLAI               = 2;
        minLAI               = 1;
    %----------------------------------------------------------------------
end
WAVE                         = WAVE.^3;
WAVE                         = WAVE-min(WAVE);
WAVE                         = WAVE./max(WAVE);
LAI                          = (1.2*maxLAI-minLAI).*WAVE+ minLAI;
LAI(LAI>maxLAI)              = maxLAI;
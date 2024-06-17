function result=separateEffectsPE(simu,source,heights,varargin)
%result=separateEffectsPE(simu,source,varargin)
%A function that apply the Parabolic Equation model with different
%phenomenenon(t and w proifiel, turbulence...) and substract them to return
%teh effect of each ones.
% Inputs:
%   simu: a structure containing the parameter of the simulation:
%    Required fields:
%       simu.Zmax:       the maximum altitude of the simulation(m)
%
%       simu.Xmax:       the maximum distance of the simulation in the x(or r)
%                         direction(m)
%
%       simu.Rstep:       the length of a step in  the x(or r) direction(m)
%
%
%       simu.Zstep:        The height the vertical layers(m)
%
%       simu.Zabs:        The height at which the absorbing layer start(m)
%
%       simu.Profile:       The name of a csv file containing
%                             a complete profile to use
%
%       simu.ZGround:        The normalized impedance of the ground
%
%    source: a structure containing the parameter of the sound source:
%    Required fields:
%       source.freq:        The frequency of the sound emmited by the
%                           source(Hz)
%       source.Zsource:        The height of the sound source(m)
%
%       source.Amplitude:      The amplitude at 10cm of the source(dB)
%
%
%
%   heights: the list of heights at wich the efffects must be computed(m)
%
%Optional inputs:
%    MaxAlt             : teh maximulm altitude for teh averaging (m)
%
%    Turbulence       : default=1, if equal to 0 no turbulent simulation is
%    computed.
%
%   Outputs:
%       result: a structure of length n, the number of heights required
%
%       Fields:
%
%           result(k).h: height for this line(m)
%
%           result(k).r : the range vector(m)
%
%           result(k).complete: the vector of losses for teh complete
%               simulation(dB)
%
%           result(k).turbulence: the vector of losses caused by turbulence(dB)
%
%           result(k).ground: the vector of losses caused by the ground (dB)
%
%           result(k).temperature_pr: the vector of losses caused by the temperature profile (dB)
%
%           result(k).wind_pr: the vector of losses caused by the wind profile (dB)
%
%           result(k).wind_pr_opposite: the vector of losses caused by the wind profile with teh wind direction shifted by 180 deg (dB)
%
%           result(k).geometric: the vector of losses caused by geometric spreading (dB)
%
%           result(k).classical: the vector of losses caused by classical attenuation (dB)
%
%           result(k).non_classical: the vector of losses caused by non classical attenuation(dB)
%
%

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"simu");
addRequired(p,"source");
addRequired(p,"heights");

default_maxalt=20;
addParameter(p,'MaxAlt',default_maxalt);%max altitude for averaging$
default_turb=1;
addParameter(p,'Turbulence',default_turb);%should turbulence losses be computed?
parse(p,simu,source,heights,varargin{:})

%%%%%%%%
maxalt=p.Results.MaxAlt;
is_turb=p.Results.Turbulence;

pe_noground=ParabolicEquation(simu,source,'NoGround',1,'Averaged',1,'AverageLim',maxalt);
pe_nopr=ParabolicEquation(simu,source,'Averaged',1,'AverageLim',maxalt);
simunoWind=simu;
simunoWind.WSpeed=0;
pe_onlyT=ParabolicEquation(simunoWind,source);

%Only wind profile of simu taken into account
pe_onlyW=ParabolicEquation(simu,source,'Averaged',2,'AverageLim',maxalt);

%Wind direction opposite of the simu input
simuInvertWind=simu;
simuInvertWind.WDirection=simu.WDirection+180;
pe_onlyInvertW=ParabolicEquation(simuInvertWind,source,'Averaged',2);
pe_all_refraOppositeWind=ParabolicEquation(simuInvertWind,source);

pe_all_refra=ParabolicEquation(simu,source);
if is_turb
    pe_complete=ParabolicEquation(simu,source,'Turbulence',simu.turbulence);
else
    pe_complete=pe_all_refra;
end

%% Substracting to separate effects
zneg=pe_noground.Z<0;
gSPL=flipud(pe_noground.SPL);
gSPL(zneg,:)=[];
gSPL=flipud(gSPL);
ground=pe_nopr.SPL-gSPL;
temperature_pr=pe_onlyT.SPL-pe_nopr.SPL;
wind_pr0=pe_onlyW.SPL-pe_nopr.SPL;
wind_pr=pe_all_refra.SPL-pe_onlyT.SPL;

wind_pr_opposite=pe_all_refraOppositeWind.SPL-pe_onlyT.SPL;
wind_pr_opposite0=pe_onlyInvertW.SPL-pe_nopr.SPL;

turbulence=pe_complete.SPL-pe_all_refra.SPL;
complete=pe_complete.SPL;

%%
n=length(heights);
for k=1:n
    [~,y]=min(abs(heights(k)-pe_all_refra.Z));%index of 1.5m altitude
    y=length(pe_all_refra.Z)-y;%flip upside down
    result(k).h=heights(k);
    h=heights(k);
    r=pe_all_refra.R;
    result(k).r=r;
    result(k).complete=source.Amplitude-complete(y,:);
    result(k).turbulence=turbulence(y,:);
    result(k).ground=ground(y,:);
    result(k).temperature_pr=temperature_pr(y,:);
    result(k).wind_pr=wind_pr(y,:);
    result(k).wind_pr0=wind_pr0(y,:);
    result(k).wind_pr_opposite=wind_pr_opposite(y,:);
    result(k).wind_pr_opposite0=wind_pr_opposite0(y,:);

    d=sqrt(r.^2+(source.Zsource-heights(k))^2);
    result(k).d=d;
    result(k).geometric=20*log10(0.1./d);
    atm=makeMartianAtm('Simu',simu,'Altitude',0.5*(source.Zsource+h));
    att=attenuationModel(atm,'FreqArray',[source.freq]);
    dbM_c=-20*att.alpha_c/log(10);
    dbM_r=-20*att.alpha_r/log(10);
    result(k).classical=dbM_c*d;
    result(k).non_classical=dbM_r*d;

end

function [grid ,profiles]=make_kgridPE(simu,source)
%grid=make_grid(simu,source)
%A function that return the grid of wavenumber(k) used in the Parabolic Equation model.
% Inputs:
%   simu: a structure containing the parameter of the simulation:
%    Required parameters:
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
%       simu.WSpeed:       If this field is not equal to 0 the wind profiel
%                           will be added to the speed of sound
%
%       simu.WDirection:       The direction (in deg) of the wind. 0 deg is
%                               the direction of the simulation, 180 the opposite direction.
%
%
%   Optionnal parameter:
%       simu.NoGround: equal to 1 if the simulation must be done without
%       ground
%    source: a structure containing the parameter of the sound source:
%    Required parameters:
%       source.freq:        The frequency of the sound emmited by the
%                           source(Hz)
%
%


zvect=0:simu.Zstep:simu.Zmax;
m=length(zvect);
rvect=0:simu.Rstep:simu.Xmax;

omega=source.freq*2*pi;

%Determining which attenuation and speed of sound should be chosen based on the source frequency.
Profile=readtable(simu.Profile);
power10=round(log10(source.freq)); % order of magnitude of the sound source frequency (needed to choose the adequate sos )
lowpower10=floor(log10(source.freq));%first power of ten smaller than the required frequency
uppower10=floor(log10(source.freq)+1);%first power of ten greater than the required frequency
switch lowpower10
    case 0 %source frequency between 1 and 10 Hz
        Cprofile=Profile.C_1Hz+(log10(source.freq)-0)*(Profile.C_10Hz-Profile.C_1Hz)/1;%1 being the log distance between 1 and 10 Hz
        logAttProfile=log10(Profile.Att_1Hz)+(log10(source.freq)-0)*(log10(Profile.Att_10Hz)-log10(Profile.Att_1Hz))/1;%same thing but the affine approximation is in teh log-log diagram
        Attprofile=10.^logAttProfile;
    case 1 %source frequency between 10 and 100 Hz
        Cprofile=Profile.C_10Hz+(log10(source.freq)-lowpower10)*(Profile.C_100Hz-Profile.C_10Hz)/1;%1 being the log distance between 1 and 10 Hz
        logAttProfile=log10(Profile.Att_10Hz)+(log10(source.freq)-0)*(log10(Profile.Att_100Hz)-log10(Profile.Att_10Hz))/1;%same thing but the affine approximation is in teh log-log diagram
        Attprofile=10.^logAttProfile;
    case 2 %source frequency between 100 and 1 000 Hz
        Cprofile=Profile.C_100Hz+(log10(source.freq)-lowpower10)*(Profile.C_1000Hz-Profile.C_100Hz)/1;%1 being the log distance between 1 and 10 Hz
        logAttProfile=log10(Profile.Att_100Hz)+(log10(source.freq)-lowpower10)*(log10(Profile.Att_1000Hz)-log10(Profile.Att_100Hz))/1;%same thing but the affine approximation is in teh log-log diagram
        Attprofile=10.^logAttProfile;
    case 3 %source frequency between 1 000 and 10 000 Hz
        Cprofile=Profile.C_1000Hz+(log10(source.freq)-lowpower10)*(Profile.C_10000Hz-Profile.C_1000Hz)/1;%1 being the log distance between 1 and 10 Hz
        logAttProfile=log10(Profile.Att_1000Hz)+(log10(source.freq)-lowpower10)*(log10(Profile.Att_10000Hz)-log10(Profile.Att_1000Hz))/1;%same thing but the affine approximation is in teh log-log diagram
        Attprofile=10.^logAttProfile;
    case 4 %source frequency greater than 10 000 Hz
        Cprofile=Profile.C_10000Hz;
        Attprofile=Profile.Att_10000Hz;
end

%Adding wind

if simu.WSpeed~=0
    Cwind=cos(simu.WDirection*pi/180)*Profile.wind;
else
    Cwind=Cprofile*0;
end
%Computing every parameter for the layers.
Attlayer=interp1(Profile.Altitude,Attprofile,zvect);
Tlayer=interp1(Profile.Altitude,Profile.t,zvect);
Cthermlayer=interp1(Profile.Altitude,Cprofile,zvect);
Windlayer=interp1(Profile.Altitude,Cwind,zvect);
%scaling wind layer
if simu.WSpeed>1%WSpeed is the magnitude of the wind(>0)
    [~,id15]=min(abs(1.5-zvect));
    Windlayer=Windlayer*simu.WSpeed/abs(Windlayer(id15));
end

%Averaging if needed
id=find(zvect<=simu.AverageLim);
id=id(end);
if simu.Averaged==1%everything averaged
    Attlayer(:)=mean(Attlayer(1:id));
    Tlayer(:)=mean(Tlayer(1:id));
    Windlayer(:)=mean(Windlayer(1:id));
    Cthermlayer(:)=mean(Cthermlayer(1:id));
elseif simu.Averaged==2%only t averaged (and consequently Ctherm, and att)
    Attlayer(:)=mean(Attlayer(1:id));
    Tlayer(:)=mean(Tlayer(1:id));
    Cthermlayer(:)=mean(Cthermlayer(1:id));
end
%Computing Clayer, the layer of effective speed of sound used to compute
%Kgrid
Clayer=Cthermlayer+Windlayer;

%Storing every profile for future display
profiles.Ceff=Clayer;
profiles.Ctherm=Cthermlayer;
profiles.T=Tlayer;
profiles.Wind=Windlayer;

%Compute the attenuation coeff for the artificial absorbing layer added
%between Zabs and Zmax
At=2;%to be tuned
if source.freq>200
    At=source.freq/1000*10;
end
abs_coeff=1i*At*(zvect-simu.Zabs).^2/(simu.Zmax-simu.Zabs)^2;
absMat=repmat(flipud(abs_coeff.'),1,length(rvect));
%Non turbulent case
Cmat=repmat(flipud(Clayer'),1,length(rvect));%flipud because in matlab 1 is the indexe of the top row whereas in the convention we use here it 1 is the indexe of the first(ie bottopm) layer
Attmat=repmat(flipud(Attlayer'),1,length(rvect));
%

grid=omega*Cmat.^(-1)+1i*Attmat;%k=omega/c+i*alpha

[~,ind_abs]=min(abs(zvect-simu.Zabs));%index of Zabs

%adding the upper bound attenuation layer
grid(1:m-ind_abs,:)=grid(1:m-ind_abs,:)+absMat(1:m-ind_abs,:);

if isfield(simu,'NoGround')
    if simu.NoGround
    %Make a lower bound attenuation layer
    neg_zvect=-(simu.Zmax-simu.Zabs):simu.Zstep:0;
    abs_coeff_low=1i*At*(neg_zvect).^2/(simu.Zmax-simu.Zabs)^2;
    grid_low=repmat(grid(m,:),length(neg_zvect),1);
    absMat_low=grid_low+repmat(flipud(abs_coeff_low.'),1,length(rvect));
    
    grid=[ grid;absMat_low(2:end,:)];
    end
end



end

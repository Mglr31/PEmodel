function result=prousGroundAndi(simu,source,varargin)
%result=prousGroundAndi(simu,source,varargin)

%
% Computes the superposition of the direct and reflected fields above porous
% ground. 
%
% The model assumes LOCAL REACTION i.e. semi-infinite homogeneous porous medium.
% > for (semi-)rigid backed porous layer, extended reaction must be used.
%
% Sources:
% Attenborough, "Predicting Outdoor Sound" (book)
% Champoux-Allard 1991 paper (for porous properties)
%Thsi function use the same format of output and input than
%ParabolicEquation to facilitate comparison
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
%       simu.ZGround:        The impedance of the ground
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
%Optional inputs:
%    mic: a structure containing the positions of the sound receptors:
%    Required fields:
%       mic.R: an array containing the positions of the mics along the r axis(m)
%
%       mic.Z:an array containing the positions of the mics along the z axis(m)
%
%   'Name',name: the name of the simulation(string)
%
%   'Turbulence',turbulence: a structure representing the kind of
%                            atmospheric turbulence that will be added to
%                            the simulation. Default value:[]: no
%                            turbulence.
%
%               Required fields:
%
%                   turbulence.wavenumbers: an array of the wavenumbers used to compute the turbulent field
%
%                   turbulence.type: the type of turbulence. Possible values:"gaussian","karman","experimental"
%                                      If the type is gaussian or karman,
%                                      turbulence needs to have two
%                                      additional fields:mu0 and a.
%                                      If the type is experimental,
%                                      turbulence need to have the field
%                                      F,F being the experimental spectral
%                                      density of mu.
%
%                  turbulence.NofRealisation: the number of random
%                                           realisation of the turbument fields that will be
%                                           computed for this simulation. The result is the
%                                           logarithmic mean of these realisation.
%
% Outputs:
%   result: A structure containing all the results of a PE simulation(it is
%           the output of the function ParabolicEquation).
%   Required fields:
%       result.SPL: the grid of sound pressure levels(db?)
%
%       result.P: the grid of complex amplitude of the pressure(Pa?)
%
%       result.K : the grid of wavenumber
%
%       result.Z: the vector of z coordinates(altitude,m)
%
%       result.R: the vector of r coordinates(range,m)
%
%       result.Freq_Source: the frequency of the source that was simulated(Hz)
%
%       result.Ztop: the maximum height of the physical
%                   simulation(excluding the absorbing layer, in m)

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"simu");
addRequired(p,"source");
default_mic=[];
addParameter(p,'mic',default_mic);
default_name="Andi porous ground dual ray simulation";
validname= @(x) isstring(x);
addParameter(p,'Name',default_name,validname);
default_turb=[];
addParameter(p,'Turbulence',default_turb);
parse(p,simu,source,varargin{:})
%%%%%%%%
result.Name=p.Results.Name;
% range and height of receiver:
xx = 0:simu.Rstep:simu.Xmax;
zz = 0:simu.Zstep:simu.Zmax;
m=length(zz);
[x,z] = meshgrid(xx,zz);
zs=source.Zsource;

% direct and reflected Tx-Rx distances:
r1 = sqrt(x.^2+(z-zs).^2);
r2 = sqrt(x.^2+(z+zs).^2);

% incidence angle for specular reflection:
theta = acos((zs+z)./r2);

%> find iheight and irange for detector height and range
% dx = xmax/(Nx-1);
% dz = zmax/(Nz-1);
% 
% irange = floor(xdet/dx+1); xdet = xx(irange)
% iheight = floor(zdet/dz+1); zdet = zz(iheight)


% Mars atm:
pr=readtable_header(simu.Profile);
T0 = mean(pr.t);
p0 = mean(pr.p);
rho0 = mean(pr.rho);

load("modelMars.mat");
model=modelMars;

    model.T=T0;
    model.P=p0;
    molarFrac=[mean(pr.n2) mean(pr.co2) mean(pr.ar) mean(pr.o2)];
    model.MolarFrac=molarFrac;
[gas_t,vib_t]=makeModelTables(model);%create and save tables corresponding to the model molecule and modes parameters
model.gasT=gas_t;
model.vibT=vib_t;
model=addDerivativePArameter(model);

mu0 = model.etha; 
kappa0 = model.kappa;
cp0 = model.cp;  % 33.35 J/mol K converted to J/kgK
cv0 = model.cv;  % 25.04 J/mol K converted to J/kgK
c0 = mean(pr.C_100Hz);

% Mars porous specs (cavalier values...)
% poro = 0.365;
% tort = 1.87;
% res = 21250;    % viscous flow resistivity [Pa s/m2]
%> from Sizemore et al, "Laboratory measurements of tortuosity and
% permeability in Mars analog soils":
poro = 0.365*(1+0.0);
%poro=0.6;
tort = 1.87*(1+0.0);
perm = 1.15e-12*(1+0.0);    % permeability = mu/resistivity [m2]
res = mu0/perm;    % viscous flow resistivity [Pa s/m2]



%% Complex density and compressibility of porous surface

f0=source.freq;
w0 = 2*pi*f0; 
k0 = w0/c0;
a_ext0=mean(pr.Att_100Hz);
rho1 = complex_density_JCAL(f0,poro,tort,res,p0,rho0,mu0,kappa0,cp0,cv0);
com1 = complex_compressibility_JCAL(f0,poro,tort,res,p0,rho0,mu0,kappa0,cp0,cv0);

k1 = w0*sqrt(rho1.*com1);
c1 = w0./real(k1);
beta = rho0*c0/(rho1*c1);  
beta=1/simu.ZGround;

%% Complex-valued reflection coefficient of porous surface 

% Plane-wave reflection coefficient:
R0 = (cos(theta)-beta)./(cos(theta)+beta);
Z=simu.ZGround;
R0=(Z-1)/(Z+1);
% "w0" and "F(w0)" in Attenborough:
aw = (1+1i)/2*sqrt(k0*r2).*(cos(theta)+beta);  
    
    Fw = 1+1i*sqrt(pi).*aw.*erfcx_2(-1i*aw,0.5); 
    %> 2nd arg is the accuracy criterion (see erfcx_2 function in this folder)

% Modified reflection coefficient: 
%> (Weyl-Van der Pol approximation. Replace it by solving the Weyl 
%  integral over dq numerically)
% R0 = 0;
% Fw = 0;
Q0 = R0 +(1-R0).*Fw; % second term is the spherical-wave correction
% Q0 = 0;
Q0=R0;
displayExpNotation(Q0);

%% Received field

        p_direct = exp(1i*k0*r1)./(4*pi*r1);   
        p_echo = Q0.*exp(1i*k0*r2)./(4*pi*r2);
        % > add attenuation factors:
        p_direct = p_direct.*exp(-a_ext0*r1);
        p_echo = p_echo.*exp(-a_ext0*r2);
        p = p_direct + p_echo;



% intensity (dB):
pref = 20e-6;   % 20 uPa
I = 20*log10(abs(p/pref));

result.SPL=I;
[~,indr]=min(abs(xx-0.1));
[~,indz]=min(abs(zz-source.Zsource));
%indz=m-indz;%the first row of the grid correspond to the maximum altitude
spl0=result.SPL(indz,indr);
%Scaling the entire SPL grid
result.SPL=flipud(result.SPL-(spl0-source.Amplitude));

result.Z=zz;
result.R=xx;
result.Ztop=simu.Zabs;

end

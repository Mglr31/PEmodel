function init=make_initPE(simu,source,ka)
%init=q0PE(simu,ka,zs)
%A function that return a vector that is the starting
%field in the PE method
% Inputs:
%   simu: a structure containing the parameter of the simulation:
%    Required parameters:
%       simu.Zmax:       the maximum altitude of the simulation(m)
%
%       simu.Zstep:        The height the vertical layers(m)
%                            a complete profile to use
%
%       simu.ZGround:        The impedance of the ground
%
%    source: a structure containing the parameter of the sound source:
%    Required parameters:
%       source.Zsource:        The height of the sound source(m)
%                           
%   ka: the wave number of the simulation (non refractive atmosphere case)
%
A0=1.3717;
A2=-0.3701;
B=3;
zs=source.Zsource;
Z=simu.ZGround;
z=simu.zvect;
%Direct starting field
qdirect=sqrt(1i*ka)*(A0+A2*ka^2*(z-zs).^2).*exp(-ka^2*(z-zs).^2/B);
%Reflected starting field
C=(Z-1)/(Z+1);%reflection coefficient(plane wave)
qreflect=C*sqrt(1i*ka)*(A0+A2*ka^2*(z-zs).^2).*exp(-ka^2*(z+zs).^2/B);
if ~simu.NoGround
    init=flipud(qdirect+qreflect).';
else
    init=flipud(qdirect).';
end
a=1;
%init=flipud(qdirect.');
end
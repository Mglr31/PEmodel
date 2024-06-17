function result=get_audible_range(simu,source,threshold,varargin)
%result=get_audible_range(simu,source,varargin)
%A function that return the maximum range at which the sound source "source" is
%audible. The source is audible if the received amplitude is greater than
%threshold.
%
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
%
%       'Averaged':     default 0, if equal to 1 the atmospheric profile will
%                       be averaged to suppress the effect of sound
%                       refraction,if 2 only the T profile will be averaged
%
%       "Light":    default 0, if equal to 1 mu fields wont be saved and
%                   spl will be saved in songle format so the result is lighter
%
%      "Resolution":    default 0, if different than zero the resolution of
%                the output SPL grid is diminished "Resolution" times. Ex:
%                if one wants to compute a millimetric grid but store a
%                centimetric grid: "Resolution",10
%
%       "NoGround": by default equal to 0, equal to 1 if one want a
%       simulation without ground
%
%   threshold : the audibility threshold in dB
%
% Outputs:
%   result: a structure containig the maximum range of audibility,i.e the
%   range for which the received SPL is equal to threshold
%   Required fields:
%       result.r: the range vector: from 0 to Xmax
%
%       result.z: the altitude vector from 0 to Zmax
%
%       result.xl,yl, xr,yr : the (x,y) coordinate of the range limit for
%       the left and right simulation
%
%       result.PEright
%
%       result.PEleft
%
%       result.Zsource: Height of the source(source.Zsource) in m
%

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"simu");
addRequired(p,"source");
addRequired(p,"threshold");
default_mic=[];
addParameter(p,'mic',default_mic);
default_name="Audibility area";
validname= @(x) isstring(x);
addParameter(p,'Name',default_name,validname);
default_turb=[];
addParameter(p,'Turbulence',default_turb);
default_av=0;
addParameter(p,'Averaged',default_av);
default_light=0;
addParameter(p,'Light',default_light);
default_res=0;
addParameter(p,'Resolution',default_res);
default_ng=0;
addParameter(p,'NoGround',default_ng);
parse(p,simu,source,threshold,varargin{:})


%PE simulation
%Right
%%
PEright=ParabolicEquation(simu,source,'Name',"PE right");

spl=flipud(PEright.SPL);
is_audible=spl>threshold;
shape=spl;
shape(~is_audible)=0;
shape(is_audible)=1;
near_ground=shape(1:10,1:10);%area close to r=0 and z=0
aud_near_ground=sum(near_ground(:))/100;


ed=edge(shape,'Canny');
sz=size(ed);
k=find(ed);
[row,col] = ind2sub(sz,k);

if aud_near_ground>0.5%if there is sound near the ground
    d02=(sz(2)-col).^2+1e6*(0-row).^2;%weighted square of the distance to the starting point(lower right corner)
else
     d02=1e6*(0-col).^2+(0-row).^2;%weighted square of the distance to the starting point(lower of the left points)
end
     [~,closest]=min(d02);
x0=col(closest);y0=row(closest);%closest point

if aud_near_ground>0.5%if there is sound near the ground
    dend2=(0-col).^2+(0-row).^2;%square of the distance to the ending point(lower left corner)
else
    dend2=1e6*(0-col).^2+(sz(1)-row).^2;%square of the distance to the ending point(higher of the leftmost points)
end
[~,closest]=min(dend2);
xend=col(closest);yend=row(closest);%closest point

[xr,yr]=draw_line(col,row,x0,y0,xend,yend);
xr=xr*simu.Rstep;yr=yr*simu.Zstep;
%%
%Left
simu.WDirection=(simu.WDirection+180);
PEleft=ParabolicEquation(simu,source,'Name',"PE left");
spl=flipud(PEleft.SPL);
is_audible=spl>threshold;
shape=spl;
shape(~is_audible)=0;
shape(is_audible)=1;
near_ground=shape(1:10,1:10);%area close to r=0 and z=0
aud_near_ground=sum(near_ground(:))/100;

ed=edge(shape,'Canny');
sz=size(ed);
k=find(ed);
[row,col] = ind2sub(sz,k);

if aud_near_ground>0.5%if there is sound near the ground
    d02=(sz(2)-col).^2+1e6*(0-row).^2;%weighted square of the distance to the starting point(lower right corner)
else
     d02=1e6*(0-col).^2+(0-row).^2;%weighted square of the distance to the starting point(lower of the left points)
end
[~,closest]=min(d02);
x0=col(closest);y0=row(closest);%closest point

if aud_near_ground>0.5%if there is sound near the ground
    dend2=(0-col).^2+(0-row).^2;%square of the distance to the ending point(lower left corner)
else
    dend2=1e6*(0-col).^2+(sz(1)-row).^2;%square of the distance to the ending point(higher of the leftmost points)
end

[~,closest]=min(dend2);
xend=col(closest);yend=row(closest);%closest point

[xl,yl]=draw_line(col,row,x0,y0,xend,yend);
xl=xl*simu.Rstep;yl=yl*simu.Zstep;

result.xr=xr;
result.xl=xl;

result.yr=yr;
result.yl=yl;

result.z=PEright.Z;
result.r=PEright.R;
result.Zsource=source.Zsource;
result.Name=p.Results.Name;
% %% Plotting(for debug)
%%
% a

end


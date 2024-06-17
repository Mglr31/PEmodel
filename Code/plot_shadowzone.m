function [f,f2, shadowzone]=plot_shadowzone(refractedres,homogeneousres,varargin)
%f=plot_shadowzone(refractedres,homogeneousres)
% A function that plot the shadow zone profile, ie the curve under wich the
% spl from the refracted PE simu is lower than the one from the homogeneous
% PE simu. ALso plot its intensity, ie the delta in dB.
%
%Inputs:
%   refractedres,homogeneousres: A structure containing all the results of a PE simulation. This stucture have a length of 1 or 2
%
%   Required fields:
%       -res(i).SPL: the grid of sound pressure levels(db?)
%
%       -res(i).Z: the vector of z coordinates(altitude,m)
%
%       -res(i).R: the vector of r coordinates(range,m)
%
%
%       -res(i).Ztop: the maximum height of the physical
%                   simulation(excluding the absorbing layer, in m)
%
%
% Optional input
%   'Source', source: a structure describing the source, so it can be plotted
%               on the figures
%
%   'Delta': by default equal to 0dB , described how many dB smaller is the
%               refractive Sound field at the shadow zone upper limit.




%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"refractedres");
addRequired(p,"homogeneousres");
default_source=[];
addParameter(p,'Source',default_source);
default_delta=0;
addParameter(p,'Delta',default_delta);

parse(p,refractedres,homogeneousres,varargin{:})
%%%%%%%%
source=p.Results.Source;
delta=p.Results.Delta;


spl=flipud(refractedres.SPL+delta-homogeneousres.SPL);%field of teh spl difference between refra and homo
intensity=spl;
R=refractedres.R;
Z=refractedres.Z;
n=length(R);
mat=spl>=0;
profileindex=zeros(1,n);
for i=1:n
    index=find(mat(:,i),1,'first');%find the first position for which mat=1,ie we are no more in the SZ
    if isempty(index)
        index=1;
    end
    if sum(mat(:,i))==0%if no point where refra is above homo
        index=length(mat(:,i));%index=index_max, SZ extent is to the top
    end
    profileindex(i)=index;
    intensity(profileindex(i)+1:end,i)=0;
end

profile= Z(profileindex);
profile(1)=0;%no SZ just under teh source/fixing some bug

f=figure;
plot(R,profile,'LineWidth',2,'DisplayName',"Shadow zone")
hold on
if ~isempty(source)
    plot(source.Xsource,source.Zsource,'*','MarkerSize',10,'DisplayName',"Sound source")
end


legend show
grid on
xlabel("Range (m)")
ylabel("Height (m)")
ylim([0,refractedres.Ztop])

set(gca,'FontSize',28)

shadowzone.R=R;
shadowzone.Z=Z;
shadowzone.profile=profile;
intensity(intensity==0)=NaN;
shadowzone.intensity=intensity;

f2=figure;
surf(shadowzone.R, shadowzone.Z, abs(intensity), 'EdgeColor', 'none');
hold on
plot3(R,profile,1000*ones(1,length(R)),'k','LineWidth',4);
if ~isempty(source)
    plot3(source.Xsource,source.Zsource,1000,'b*','MarkerSize',15,'LineWidth',3,'DisplayName',"Sound source")
end
shading(gca, 'flat');
view(2);
%colormap(turbo);
h=colorbar;

clim([0 max(abs(intensity),[],"all")]);
ylim([0 refractedres.Ztop])
xlabel("Range (m)")
ylabel("Height (m)")
set(gca,'FontSize',28)
ylabel(h, 'Refraction losses (dB)','FontSize',24)
grid off
title(refractedres.Name)



end

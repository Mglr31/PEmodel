function atm=makeMartianAtm(varargin)
%atm=makeMartianAtm(varargin)
%A function that make the atm struct corresponding to the Ls,
%LTST,Site,Longitude,Latitude and altitude passed as argument.
%Can alos pass a simu with a profile as an argument.
%% Input parser
p = inputParser;

%Optionnal input
default_simu=[];
addParameter(p,'Simu',default_simu);
default_LTST=12;
addParameter(p,'LTST',default_LTST);
default_Ls=180;
addParameter(p,'Ls',default_Ls);
default_latitude=18.45;
addParameter(p,'Latitude',default_latitude);
default_longitude=77.4;
addParameter(p,'Longitude',default_longitude);
default_site="none";
addParameter(p,'Site',default_site);
default_altitude=1;
addParameter(p,'Altitude',default_altitude);

parse(p,varargin{:});

ltst=p.Results.LTST;
ls=p.Results.Ls;
long=p.Results.Longitude;
lat=p.Results.Latitude;
site=p.Results.Site;
alt=p.Results.Altitude;

if site=="Perseverance"
    lat=18.45;log=77.4;
end


%%
if ~isempty(p.Results.Simu)
    titre=p.Results.Simu.Profile;
else

titre=strcat("profile",string(ltst),"_",string(lat),...
            " _",string(long),"_",string(ls),".csv");
end
try
    pr=readtable_header(titre);
catch
    disp("No corresponding profile found")
    return
end

load("modelMars.mat");
atm=modelMars;

atm.T=interp1(pr.Altitude,pr.t,alt);
atm.P=interp1(pr.Altitude,pr.p,alt);
try
atm.MolarFrac=[interp1(pr.Altitude,pr.vmr_n2,alt),...
    interp1(pr.Altitude,pr.vmr_co2,alt),...
    interp1(pr.Altitude,pr.vmr_ar,alt),...
    interp1(pr.Altitude,pr.vmr_o2,alt)];
catch
    atm.MolarFrac=[interp1(pr.Altitude,pr.n2,alt),...
    interp1(pr.Altitude,pr.co2,alt),...
    interp1(pr.Altitude,pr.ar,alt),...
    interp1(pr.Altitude,pr.o2,alt)];
end
atm.model_csv=[];
atm.experiment_csv=[];
atm.description="Atmosphere "+pr.Properties.Description(9:end-1)+" and an altitude of "+string(alt)+"m";
end
figure
surf(rvect, zvect, flipud(splgrid), 'EdgeColor', 'none')
%surf(x, y, flipud(z), 'EdgeColor', 'none')
colormap(turbo);
h=colorbar;
ylabel(h, ' (log(Pa))','FontSize',14)
caxis([50 150]);
ylim([0 simu.Zabs])
shading(gca, 'flat');
view(2)


%%
simu.Xmax=100;
source.freq=84;
simu.ZGround=4.3959 + 4.3918i;
simuconstant.Xmax=1000;
simurefra.Xmax=1000;
result1=ParabolicEquation(simu0,source,'Name',"0h LTST");
result=ParabolicEquation(simu12,source,'Name',"12h LTST");
results(3)=ParabolicEquation(simurefra,source,'Name',"very refractive atm");
results(4)=ParabolicEquation(simuconstant,source,'Name',"non refractive atm");

%%
plot_PE(r,["mu"])
plot_slice(result,[1.5 ])
plot_PE(result1,["spl"])
%%
simurefra.Zstep=0.1;
simurefra.Xmax=500;
simurefra
r=ParabolicEquation(simurefra,source,"Turbulence",turbulence);

r2=ParabolicEquation(simu,sr);
simuconstant.Xmax=500
r3=ParabolicEquation(simuconstant,source);
%% Real profile
source
load("simu12.mat");
simu12.Xmax=300;

r(1)=ParabolicEquation(simu12,source,'Name',"12h LTST");
simu12.WSpeed=1;
simu12.WDirection=180;
r(2)=ParabolicEquation(simu12,source,'Name',"12h LTST opposite wind ");

load("simu0.mat");
simu0.Xmax=300;

r(3)=ParabolicEquation(simu0,source,'Name',"0h LTST");
simu0.WSpeed=1;
simu0.WDirection=180;
r(4)=ParabolicEquation(simu0,source,'Name',"0h LTST opposite wind ");
simu12.WDirection=0;

r(5)=ParabolicEquation(simu12,source,'Name',"12h LTST favorable wind ");

%% plot

plot_PE(r(1),["spl"])
plot_PE(r(2),["spl"])
plot_PE(r(3),["spl"])
plot_PE(r,["spl"])

plot_slice(re,[15 ])

%%
r=ParabolicEquation(simuconstant,source,'Name',"PE constant");
rref=ParabolicEquation(simurefra,source,'Name',"Pe Refra");

randi=prousGroundAndi(simuconstant,source);
r2=randi;
r2.SPL=abs(r.SPL-randi.SPL);
plot_PE(r2,["spl"])

%% Section 3.4: Parabolic Equation refracting atm with wind 12 LTST
load("simu18.mat");
load("source.mat");

source.freq=100;
source.Amplitude=100;
source.Z=5;
simuWind=simu18;
simuWind.WSpeed=1;
simuBackWind=simuWind;
simuBackWind.WDirection=180;
r18=ParabolicEquation(simu18,source,'Name',"No wind");
rwind=ParabolicEquation(simuWind,source,'Name',"Forward wind");
rbwind=ParabolicEquation(simuBackWind,source,'Name',"Backward wind");
%%
plot_PE(r18,["spl"])
ca=gca;
set(ca.YAxis ,'visible','off')
plot_PE(r18,["windprofile"])
xlim([-6 6])

plot_PE(rwind,["spl"])
ca=gca;
set(ca.YAxis ,'visible','off')
plot_PE(rwind,["windprofile"])
xlim([-6 6])

plot_PE(rbwind,["spl"])
ca=gca;
set(ca.YAxis ,'visible','off')
plot_PE(rbwind,["windprofile"])
xlim([-6 6])

%%
res(1)=r18;
res(2)=rwind;
res(3)=rbwind;
f=plot_slice(res,[1.5]);
close(f)

%% Section 3.5: Parabolic Equation refracting turbulent atm  12 LTST
load("simu12.mat");
load("source.mat");
load("turbulence.mat")
%Changing simulation parameter
simu=simu12;
simu.Zabs=30;
simu.Zmax=50;
simu.Xmax=300;
%Source
source.freq=100;
source.Amplitude=100;
source.Zsource=5;
%Making turbulence
profile=readtable(simu.Profile);
profile(profile.Altitude>simu.Zabs,:)=[];%not shown in the result
c0=mean(profile.C_100Hz);
T0=mean(profile.t);
turbulence.c0=c0;
turbulence.T0=T0;
%T0=260;
kmin=1/30;%size of the biggest eddies
turbulence.L=1/kmin;
Turbulence.Lt=L;
turbulence.Lv=L;
kmax=1/7e-3;%7mm: kolmogorov scale on Mars(Petrosyan);
turbulence.wavenumbers=logspace(floor(log10(kmin)),floor(log10(kmax)),100);

%turbulence.wavenumbers=linspace(kmin,kmax,100);


turbulence.type="karmanwind";
turbulence.NofRealisation=20;

sigmaT=9;%air T std: 3K
turbulence.sigmaT=sigmaT;
sigmaV=0;%3m/s
turbulence.sigmaV=sigmaV;
turbulence.mu0=0.5*sigmaT/T0;

turbulence.a=0.1;%TBD

%Running simulation
r12=ParabolicEquation(simu,source,'Name',"Non turbulent atmosphere");
save("r12","r12");
r12av=ParabolicEquation(simu,source,'Name',"Hot and homogeneous atmosphere",'Averaged',1);
save("r12av","r12av");

rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent, a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_01_3_3=rturb;
save("rturbVKW_01_3_3","rturbVKW_01_3_3");
turbulence.a=1;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_1_3_3=rturb;
save("rturbVKW_1_3_3","rturbVKW_1_3_3");
turbulence.a=3;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_3_3_3=rturb;
save("rturbVKW_3_3_3","rturbVKW_3_3_3");

turbulence.type="karman";
sigmaT=6;%air T std: 3K
turbulence.sigmaT=sigmaT;
sigmaV=3;%3m/s
turbulence.sigmaV=sigmaV;
turbulence.mu0=0.5*sigmaT/T0;

turbulence.a=0.1;%TBD
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent, a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVK_01_6_0=rturb;
save("rturbVK_01_6_0","rturbVK_01_6_0");
turbulence.a=1;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVK_1_6_0=rturb;
save("rturbVK_1_6_0","rturbVK_1_6_0");
turbulence.a=3;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVK_3_6_0=rturb;
save("rturbVK_3_6_0","rturbVK_3_6_0");

turbulence.type="karmanwind";
sigmaT=6;%air T std: 3K
turbulence.sigmaT=sigmaT;
sigmaV=0;%3m/s
turbulence.sigmaV=sigmaV;
turbulence.mu0=0.5*sigmaT/T0;
turbulence.a=0.1;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent, a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_01_6_0=rturb;
save("rturbVKW_01_6_0","rturbVKW_01_6_0");
turbulence.a=1;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_1_6_0=rturb;
save("rturbVKW_1_6_0","rturbVKW_1_6_0");
turbulence.a=3;
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence,'Name',strcat("Turbulent(",turbulence.type,") a= ",string(turbulence.a), ", mu0= ",string(turbulence.mu0),", nb real: ",string(turbulence.NofRealisation)));
rturbVKW_3_6_0=rturb;
save("rturbVKW_3_6_0","rturbVKW_3_6_0");
%%
plot_PE(r12,["spl"])
ylim([0 20])



plot_PE(rturb2,["mu","spl"])
ylim([0 20])




close all
%%
clear rtot
rtot(1)=r12av;
rtot(2)=r12;
rtot(3)=rturbVK_01_6_0;
rtot(4)=rturbVK_1_6_0;
rtot(5)=rturbVK_3_6_0;
rtot(6)=rturbVKW_01_6_0;
rtot(7)=rturbVKW_1_6_0;
rtot(8)=rturbVKW_3_6_0;
rtot(9)=rturbVKW_01_3_3;
rtot(10)=rturbVKW_1_3_3;
rtot(11)=rturbVKW_3_3_3;
f=plot_slice(rtot,[1.5]);

%%





res(1)=r12;
res(2)=rturb;
res(3)=r12av;
res(4)=rturb2;
r(1)=r12av;
r(2)=r12;
r(3)=rturb;
f=plot_slice(r,[1.8],"AllSPL",1);
%%
f=figure
loglog(t.wavenumbers,F(2,:),'DisplayName',string(z(2)))
hold on
loglog(t.wavenumbers,F(100,:),'DisplayName',string(z(100)))
hold on
loglog(t.wavenumbers,F(200,:),'DisplayName',string(z(200)))
grid on 
legend show






%Changing simulation parameter
simu=simu12;
simu.Zabs=30;
simu.Zmax=50;
simu.Xmax=300;
%adding wind
simu.WDirection=180;
simu.WSpeed=1;
%Source
source.freq=100;
source.Amplitude=100;
source.Zsource=3;
%Making turbulence
z=0:simu.Zstep:simu.Zmax;

turbulence=makeTurbulence(simu.Profile,z,20);


turbulence.type="karmanwind";
turbulence.NofRealisation=10;
turb145=makeTurbulence(simu.Profile,z,20);
turb145.type="karmanwind";
turb145.sigmaT=4;%Asier value
turb145.sigmaV=3;%Viúdez-Moreiras value
turb145.NofRealisation=10;
%%
rturb=ParabolicEquation(simu,source,"Turbulence",turb145,'Name',"Turbulent atmosphere, z-dependent turbulence");
%%

plot_PE(rturb,["spl","mu"])
ylim([0 20])

%% PE
load("simu12.mat");
load("source.mat");
load("turbulence.mat")
%Changing simulation parameter
simu=simu12;
simu.Zabs=30;
simu.Zmax=50;
simu.Xmax=300;
%adding wind
simu.WDirection=180;
simu.WSpeed=1;
%Source
source.freq=100;
source.Amplitude=100;
source.Zsource=3;
%Making turbulence
z=0:simu.Zstep:simu.Zmax;

%Homogeneous turbulent parameter at z=10m
turbulence10=makeTurbulence(simu.Profile,10,20);
turbulence10.type="karmanwind";
turbulence10.NofRealisation=20;


%Inhomogeneous turbulent parameter at z=1.45m sigma=max(sigma in situ)
[~,index]=min(abs(z-1.45));%index at wivh z=1.45m
turb145=makeTurbulence(simu.Profile,10,20);
turb145.type="karmanwind";
% turb145.sigmaT=4/turb145.sigmaT(index)*turb145.sigmaT;%Asier value
% turb145.sigmaV=3/turb145.sigmaV(index)*turb145.sigmaV;%Viúdez-Moreiras value
turb145.NofRealisation=100;
% turb145=remakeTurbulence(turb145);

%Inhomogeneous turbulent parameter 0<z<20m
turbinhomo=makeTurbulence(simu.Profile,10,20);
turbinhomo.type="VKostachev";
turbinhomo.NofRealisation=100;

turb10=turb145;
turb10.sigmaT=10*turb10.sigmaT;
turb10.sigmaV=10*turb10.sigmaV;

%% Running simulation
tr12=ParabolicEquation(simu,source,'Name',"Refractive and non turbulent atmosphere");
tr12av=ParabolicEquation(simu,source,'Name',"Non refractive and non turbulent atmosphere",'Averaged',1);
trturbS=ParabolicEquation(simu,source,"Turbulence",turb145,'Name',"Salomon");
trturb10=ParabolicEquation(simu,source,"Turbulence",turb10,'Name',"10*salomon");

%%
trturbO=ParabolicEquation(simu,source,"Turbulence",turbinhomo,'Name',"Ostachev");
%%

tres(1)=tr12av;
tres(2)=tr12;
tres(3)=trturbS;
tres(4)=trturbO;

%%

f=plot_slice(tres,[1.5]);

%% Testing the effect of incresing the spectral resolution on the computation time
L=20;
turb100=makeTurbulence(simu.Profile,z,20);
turb100.type="karmanwind";
turb100.NofRealisation=5;
kmin=1/L;
kmax=1/7e-3;%7mm: kolmogorov scale on Mars(Petrosyan);

turb100.wavenumbers=linspace(kmin,kmax,100);
r100=ParabolicEquation(simu,source,"Turbulence",turb100,'Name',"100");
turb200=turb100;
turb200.wavenumbers=linspace(kmin,kmax,200);
r200=ParabolicEquation(simu,source,"Turbulence",turb200,'Name',"200");

turb500=turb100;
turb500.wavenumbers=linspace(kmin,kmax,500);
r500=ParabolicEquation(simu,source,"Turbulence",turb500,'Name',"500");


turb1000=turb100;
turb1000.wavenumbers=linspace(kmin,kmax,1000);
r1000=ParabolicEquation(simu,source,"Turbulence",turb1000,'Name',"1000");

%%
figure
x=[100 200 500 1000];
y=[r100.ComputationTime,r200.ComputationTime,r500.ComputationTime,r1000.ComputationTime]/5;

plot(x,y)
hold on
plot(x,0.0091*x+7.437)

p=polyfit(x,y,1);
%% Comparing computation time parrallel vs regular
rturb=ParabolicEquation(simu,source,"Turbulence",turbulence10,'Name',"Serial computing");
rturbparra=ParabolicEquationParra(simu,source,"Turbulence",turbulence10,'Name',"Parrallel computing");

res(1)=rturb;res(2)=rturbparra;
plot_slice(res,[1.5])


%% Comparing everything
%%% Inluence of kvect(100-1000-10 000 points)
kmax=1/7e-3;%7mm: kolmogorov scale on Mars(Petrosyan);

turb100=makeTurbulence(simu.Profile,10,20);
turb100.sigmaV=0;
turb100.type="VKostachev";
turb100.NofRealisation=20;


turb100.wavenumbers=linspace(0,kmax,100);

turb1000=turb100;
turb1000.wavenumbers=linspace(0,kmax,1000);

turb1000L=turb100;
turb1000L.wavenumbers=linspace(1/20,kmax,1000);

turb3000=turb100;
turb3000.wavenumbers=linspace(0,kmax,3000);

r100=ParabolicEquationParra(simu,source,"Turbulence",turb100,'Name',"100");
r1000=ParabolicEquationParra(simu,source,"Turbulence",turb1000,'Name',"1 000");
r3000=ParabolicEquationParra(simu,source,"Turbulence",turb3000,'Name',"10 000");

r1000L=ParabolicEquationParra(simu,source,"Turbulence",turb1000L,'Name',"1 000 L");

x=[100 1000 3000];
y=[r100.stdMu r1000.stdMu r3000.stdMu];
figure;plot(x,y/stdMu)
hold on
plot(1000, r1000L.stdMu/stdMu,'*')
grid on
xlabel("Number of modes")
ylabel("Std/expectedStd")
title("Ostachev spectra only T")

%%
r(1)=r100;r(2)=r1000;r(3)=r10000;

plot_slice(r,[1.5])

plot_PE(r(1),["randmu"])
plot_PE(r(2),["randmu",])
plot_PE(r(3),["randmu"])


allMumat=cat(3,r100.allMu{:});
stdMu100=std(allMumat,[],"all");
allMumat=cat(3,r1000.allMu{:});
stdMu1000=std(allMumat,[],"all");
allMumat=cat(3,r3000.allMu{:});
stdMu3000=std(allMumat,[],"all");

stdMut=0.5*turb100.sigmaT/turb100.T0;
stdMuv=turb100.sigmaV/turb100.c0;
stdMu=sqrt(stdMut^2+stdMuv^2);

%%
kmax=1/7e-3;%7mm: kolmogorov scale on Mars(Petrosyan);
nbp=100*(1:1:100);

for i=1:length(nbp)
    k=linspace(1e-3,kmax,nbp(i));
    dk=k(2);
    k2=logspace(-3,2,nbp(i));
    k2=[ k2 kmax];
    F=muspectrum2DKarmanwind(k,t.Lt,t.Lv,t.sigmaT,t.sigmaV,t.T0,t.c0,t.L);
    int(i)=dk*sum(F);
    F2=muspectrum2DKarmanwind(k2,t.Lt,t.Lv,t.sigmaT,t.sigmaV,t.T0,t.c0,t.L);
    int2(i)=trapz(k2,F2);
end
%%
figure
semilogx(nbp,int,'DisplayName',"Sum")
hold on
semilogx(nbp,int2,'DisplayName',"Trapeze")
xlabel("Number of point used to compute the integral")
ylabel("Value of the integral")
grid on
set(gca,'FontSize',24)
legend show

figure
plot(nbp,int./int2)
figure
semilogx(k2,F2)
grid on
k=linspace(1e-3,kmax,nbp(100));
xlim([0 ,2*k(2)])
%% Amplitude distribution

x=[10 50 100 200];
y=[637 77 3 2];
figure
loglog(x,y)
%%

%load("rturbinhomo500.mat");

m=length(rturbinhomo.R);
R=rturbinhomo.R;
a=zeros(1,m);
for i=1:m
    i
    var=exctract1pointSoundFluctuation(rturbinhomo,1.5,R(i));
    try
    fitd=fitdist(var.NormA','gamma');
    a(i)=fitd.a;
    catch
        a(i)=NaN;
    end

end

%%
figure
semilogy(R,a)
grid on
xlabel("Range(m)")
ylabel("a parameter in the gamma distribution")
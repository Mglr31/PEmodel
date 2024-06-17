load("simu6.mat");
load("simu12.mat");
load("source.mat");

source.freq=100;
source.Amplitude=100;
source.Zsource=5;
r6=ParabolicEquation(simu6,source,'Name',strcat("Z= ",string(simu6.ZGround)," | reference"));
%r12=ParabolicEquation(simu12,source,'Name',"Hot and refractive atmosphere");
%r12av=ParabolicEquation(simu12,source,'Name',"Hot and homogeneous atmosphere",'Averaged',1);

atm=getAtmFromProfile(simu6.Profile,source.freq);
atm12=getAtmFromProfile(simu12.Profile,source.freq);
Z=getGroundImpedance(cannonicalGround,atm,source.freq);
Z12=getGroundImpedance(cannonicalGround,atm12,source.freq);

ground=cannonicalGround;
ground.tort=0.6;
Ztort=getGroundImpedance(ground,atm,source.freq);
simuZ=simu6;
simuZ.ZGround=Z;

simuZ12=simu6;
simuZ12.ZGround=Z12;

simuZtort=simu6;
simuZtort.ZGround=Ztort;

rZ=ParabolicEquation(simuZ,source,'Name',strcat("Z= ",string(simuZ.ZGround)," | 6h"));
rZ12=ParabolicEquation(simuZ12,source,'Name',strcat("Z= ",string(simuZ12.ZGround),"| 12h "));
rZtort=ParabolicEquation(simuZtort,source,'Name',strcat("Z= ",string(simuZtort.ZGround),"| tort=0.6 "));

%%

plot_PE(r6,["spl"])
ylim([0 20])

plot_PE(rZ,["spl"])
ylim([0 20])

plot_PE(rZ12,["spl"])
ylim([0 20])

plot_PE(rZtort,["spl"])
ylim([0 20])
%%
rcomp=r6;
rcomp.SPL=abs(r6.SPL-rZ.SPL);
plot_PE(rcomp,["spl"])
clim([0 2])
ylim([0 20])
rcomp2=rZ12;
rcomp2.SPL=abs(r6.SPL-rZ12.SPL);
plot_PE(rcomp2,["spl"])
clim([0 2])
ylim([0 20])

rcomp3=rZtort;
rcomp3.SPL=abs(r6.SPL-rZtort.SPL);
plot_PE(rcomp3,["spl"])
clim([0 2])
ylim([0 20])

rcomp4=rZtort;
rcomp3.SPL=abs(r6.SPL-rZtort.SPL);
plot_PE(rcomp3,["spl"])
clim([0 2])
ylim([0 20])
%%
r(1)=r6;
r(2)=rZ;
r(3)=rZ12;

plot_slice(r,[5])

%%
randi=prousGroundAndi(simu6,source);
randi06=prousGroundAndi(simu6,source);

plot_PE(randi,["spl"])
ylim([0 20])
plot_PE(randi06,["spl"])
ylim([0 20])

rcomp=randi;
rcomp.SPL=abs(randi.SPL-randi06.SPL);
plot_PE(rcomp,["spl"])
clim([0 10])
ylim([0 20])

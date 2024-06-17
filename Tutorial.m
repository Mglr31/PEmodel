%simulation grid with atmospheric data representing the Martian atmosphere
%at 12 local time at Perseverance site
load("simu12.mat")
simu12.Zabs=20;%Max alt with realistic results
simu12.Zmax=40;%Max alt of the grid
simu12.WSpeed=1;%adding wind to the simulation
%84 Hz, 80 dB source at 5 m above the ground 
load("source.mat")
source.Amplitude=100;%changing amplitude to 100 dB
source.freq=100;
%% Lauching PE simulation
pe=ParabolicEquation(simu12,source);

%% Displaying the resulting sound field
plot_PE(pe,["spl"])

%% Displaying a slice of SPL at 3 m
plot_slice(pe,3)

%% Separate the different effects responsible for the losses
effects=separateEffectsPE(simu12,source,1.5,'Turbulence',0);%no effect of turbulence
plot_effectsPE(effects,1.5,"allnosum");
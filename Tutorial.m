%Simulation grid with atmospheric data representing the Martian atmosphere
%at 12 local time at Perseverance site
load("simu12.mat")
simu12.Zabs=20;%Max alttitude of the sound field
simu12.Zmax=40;%Max alttitude of the grid
simu12.WSpeed=1;%adding wind to the simulation
simu12.ZGround=10+10i;%Ground acoustic impedance
%loading the source
load("source.mat")
source.Amplitude=100;%changing amplitude to 100 dB
source.freq=100;
%% Lauching PE simulation
pe=ParabolicEquation(simu12,source);

%% Displaying the resulting sound field
plot_PE(pe,["spl"])

%% Displaying a slice of SPL at 3 m
plot_slice(pe,3)




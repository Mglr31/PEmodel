function atm=getAtmFromProfile(profile_name,freq)
%atm=getAtmFromProfiel(profile_name)
%A function that return an atm structure(used in the getGroundImpedance
%function) from the MCD profile "profile_name" for the frequency freq(Hz)

pr=readtable_header(profile_name);
atm.T=pr.t(1);
atm.P=pr.p(1);
atm.listOfMole=["N2"		"CO2"	"Ar"	"O2"];
atm.MolarFrac=[pr.n2(1)	pr.co2(1)	pr.ar(1) pr.o2(1)];
atm.listOfMode=["nu"	"nu1"	"none"	"nu"];

power10=round(log10(freq)); % order of magnitude of the sound source frequency (needed to choose the adequate sos )
switch power10
    case 0
        atm.c=pr.C_1Hz(1);
    case 1
        atm.c=pr.C_10Hz(1);
    case 2
        atm.c=pr.C_100Hz(1);
    case 3
        atm.c=pr.C_1000Hz(1);

    case 4
       atm.c=pr.C_10000Hz(1);
end



end
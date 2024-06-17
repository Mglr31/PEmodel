%sources heights
heights=[1 5 10 20 30 40 50 80 100 130 160];
heights=linspace(1,150,20)
%rover coordinates
r=80;
z=1.8;

load("simu12.mat")
simu=simu12;

load("source.mat")
source.freq=100;
source.Amplitude=100;
n=length(heights);
simu.Xmax=r+1;
%%

for i=1:n
    
    h=heights(i);
    source.Zsource=h;
    if simu.Zmax<1.3*h
        simu.Zmax=h*1.5;
        simu.Zabs=simu.Zmax*0.8;
    end
        disp(simu.Zmax)
        disp(h)
        tic
        pe=ParabolicEquation(simu,source);
        toc
        spl(i)=get_spl(pe,r,z);
        i
   
end


%% Plotting

figure
plot(spl,heights,'LineWidth',2,'HandleVisibility','off')
xline(5,'LineWidth',2,'DisplayName',"Audibility threshold")
ylabel("Source height(m)")
xlabel("Received sound(dB)")
legend show
set(gca,'FontSize',22)

%% Audibility area
heights=[1  10  50 150];
%rover coordinates


load("simu12.mat")
simu=simu12;

load("source.mat")
source.freq=100;
source.Amplitude=100;
n=length(heights);
%%

for i=1:n
    
    h=heights(i);
    source.Zsource=h;
    if simu.Zmax<1.3*h
        simu.Zmax=h*1.5;
        simu.Zabs=simu.Zmax*0.8;
    end
        disp(simu.Zmax)
        disp(h)
        tic
        results(i)=get_audible_range(simu,source,20,'Name',"Source height: "+string(h)+" m");
        toc
        i
   
end

%%
plot_audible_range(a)


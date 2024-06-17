function f=plot_PE(result,plotOption)
%plot_PE(result,plotOption)
%A function that plot the results of a Parabolic Equation simulation
%Inputs:
%   result: A structure containing all the results of a PE simulation(it is
%           the output of the function ParabolicEquation).
%   Required fields:
%       result.SPL: the grid of sound pressure levels(db?)
%
%       result.P: the grid of complex amplitude of the pressure(Pa?)
%
%       result.K : the grid of wavenumber
%
%       result.Z: the vector of z coordinates(altitude,m)
%
%       result.R: the vector of r coordinates(range,m)
%
%       result.Freq_Source: the frequency of the source that was simulated(Hz)
%
%       result.Ztop: the maximum height of the physical
%                   simulation(excluding the absorbing layer, in m)
%
%   plotOption: an array of string defining the figure that are to be
%               plotted.(Ex:plotOption=["spl","c"])
%   Possible values:
%
%       "spl": The grid of the sound pressure level(db)
%
%       "p": The grid of the amplitude of the pressure(abs of the complex
%            amplitude,Pa)
%
%       "phase" : The grid of the phase of the pressure(angle of the complex
%            amplitude,deg)
%
%       "c": the grid of the effective speed of sound(inverse of real part of the wave
%               number divide by the pulsation; k=omega/c+i*alpha; in m/s)
%
%       "alpha": the grid of the attenuation coefficient used in the
%               simulation(imaginary part of the wave
%               number;k=omega/c+i*alpha; in m-1;
%
%       "cprofile": the profile of speed of sound(in m/s) works for layered
%                   atmosphere only
%
%       "tprofile": the profile of temperature(K)
%
%       "windprofile": the profile of wind speed(m/s)
%
%       "allSPL":  All the grids of the sound pressure level(case of a turbulent atm,db)
%
%       "mu":  The grid of the refractive index mu(case of a turbulent atm, no unit )
%
%      "randmu":  A random grid of the refractive index mu (case of a turbulent atm,  no unit )
%
%       "SI" : the grid of teh scintillation index (case of a turbulent atm,  no unit)
%
%       "lnpVariance" : the grid of the variance of the log-normal pressure index (case of a turbulent atm,  no unit)
if plotOption=="all"
    plotOption=["spl","p","alpha","c","cprofile"];
end
n=length(plotOption);
z=result.Z;
ztop=result.Ztop;
r=result.R;


for i =1:n
    f=figure();
    if plotOption(i)=="spl"
        surf(r, z, flipud(result.SPL), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' SPL (db)','FontSize',24)
        caxis([-10 90]);
        %caxis([result.SPLmin 100])
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        if isfield(result,'Name')
            title(result.Name)
        end
        if isfield(result,'Mic')
            m=result.Mic;
            hold on
            plot3([0],[3],500,'b*','MarkerSize',20,'LineWidth',3,'DisplayName',"Source")
            hold on
            plot3([m.X],[m.Z],500,'k+','MarkerSize',15,'LineWidth',3,'DisplayName','Microphone')

        end
        if isfield(result,'shadowzone')
            sz=result.shadowzone;
            if ~isempty(sz)
                hold on
                plot3([sz.R],[sz.profile],500*ones(1,length(sz.R)),'Color',[0.7 0.7 0.7],'MarkerSize',20,'LineWidth',4,'DisplayName',"Source")
                hold on
            end
        end

    elseif plotOption(i)=="p"
        surf(r, z, flipud(abs(result.P)), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' Pressure amplitude (Pa)','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title(result.Name)
    elseif plotOption(i)=="phase"
        surf(r, z, flipud(angle(result.P)*180/pi), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' Complex pressure phase (deg)','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title(result.Name)

    elseif plotOption(i)=="wavefront"
        phase=flipud(angle(result.P)*180/pi);
        wavefront=phase<10&phase>-10;%phase=0
        surf(r, z, double(wavefront), 'EdgeColor', 'none')
        map=[1 1 1
            0 0 0];
        colormap(map);
        h=colorbar;
        set(h,'Visible','off')
        ylabel(h,' Wavefronts','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title(result.Name)
    elseif plotOption(i)=="c"
        omega=2*pi*result.Freq_Source;
        kreal=real(result.K);
        surf(r, z, flipud(omega*kreal.^(-1)), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' Effective speed of sound (m/s)','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")

    elseif plotOption(i)=="alpha"
        surf(r, z, flipud(imag(result.K)), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' Attenuation coefficient (m^{-1})','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")

    elseif plotOption(i)=="cprofile"
        omega=2*pi*result.Freq_Source;
        kreal=real(result.K);
        cprofile=flipud(omega*kreal(:,1).^(-1));
        plot(cprofile,z,'LineWidth',2)
        xlabel(" Effective speed of sound (m/s)")
        ylabel("Altitude (m)")
        grid on


    elseif plotOption(i)=="alphaprofile"
        %alpha profile used in the simulation including the upper (and sometiems lower) attenuation
        %layer
        kimag=imag(result.K);
        aprofile=flipud(kimag(:,1));
        plot(aprofile,z,'LineWidth',2)
        xlabel(" Attennuation coefficient (m^{-1})")
        ylabel("Altitude (m)")
        grid on


    elseif plotOption(i)=="tprofile"

        plot(result.profiles.T,z,'LineWidth',2)
        xlabel(" Air T (K)")
        ylabel("Altitude (m)")
        grid on
        ylim([0 ztop])

    elseif plotOption(i)=="windprofile"

        plot(result.profiles.Wind,z,'LineWidth',2)
        xlabel(" WS (m/s)")
        ylabel("Altitude (m)")
        grid on
        ylim([0 ztop])
    elseif plotOption(i)=="allSPL"
        N=length(result.allSPL);
        for k=1:2
            figure
            surf(r, z, flipud(result.allSPL{k}), 'EdgeColor', 'none')
            colormap(turbo);
            h=colorbar;
            ylabel(h, ' SPL (db)','FontSize',22)
            caxis([-10 90]);
            ylim([0 ztop])
            shading(gca, 'flat');
            view(2)
            xlabel("Range (m)")
            ylabel("Altitude (m)")
        end

    elseif plotOption(i)=="mu"

        %surf(r, z, result.Mu/std(result.Mu(:)), 'EdgeColor', 'none')
        surf(r, z,abs(result.Mu), 'EdgeColor', 'none')

        colormap(turbo);
        h=colorbar;
        ylabel(h, ' |\mu|','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title({"Refractive index fluctuation";result.Name})
        title(result.Name)
    elseif plotOption(i)=="randmu"
        l=length(result.allMu);
        index=floor(rand(1)*(l-1))+1;
        %surf(r, z, result.Mu/std(result.Mu(:)), 'EdgeColor', 'none')

        Mu=result.allMu{index};
        surf(r, z,abs(Mu), 'EdgeColor', 'none')

        colormap(turbo);
        h=colorbar;
        ylabel(h, ' |\mu|','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range(m)")
        ylabel("Altitude(m)")
        title({"Refractive index fluctuation";result.Name})
        title(result.Name)

    elseif plotOption(i)=="SI"
        surf(r, z, flipud(result.SI), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, ' Scintillation index','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title(result.Name)
    elseif plotOption(i)=="lnpVariance"
        surf(r, z, flipud(result.sigmaX2), 'EdgeColor', 'none')
        colormap(turbo);
        h=colorbar;
        ylabel(h, 'Log-normal pressure variance','FontSize',22)
        ylim([0 ztop])
        shading(gca, 'flat');
        view(2)
        xlabel("Range (m)")
        ylabel("Altitude (m)")
        title(result.Name)
    end
end

set(gca,'FontSize',26)
end


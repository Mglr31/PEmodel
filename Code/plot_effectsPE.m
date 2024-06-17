function f=plot_effectsPE(results,height,plotOptions)
%f=plot_effectsPE(results,plotOptions)
%A function that plots the effects conatined in results, an output of the
%function separateEffectsPE
%

n=length(results);
res=[];
for i=1:n
    if results(i).h==height
        res=results(i);
        break
    end
end

if isempty(res)
    disp("No height corresponding the the required height in the results structure")
    return

else

    if plotOptions=="all"
        plotOptions=["ground","geometric","attenuation",...
            "temperature","wind","turbulence","complete","sum"];
    end
    if plotOptions=="all0"
        plotOptions=["ground","geometric","attenuation",...
            "temperature","wind0","turbulence","complete","sum"];
    end
    if plotOptions=="allnosum"
        plotOptions=["ground","geometric","attenuation",...
            "temperature","wind","turbulence","complete"];
    end
    f=figure;
    r=res.r;

    d=res.d;

    loc=["#0072BD","#D95319","#EDB120",	"#7E2F8E",	"#77AC30"	"#4DBEEE",	"#A2142F","#000000","#0072BD","#00FF00","#00FF00","#77AC30"];


    if ismember("ground",plotOptions)
        plot(r,-res.ground,'--','Color',loc(2),'LineWidth',2,'DisplayName','Ground')
        hold on
    end
    if ismember("geometric",plotOptions)
        plot(r,-20*log10(0.1./d),'Color',loc(1),'LineWidth',2,'DisplayName','Geometric spreading')
        hold on
    end
    if ismember("classical",plotOptions)
        plot(r,-res.classical,'Color',loc(2),'LineWidth',2,'DisplayName','Classical attenuation')
        hold on
    end
    if ismember("non classical",plotOptions)
        plot(r,-res.non_classical,'Color',loc(3),'LineWidth',2,'DisplayName','Non classical attenuation')
        hold on
    end
    if ismember("attenuation",plotOptions)
        plot(r,-res.non_classical-res.classical,'Color',loc(3),'LineWidth',2,'DisplayName','Attenuation')
        hold on
    end
    if ismember("temperature",plotOptions)
        plot(r,-res.temperature_pr,'--','Color',loc(4),'LineWidth',2,'DisplayName','Temperature profile')
        hold on
    end
    if ismember("wind",plotOptions)
        plot(r,-res.wind_pr,'--','Color',loc(5),'LineWidth',2,'DisplayName','Wind speed profile')
        hold on
    end
    if ismember("wind0",plotOptions)
        plot(r,-res.wind_pr0,'--','Color',loc(5),'LineWidth',2,'DisplayName','Wind speed profile')
        hold on
    end
    if ismember("turbulence",plotOptions)
        plot(r,-res.turbulence,'-.','Color',loc(6),'LineWidth',2,'DisplayName','Turbulence')
        hold on
    end
    if ismember("sum",plotOptions)
        if ismember("wind0",plotOptions)
            sum=-20*log10(0.1./d)-res.turbulence-res.wind_pr0-res.temperature_pr-res.non_classical-res.classical-res.ground;
        elseif ismember("wind",plotOptions)
            sum=-20*log10(0.1./d)-res.turbulence-res.wind_pr0-res.temperature_pr-res.non_classical-res.classical-res.ground;
        end
        plot(r,sum,'Color',loc(7),'LineWidth',3,'DisplayName','Sum of all causes')
        hold on
    end
    if ismember("complete",plotOptions)
        plot(r,res.complete,'Color',loc(8),'LineWidth',3,'DisplayName','Losses from full simulation')
        hold on
    end
end
ymax=res.complete(end-1);
ymin=min(-res.turbulence);

ylim([ymin-5 ymax+62])
xlim([0 max(r)])
legend('Location','northwest')
grid on
xlabel("Range (m)")
ylabel("Losses (dB)")
set(gca,'FontSize',24)

end

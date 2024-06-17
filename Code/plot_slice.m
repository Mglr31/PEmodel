function f=plot_slice(results,heights,varargin)
%plot_slice(results,heights)
% A function that plot slices of the Parabolic Equations simulation
% contained in results at all the heights listed in heights
%Inputs:
%   results: A structure containing all the results of a PE simulation. This stucture can have a length >1 if multiple results need to be plotted simultaneously(it is
%           the output of the function ParabolicEquation).
%   Required fields:
%       results(i).SPL: the grid of sound pressure levels(db?)
%
%       results(i).Z: the vector of z coordinates(altitude,m)
%
%       results(i).R: the vector of r coordinates(range,m)
%
%
%       results(i).Ztop: the maximum height of the physical
%                   simulation(excluding the absorbing layer, in m)
%
%       Optional fields:
%       results(i).Zsource: altitude of the source(m).
%
%      heights: an array of all the heights at which a slice must be
%               plotted(m).
%
% Optionnal input:
%       "AllSPL",1  : plot all the computed spl in grey in the case of a turbulent
%       simulation
%
%
%       "Alpha',1: plot the SPL for a simple exponential decay with tha
%       alpha at this height
%       

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"results");
addRequired(p,"heights");
default_as=0;
addParameter(p,'AllSPL',default_as);
default_alpha=0;
addParameter(p,'Alpha',default_alpha)
parse(p,results,heights,varargin{:})
%%%%%%%%

n=length(heights);
m=length(results);

for i=1:n% for each height

    f=figure();
    if isfield(results(1),"Zsource")%if the height of the source was given
        %computing where the wide angle approximation holds
        zs=results(i).Zsource;
        rvalid=abs(zs-heights(i))/tan(70*pi/180);
        xline(rvalid,'DisplayName',"Limit of the wide-angle approximation")
        hold on
    end
    for j=1:m%for each simulation
        spl=flipud(results(j).SPL);
        if heights(i)<=results(j).Ztop
            [~,index]=min(abs(heights(i)-results(j).Z));
            
            if p.Results.AllSPL && ~isempty(results(j).allSPL)
               

                h=length(results(j).allSPL);
                grey=[0.8,0.8,0.8,0.3];
                energ=zeros(h,length(results(j).SPL(1,:)));
                line=energ;
                for k=1:h
                    aspl=flipud(results(j).allSPL{k});
                    energ(k,:)=10.^(aspl(index,:)/10);
                    line(k,:)=aspl(index,:);
                    plot(results(j).R,aspl(index,:),'Color',grey,'LineWidth',2,'HandleVisibility',"off")
                    hold on
                end
                 plot(results(j).R,spl(index,:),'LineWidth',2,'DisplayName',"Average SPL")
            hold on
            %sigma=std(energ,[],1);
            sigma=std(line,[],1);%standard deviation of all the spl at each range (in dB)
           % e=10.^(spl(index,:)/10);
            plot(results(j).R,spl(index,:)+sigma,'--k','LineWidth',2,'DisplayName',"  +1\sigma and - 1\sigma")
            hold on
            plot(results(j).R,spl(index,:)-sigma,'--k','LineWidth',2,'HandleVisibility',"off")
            else
            plot(results(j).R,spl(index,:),'LineWidth',2,'DisplayName',results(j).Name)
            hold on
            if p.Results.Alpha
                K=flipud(results(j).K);
                alpha=imag(K(index,1));
                a0=results(j).source.Amplitude;
                alpha_db=-20*alpha/log(10);
                r=results(j).R;
                plot(r,a0+20*log10(0.1./r)+alpha_db*r,'LineWidth',2,'DisplayName',results(j).Name+" :att curve")
            end
            end
        else
            disp("The desired height is higher than the max altitude of the simulation. The slice was not plotted.")
        end
    end
    legend show
    grid on
    xlabel("Range (m)")
    ylabel("SPL (dB)")
    title(strcat("Slice at a height of ",string(heights(i))," m"))
    set(gca,'FontSize',26)


end







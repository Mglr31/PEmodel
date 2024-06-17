function f=plot_Vslice(results,ranges)
%plot_Vslice(results,ranges)
% A function that plot vertical slices of the Parabolic Equations simulation
% contained in results at all the heights listed in heights
%If the length of result is 2, it plots the difference between the first PE
%result and the second one(result(1)-result(2))
%Inputs:
%   results: A structure containing all the results of a PE simulation. This stucture have a length of 1 or 2 
%           
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
%      ranges: an array of all the ranges at which a vertical slice must be
%               plotted(m).


n=length(ranges);
m=length(results);
f=figure();
for i=1:n% for each range

%     
%     if isfield(results(i),"Zsource")%if the height of the source was given
%         %computing where the wide angle approximation holds
%         zs=results(i).Zsource;
%         rvalid=abs(zs-heights(i))/tan(70*pi/180);
%         xline(rvalid,'DisplayName',"Limit of the wide-angle approximation")
%         hold on
%     end

    if m==1
        spl=flipud(results.SPL);
    elseif m==2
        spl=flipud(results(1).SPL-results(2).SPL);
    end

        
            [~,index]=min(abs(ranges(i)-results(1).R));
            plot(spl(:,index),results(1).Z,'LineWidth',3,'DisplayName',strcat("At a range of ",string(ranges(i))," m"))
            hold on
        
    end
    legend show
    legend('Location','northwest')
    grid on
    if m==1
        xlabel("SPL(dB)")
    elseif m==2
        xlabel("Delta SPL (dB)")
    end

    
    ylabel("Height(m)")
    ylim([0,results(1).Ztop])
    %title("Vertical slices")
    set(gca,'FontSize',26)


end







function slice=extract_slice(results,height,varargin)
%slice=extract_slice(results,heights,varargin)
% A function that extract slices of the Parabolic Equations simulation
% contained in results at the altitude "height" in m 
%Output:
%   slice: A structure containing the slices of teh PE simulations.
%   Required fields:
%       slice(i).SPL: the spl slice for results(i)
%
%
%       slice(i).R: the vector of r coordinates(range,m) for results(i)
%
%
%
      

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"results");
addRequired(p,"height");


parse(p,results,height,varargin{:})
%%%%%%%%

m=length(results);


    for j=1:m%for each simulation
        spl=flipud(results(j).SPL);
        if height<=results(j).Ztop
            [~,index]=min(abs(height-results(j).Z));
            slice(j).R=results(j).R;
            slice(j).SPL=spl(index,:);
        end
    end


end







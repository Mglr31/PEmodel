function result=quantify_turbulence_effect(rTurb,rNonTurb,rNonRefra,height)
%result=quantify_tubulece_effect(rTurb,rNonTurb,rNonRefra)
%
%A fucntion taht quantify the effect of the turbulence on the reduction of
%the shadow zone depth at a specific height.
%
%Inputs:
%   
%   rTurb : The PE simulation in the refractive+turbulent case 
%   
%   rNonTurb : The PE simulation in the refractive+non-turbulent case 
%   
%   rNonRefra : The PE simulation in the non-refractive+non-tubrbulent case 
%
%   -> the 3 simulations must have been made on the same grid
%
%   height : the height at which the SZ is mesured
%
%
%Outputs: reuslt a structure with the fields
%
%   .SZdepth : the depth of the SZ in the non-turbulent case 
%
%   .Is_negSZ : 1 is the depth of the nonturbSZ is negative, 0 else
%
%   .TurbSZdepth : the depth of the SZ in the turbulent case 
%
%   .SZreduction : reduction of the SZ depth brought by teh turbulence in %


turb=extract_slice(rTurb,height);
Nturb=extract_slice(rNonTurb,height);
Nrefra=extract_slice(rNonRefra,height);

result.R=turb.R;
result.SZdepth=Nrefra.SPL-Nturb.SPL;
result.Is_negSZ=(Nrefra.SPL-Nturb.SPL)<=0;
result.TurbSZdepth=Nrefra.SPL-turb.SPL;
result.SZreduction=(result.SZdepth-result.TurbSZdepth)./result.SZdepth*100;


end

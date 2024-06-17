function spl=get_spl(result,r,z)
%spl=get_spl(result,r,z)
%
%Get the spl at a raneg of r meters and a height of z meters i the PE
%simulation contained in result.

zvect=result.Z;
rvect=result.R;
[~,zidx]=min(abs(zvect-z));
[~,ridx]=min(abs(rvect-r));

SPL=result.SPL;
sz=size(SPL);
spl=SPL(sz(1)-zidx,ridx);




end
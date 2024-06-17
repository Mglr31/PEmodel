function result=ParabolicEquation(simu,source,varargin)
%result=ParabolicEquation(simu,source)
%A function that apply the Parabolic Equation model to return a field of SPL in th e(z,r) plan.
% Inputs:
%   simu: a structure containing the parameter of the simulation:
%    Required fields:
%       simu.Zmax:       the maximum altitude of the simulation(m)
%
%       simu.Xmax:       the maximum distance of the simulation in the x(or r)
%                         direction(m)
%
%       simu.Rstep:       the length of a step in  the x(or r) direction(m)
%
%
%       simu.Zstep:        The height the vertical layers(m)
%
%       simu.Zabs:        The height at which the absorbing layer start(m)
%
%       simu.Profile:       The name of a csv file containing
%                             a complete profile to use
%
%       simu.ZGround:        The normalized impedance of the ground
%
%    source: a structure containing the parameter of the sound source:
%    Required fields:
%       source.freq:        The frequency of the sound emmited by the
%                           source(Hz)
%       source.Zsource:        The height of the sound source(m)
%
%       source.Amplitude:      The amplitude at 10cm of the source(dB)
%
%
%Optional inputs:
%    mic: a structure containing the positions of the sound receptors:
%    Required fields:
%       mic.R: an array containing the positions of the mics along the r axis(m)
%
%       mic.Z:an array containing the positions of the mics along the z axis(m)
%
%   'Name',name: the name of the simulation(string)
%
%   'Turbulence',turbulence: a structure representing the kind of
%                            atmospheric turbulence that will be added to
%                            the simulation. Default value:[]: no
%                            turbulence.
%
%               Required fields:
%
%                   turbulence.wavenumbers: an array of the wavenumbers used to compute the turbulent field
%
%                   turbulence.type: the type of turbulence. Possible values:"gaussian","karman","experimental"
%                                      If the type is gaussian or karman,
%                                      turbulence needs to have two
%                                      additional fields:mu0 and a.
%                                      If the type is experimental,
%                                      turbulence need to have the field
%                                      F,F being the experimental spectral
%                                      density of mu.
%
%                  turbulence.NofRealisation: the number of random
%                                           realisation of the turbument fields that will be
%                                           computed for this simulation. The result is the
%                                           logarithmic mean of these realisation.
%
%
%       'Averaged':     default 0, if equal to 1 the atmospheric profile will
%                       be averaged to suppress the effect of sound
%                       refraction,if 2 only the T profile will be averaged
%
%       'AverageLim'    default 1000, the maximum height to take into account for the averaging 
%
%
%
%       "Light":    default 0, if equal to 1 mu fields wont be saved and
%                   spl will be saved in single format so the result is lighter
%
%      "Resolution":    default 0, if different than zero the resolution of
%                the output SPL grid is diminished "Resolution" times. Ex:
%                if one wants to compute a millimetric grid but store a
%                centimetric grid: "Resolution",10
%
%       "NoGround": by default equal to 0, equal to 1 if one want a
%       simulation without ground
%
%
%
% Outputs:
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
%
%   Optionnal fields (if Turbulence,"turbulence")
%
%       result.allSPL
%
%       result.allMu
%

%%%%%%%%
%Input parser
p=inputParser();
addRequired(p,"simu");
addRequired(p,"source");
default_mic=[];
addParameter(p,'mic',default_mic);
default_name="PE simulation";
validname= @(x) isstring(x);
addParameter(p,'Name',default_name,validname);
default_turb=[];
addParameter(p,'Turbulence',default_turb);
default_av=0;
addParameter(p,'Averaged',default_av);
default_avl=1000;
addParameter(p,'AverageLim',default_avl);
default_light=0;
addParameter(p,'Light',default_light);
default_res=0;
addParameter(p,'Resolution',default_res);
default_ng=0;
addParameter(p,'NoGround',default_ng);
parse(p,simu,source,varargin{:})

%%%%%%%%

%adding ground impedance to the simu if not already added
if ~isfield(simu,"ZGround")
    simu.ZGround=getGroundImpedanceVBO(simu.Ground,makeMartianAtm('Simu',simu),source.freq);
end

light=p.Results.Light;
reso=p.Results.Resolution;
%Vector representing the z and r direction
zvect=0:simu.Zstep:simu.Zmax;
rvect=0:simu.Rstep:simu.Xmax;

simu.Averaged=p.Results.Averaged;
simu.AverageLim=p.Results.AverageLim;
simu.NoGround=p.Results.NoGround;
if simu.NoGround
    neg_zvect=-(simu.Zmax-simu.Zabs):simu.Zstep:0;
    zvect=[neg_zvect zvect(2:end)];
end
simu.zvect=zvect;
[kgrid,profiles]=make_kgridPE(simu,source);%grid of wavenumber
result.profiles=profiles;
s=size(kgrid);%grid dimension
m=s(1);
n=s(2);

psigrid=kgrid;
psigrid(:)=0;
[~,idx_ground]=min(abs(zvect));%index of the altitude z=0
ka=kgrid(m-idx_ground+1,1);% k at the ground(or maybe at a medium altitude??)
k0=kgrid(m-idx_ground+1,1);% k at the ground(m=ground)
init=make_initPE(simu,source,ka);
psigrid(:,1)=flipud(init);

%Computing all the matrices that will be used in the matricial equation of propagation
alpha=0.5*1i/ka;
beta=flipud(0.5*1i*(kgrid(:,2).^2-ka^2)./ka);%vector of beta for the layered(non turbulent case);flipped because in kgrid the 1 index correspond to the top of the grid ie the upper atm whereas in the convention that is used in all this code index 1 is the ground
gamma=alpha/simu.Zstep^2;
D=diag(beta);

%ground condition
sigma2=-(3-2*1i*k0*simu.Zstep/simu.ZGround)^-1;
sigma1=-4*sigma2;
%Upper boundary condition
tau2=-(3+2*1i*k0*simu.Zstep)^-1;
tau1=-4*tau2;

c=zeros(1,m);
c(1)=-2;
c(2)=1;
r=zeros(1,m);
r(1)=-2;
r(2)=1;
% m1=[ones(m,1),-3*ones(m,1),ones(m,1)];
% T= spdiags(a,-1:1,m,m);
T=toeplitz(c,r);%tridiag matrix of -2 and 1
if simu.NoGround
    T(1,1)=-2+tau1;
    T(1,2)=1+tau2;
else
    T(1,1)=-2+sigma1;
    T(1,2)=1+sigma2;
end
T(m,m)=-2+tau1;
T(m,m-1)=1+tau2;

M1=speye(m)+0.5*simu.Rstep*(gamma*T+D)+1/(2*1i*ka)*(gamma*T+D);
M2=speye(m)-0.5*simu.Rstep*(gamma*T+D)+1/(2*1i*ka)*(gamma*T+D);

M1=sparse(M1);
M2=sparse(M2);
%Main process: propagation of the starting field
turbulence=p.Results.Turbulence;
if isempty(turbulence)%simulation without atm turbulence
    tic
    for j=2:n
        %psigrid(:,j)=flipud(M2\M1*flipud(psigrid(:,j-1)));
        psigrid(:,j)=flipud(tridiagonal(M2,M1*flipud(psigrid(:,j-1))));

    end
    t1=toc;

    rgrid=repmat(rvect,m,1);%a grid with the value of the r coordinate at each point
    psi2p_mat=exp(1i*ka*rgrid).*(rgrid).^(-0.5);%a matrix to go from psi to p(pressure); p=psi*exp(i.ka.r)/sqrt(r)
    pgrid=psigrid.*psi2p_mat;%grid of the complex pressure amplitude(see Salomons 2001, eq 2.4)
    pref=2e-5;%20 micro Pascal
    splgrid=10*log10(0.5*abs(pgrid).^2/pref^2);%grid of the sound pressure level;

    result.SPL=splgrid;
    result.P=pgrid;
    result.K=kgrid;

    %Some field that are empty because turbulence is not included
    result.allSPL=[];
    result.allMu=[];
    result.Mu=[];
    result.stdMu=[];
    result.turbulence=[];

else %if turbulence need to be added
    rgrid=repmat(rvect,m,1);%a grid with the value of the r coordinate at each point
    psi2p_mat=exp(1i*ka*rgrid).*(rgrid).^(-0.5);%a matrix to go from psi to p(pressure); p=psi*exp(i.ka.r)/sqrt(r)
    pref=2e-5;%20 micro Pascal
    %mugridtot=zeros(m,n);%addition of each mu grid

    tic
    N=turbulence.NofRealisation;
    allSPL=cell(1,N);%a cell array where all the spl of the different realisation will be stored
    if ~light
        allMu=cell(1,N);%a cell array where all the mu fields of the different realisation will be stored
    else
        allMu=[];
    end
    for k=1:N
        disp(strcat("Turbulence: step ",string(k)," on ",string(N)))
        mugrid=make_mugridPE(turbulence,zvect,rvect);%a random grid of refractive index fluctuation.
        if ~light
            allMu{k}=mugrid;
        end
       % mugridtot=mugrid+mugridtot;
        for j=2:n-1     %propagating the sound field for one random realisation

            mu2=0.5*(mugrid(:,j)+mugrid(:,j+1));
            psigrid(:,j)=flipud(tridiagonal(M2,M1*flipud(psigrid(:,j-1))).*exp(1i*ka*mu2*simu.Rstep));

        end
        if light
            allSPL{k}=single(10*log10(0.5*abs(psigrid.*psi2p_mat).^2/pref^2));%grid of the sound pressure level;
        else
            allSPL{k}=10*log10(0.5*abs(psigrid.*psi2p_mat).^2/pref^2);%grid of the sound pressure level;

        end

        if reso~=0
            allSPL{k}=BlockMean(allSPL{k},reso,reso);
            allMu{k}=BlockMean(allMu{k},reso,reso);
        end
    end
    t1=toc;
    result.allSPL=allSPL;
    result.allMu=allMu;
    if light==0
        allMumat=cat(3,allMu{:});
        result.stdMu=std(allMumat,[],"all");
        result.Mu=sum(allMumat,3)/N;
    else
        result.Mu=0;
        result.stdMu=0;
    end
    allSPLmat=cat(3,allSPL{:});
    clear("allSPL");
    SPL=10*log10((1/N)*sum(10.^(allSPLmat/10),3));%logarithmic average of the SPL correspond to an averaging of the pressure
    result.SPL=SPL;
    
    
    result.P=[];
    result.K=[];
    result.turbulence=turbulence;

end



if reso==0
result.Z=zvect;
result.R=rvect;
else
    N=length(zvect);
    M=length(rvect);
    result.Z=zvect(round(reso/2):reso:N);
    result.R=rvect(round(reso/2):reso:M);
end
result.Ztop=simu.Zabs;
result.Name=p.Results.Name;
result.ComputationTime=t1;
result.Freq_Source=source.freq;

result.shadowzone=[];
%Dealing with the amplitude of the source
%finding the point situated 10cm in front of the source
[~,indr]=min(abs(result.R-0.1));
[~,indz]=min(abs(result.Z-source.Zsource));
indz=length(result.Z)-indz;%the first row of the grid correspond to the maximum altitude
spl0=result.SPL(indz,indr);
%Scaling the entire SPL grid
result.SPL=result.SPL-(spl0-source.Amplitude);
l=length(result.allSPL);
for i=1:l
    result.allSPL{i}=result.allSPL{i}-(spl0-source.Amplitude);
end
%%%
%result.SPL=result.SPL-10*log(4*pi*rgrid.^2);%geometrical attenuation
%%%
result.SPLmin=result.SPL(length(result.Z)-indz,end);
result.source=source;
mic=p.Results.mic;
if ~isempty(mic)
    mic=mic{1,1};
    flag=max(mic.Z)<=simu.Zabs& max(mic.R)<=simu.Xmax;
    if ~flag
        disp("Error: The mic position must be between 0 and Xmax and between 0 and Zabs. The amplitudes at the mic positions were not added.")
    else
        l=length(mic.Z);
        for i=1:l
            z=round(mic.Z(i)/simu.Zstep);
            r=round(mic.R(i)/simu.Rstep);
            spl=flipud(splgrid);
            amp(i)=spl(z,r);  %#ok<AGROW>
        end
        mic.amp=amp;
        result.mic=mic;

    end
end
end


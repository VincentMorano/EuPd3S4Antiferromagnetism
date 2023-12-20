function EXP = EuPd3S4Exp1
%EuPd3S4Exp1EXP Recieve inputs from EuPd3S4Exp1Nuclear.m to output EXP
%structure
%   Send experimental parameters here. Then return experiment structure.
%   ResLib. Usually only fit monochromator and sample mosaic. Nominal BT7
%   monochromator and analyzer mosaics are 30'. Sample should be much less
%   than these given narrow peaks.

%-------------------------   Computation type    -----------------------------------------
EXP.method=0; %Cooper-Nathans approximation by default.
%EXP.method=1; % Popovici approximation.
EXP.moncor=1; % Intensity normalized to flux on monitor by default...
%EXP.moncor=0; % Intensity normalization to white beam flux.

%-------------------------   Monochromator and analyzer    -------------------------------
EXP.mono.tau='PG(002)';     %PG(002), PG(004), GE(111), GE(220),GE(311),BE(002) or PG(110)...
%EXP.mono.tau=3.74650;      %...or any manual input, in reciprocal angstroms
EXP.mono.mosaic=30;         %Minutes of arc, default=25
%EXP.mono.vmosaic=45;       %...For anisotropic mosaic: vertical mosaic, minutes of arc
EXP.ana.tau='PG(002)';      %PG(002)
EXP.ana.mosaic=30;          %Minutes of arc, default=25
%EXP.ana.vmosaic=25;        %...For anisotropic mosaic: vertical mosaic, minutes of arc

%-------------------------   Sample -------------------------------------------------------
EXP.sample.a=6.691070;      %A: angstroms (expt1: 6.69107, expt2: 6.68963, icsd: 6.67650, our xrd at 213 K: 6.67570(7))
EXP.sample.b=6.691070;      %B: angstroms
EXP.sample.c=6.691070;      %C: angstroms
EXP.sample.alpha=90;    %Alpha: degrees of arc
EXP.sample.beta=90;     %Beta: degrees of arc
EXP.sample.gamma=90;    %Gamma: degrees of arc
EXP.sample.mosaic=10;  %Optional sample mosaic: minutes of arc, default=30
%EXP.sample.vmosaic=60; %...For anisotropic mosaic: vertical mosaic, minutes of arc

%-------------------------   Soller and neutron guide collimation    ----------------------
EXP.hcol=[77 25 25 120];        % Horizontal collimation: FWHM minutes of arc, limited to a maximum by FWHM calculation (2/13/23 notes)
EXP.vcol=[300 300 300 300];     % Vertical collimation: FWHM minutes of arc, vertical focusing
                                % If the beam divergence is limited by a neutron guide, 
                                % the parameter value is the negative of the guide's m-value. 
                                % For example, for a 58-Ni guide (m=1.2) the corresponding 
                                % values in EXP.hcol and EXP.vcol is -1.2.

%-------------------------   Fixed neutron energy    --------------------------------------
EXP.efixed=35.0;    %Fixed neutron energy in meV
%EXP.infin=-1;      % Fixed final energy: default value
%EXP.infin=1;       % .. negative for fixed incident energy

%-------------------------   Experimental geometry    -------------------------------------
%EXP.dir1=1;        %  Monochromator scattering direction opposite to sample scattering by default
%EXP.dir1=-1;       %  ...negative if scattering directions are the same.
%EXP.dir2=1;        %  Analyzer scattering direction opposite to sample scattering by default
%EXP.dir2=-1;       %  ...negative if scattering directions are the same.
%EXP.mondir=1;      %  Monochromator scattering angle is positive by default

EXP.orient1=[1 1 0]; %First orienting vector in scattering plane
EXP.orient2=[0 0 1]; %Second orienting vector in scattering plane

%-------------------------   Analyzer reflectivity parameters ----------------------------
% Do not assign these fields if no reflectivity corrections are to be made
%EXP.ana.thickness=0.2;     %Analyzer thickness
%EXP.ana.Q=0.1287;			%Kinematic reflectivity coefficient for PG(002)


%-------------------------   Horizontally focusing analyzer  ----------------------------
% These fields are used only if mode=0
%EXP.horifoc=1;  %Horizontally-focused analyzer
EXP.horifoc=-1; %Flat analyzer
end
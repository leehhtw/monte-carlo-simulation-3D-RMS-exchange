%% Set up and perform Monte Carlo simulation of diffusion in a sphere at the center of FOV
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');

target = fullfile(root,'onesphere');
mkdir(target);

%% Create a sphere at the cneter of FOV

Rsph = 5;           % sphere radius, um
Lvox = 20;          % side length of the FOV, um
Nvox = 200;         % matrix size of the geometry
lvox = Lvox/Nvox;   % side length of each voxel

% Create the geometry, ICS = 1, ECS = 0
[Y,X,Z] = meshgrid(1:Nvox, 1:Nvox, 1:Nvox);
X = X-0.5; Y = Y-0.5; Z = Z-0.5;
BW = (X*lvox - Lvox/2).^2 + (Y*lvox - Lvox/2).^2 +(Z*lvox - Lvox/2).^2 <= Rsph^2;

% In simulations, ICS = 1, ECS = 2
BW = 2-uint8(BW);

% Save the geometry
RMS = rmsobj();
RMS.saveSubstrate(fullfile(target,'geometry.bin'),BW,round(Lvox/Nvox*1e3));

%% Plot geometry
tic;
[F,V] = isosurface(BW,1);
toc;

tic;
TR = triangulation(F,V);
toc;

tic;
figure; trisurf(TR,'edgealpha',0,'facecolor',0.7*[1 1 1]); xlim([0 Nvox]); ylim([0 Nvox]); zlim([0 Nvox]); pbaspect([1 1 1]);
toc;
material dull
camlight

%% Create b-value for signal calculation (not used in the analysis)
bvec = eye(3);      % gradient direction
bval = 1;           % b-value, ms/um2
DEL  = 2;           % diffusion time, ms
del  = 1;           % pulse width, ms
TE = DEL + del;     % echo time = diffusion time + pulse width
btab = btabgen(bvec,bval,[DEL del TE],fullfile(target,'btable.txt'));

%% Perform simulations
perm = 0.1;     % permeability, um/ms
Din = 2;        % ICS intrinsic diffusivity, um2/ms
Dex = 0.5;      % ECS intrinsic diffusivity, um2/ms


% Create the shell script to run simulations
fileID = fopen(fullfile(target,'job.sh'),'w');
rootrms = fullfile(root0,'lib','rms');
% Please go to the directory rootrms and open a terminal to compile the 
% code if you have not done so: sh compile_rms.sh
fprintf(fileID,'#!/bin/bash\n');
filename = 'geometry.bin';
command = sprintf('%s/myrms %s %s/btable.txt %s -time %u -particle %u -space 4 -permeability %.2f -dintra %.2f -dextra %.2f -mspoints 100',...
    rootrms,...
    fullfile(target,filename),...
    target,...
    fullfile(target,'output.bin'),...
    10,...
    1e6,...
    perm,...
    Din,...
    Dex);
fprintf(fileID,command);
fclose(fileID);

% Please go the directory target and run simualtions in a terminal: sh job.sh

%% Read the simulation results
% Read results
RMS = rmsobj(fullfile(target,'output.bin'),fullfile(target,'btable.txt'));

% Read substrate
[BW, vs] = RMS.readSubstrate(fullfile(target,'geometry.bin'));
vs = vs*1e-3;                       % voxel size, um

NParbin = RMS.NParbin;              % number of random walkers in each spherical shell bin
[nt, nbin] = size(NParbin);         % number of time points x number of bins
n = size(BW,1);                     % matrix size, bin width = (n*vs/2) / nbin
R = (0:nbin) * (n*vs/2)/nbin;       % radius of the spherical shells, um
V = 4/3*pi*R.^3; V = diff(V);       % volume fo the spherical shells, um^3
d = NParbin./V;                     % partical density, #/um^3
x = (0.5:nbin-0.5)*(n*vs/2)/nbin;   % position of the bin center, um

%% Calculate the permeability

kappa = zeros(nt,1);        % initialize permeability
Rs = 5;                     % sphere radius, um
[~, IR] = min(abs(x-Rs));   % bin index of the membrane position
nl  = 10;                   % number of data points to fit the density flux over the left side (inside cell)
nle = 3;                    % number of data points (next to the membrane) to be excluded from the fitting (inside cell)
nr  = 10;                   % number of data points to fit the density flux over the right side (outside cell)
nre = 4;                    % number of data points (next to the membrane) to be excluded from the fitting (outside cell)
for i = 1:nt
    di = d(i,:);                    % particle density at time t
    Il = IR-nle-nl+1:IR-nle;        % the bin index over the left side (inside cell)
    Ir = IR+nre     :IR+nre+nr-1;   % the bin index over the right side (outside cell)
    
    dl = di(Il);                    % density over the left side
    xl = x(Il);                     % bin position over the left side
    dr = di(Ir);                    % density over the right side
    xr = x(Ir);                     % bin position over the right side
    
    % Fit the density to a linear model (inside cell)
    Al = [ones(nl,1), xl(:)];
    Xl = Al\dl(:);
    
    % Fit the density to a linear model (outside cell)
    Ar = [ones(nr,1), xr(:)];
    Xr = Ar\dr(:);
    
    % Average the flux inside and outside cell
    flow = (RMS.Din*abs(Xl(2))+RMS.Dex*abs(Xr(2)))/2;
%     flow = RMS.Dex*abs(Xr(2));
    
    % Calculate particle density offset at the membrane
    offset = ([1, Rs]*Xl - [1, Rs]*Xr);
    
    % Permeability, um/ms
    kappa(i) = flow/offset; 
end

% Smooth the permeability over time
kappa = smooth(kappa,10,'sgolay');

% Plot the permeability over time
figure;
plot(RMS.TD, kappa)
h1 = yline(RMS.kappa);           % input value of permeability
h2 = yline(RMS.kappa*1.5);   % permeability scaled by the ratio of (S/V)_(voxelized geometry) to (S/V)_(smooth geometry)
set(h1,'color','k');
set(h2,'color','r');
xlim([0 max(RMS.TD)]);
ylim([0 0.5]);
legend([h1 h2],{'input permeability','scaled permeability'},'fontsize',20,'box','off');
xlabel('diffusion time, ms','fontsize',20);
ylabel('permeability, \mum/ms','fontsize',20);

%%
i = 500;
di = d(i,:);                    % particle density at time t
Il = IR-nle-nl+1:IR-nle;        % the bin index over the left side (inside cell)
Ir = IR+nre     :IR+nre+nr-1;   % the bin index over the right side (outside cell)

dl = di(Il);                    % density over the left side
dr = di(Ir);                    % density over the right side

figure;
hold on;
plot(x,  di, '.');
plot(xl, dl, 'ro');
plot(xr, dr, 'bo');
xlabel('bin position, \mum','fontsize',20);
ylabel('particle density','fontsize',20);
xline(Rs);
xlim([4 6]);
ylim([0 2000]);
box on;


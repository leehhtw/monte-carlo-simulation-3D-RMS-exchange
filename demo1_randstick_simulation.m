%% Set up and perform Monte Carlo simulation of diffusion in randomly positioned, randomly oriented sticks
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');

%% Calculate the mean length of randomly positioned, randomly oriented sticks in a box
N = 2e6;                    % number of sampled sticks
pos = rand(N,3);            % randomly positioned sticks
dir = normr(randn(N,3));    % randomly oriented sticks

% Calculate the length of each stick in the field of view
rs = randstick();
[~,~,L] = rs.stickend(pos,dir);

% Plot the histogram of stick length
figure; hist(L,100);

fprintf('The mean lenght of randomly positioned, randomly oriented sticks is %.4f * side length.\n', mean(L));
% The answer is 0.8965 * side length of the field of view

%% Generate geometry
rs = randstick();

seed = 1;       % seed of random number generation
abar = 5;       % mean distance between beadings, micrometer
astd = 2.5;     % std of distance between beadings, micrometer
lbar = 5;       % bead width, micrometer
rcsa = 0.5;     % mean radius of cross-section, micrometer
cv = [0 0.1 0.2];       % coefficient of variation of radius = std(r)/mean(r)
Lvox = 30;      % length of the field of view, micrometer
Nvox = 600;     % matrix size of the geometry

f = [0.4 0.5 0.6];        % targeted intra-cellular volume fraction

for i = 1:numel(cv)
    cvi = cv(i);
    for j = 1:numel(f)
        fj = f(j);
        % The estimation of number of randomly positioned, randomly oriented fibers
        % 0.8964 corrects the mean length, 1.4 corrects the overlapping
        N = fj*Lvox^3/(pi*rcsa^2*Lvox*0.8965)*1.4;
        N = round(N);

        % Initialize random position and random orientation
        rng(seed);                  % use the same seed for reproducibility
        pos = rand(N,3);            % cylinder position
        dir = normr(randn(N,3));    % cylinder direction

        % Generate the geometry of randomly positioned, randomly oriented fibers
        % 1 = ICS, 0 = ECS
        tic;
        BW = rs.multistick(Nvox,Lvox,pos,dir,seed,abar,astd,lbar,rcsa,cvi);
        toc;

        % Output the actual volume fraction
        fprintf('ICS volume fraction = %.2f, cf. %.2f.\n', nnz(BW==1)/numel(BW), fj);

        % In simulations, 1 = ICS, 2 = ECS
        BW = 2-uint8(BW);

        % Save geometry in the bin file
        RMS = rmsobj();
        filename = sprintf('cv%u_f%u', cvi*100, fj*100);
        mkdir(fullfile(root,filename));
        RMS.saveSubstrate(fullfile(root,filename,'fiber.bin'), BW, round(Lvox/Nvox*1e3));
    end
end

%% Visualize geometry in 3-dimension
cv = [0 0.1 0.2];         % coefficient of variation of radius = std(r)/mean(r)
f = [0.4 0.5 0.6];        % intra-cellular volume fraction
i = 1; j = 1;
proj = sprintf('cv%u_f%u', cv(i)*100, f(j)*100);

% Read geometry
RMS = rmsobj();
[BW, vs] = RMS.readSubstrate(fullfile(root,proj,'fiber.bin'));

% % Undersample for fast visualization (optional)
% BW_undersample = imresize3(BW, 0.5, 'nearest');
BW_undersample = BW;

% Translating to the triangulation takes 20 min if no undersampling
tic;
[F,V] = isosurface(BW_undersample,1);
toc;

tic;
TR = triangulation(F,V);
toc;

% Plot geometry
[n1, n2, n3] = size(BW_undersample);
tic;
figure; trisurf(TR,'edgealpha',0,'facecolor',0.7*[1 1 1]); xlim([0 n1]); ylim([0 n2]); zlim([0 n3]); pbaspect([1 1 1]);
toc;
material dull
camlight
axis off

%% Run simulations
cv = [0 0.1 0.2];         % coefficient of variation of radius = std(r)/mean(r)
f = [0.4 0.5 0.6];        % targeted intra-cellular volume fraction

fileID = fopen(fullfile(root,'job.sh'),'w');
rootrms = fullfile(root0,'lib','rms');
fprintf(fileID,'#!/bin/bash\n');
filename = 'fiber.bin';
for i = 1:numel(cv)
    cvi = cv(i);
    for j = 1:numel(f)
        fj = f(j);
        % Project name
        proj = sprintf('cv%u_f%u', cvi*100, fj*100);
        target = fullfile(root,proj);

        % Read geometry
        RMS = rmsobj();
        [BW, vs] = RMS.readSubstrate(fullfile(root,proj,'fiber.bin'));
        vs = vs*1e-3;                   % side length of each pixel, nm

        % Setup permeability
        fi = nnz(BW==1)/numel(BW);      % ICS volume fraction
        Sx = nnz(diff(BW,1,1))*vs^2;    % surface area in x-direction
        Sy = nnz(diff(BW,1,2))*vs^2;    % surface area in y-direction
        Sz = nnz(diff(BW,1,3))*vs^2;    % surface area in z-direction
        S = Sx+Sy+Sz;                   % surface area, um^2
        V = nnz(BW==1)*vs^3;            % volume, um^3
%         tex = 25;                       % exchange time, ms
%         perm = (1-f)/(tex*S/V);         % permeability, um/ms
        perm = 0.01;                    % permeability, um/ms
        Din = 1;                        % ICS intrinsic diffusivity, um2/ms
        Dex = 2;                        % ECS intrinsic diffusivity, um2/ms
        fprintf('cv%u, f%u, exchange time = %.2f ms.\n',cvi*100, round(fi*100), (1-fi)/perm/(S/V));
        
        % Design b-value for PGSE, not really used in this simulation
        gdir = eye(3);
        bval = 1;
        DEL = 2;
        del = 1;
        TE = DEL + del;
        btab = btabgen(gdir,bval,[DEL del TE],fullfile(target,'btable.txt'));

        % Create shell script for simulations
        command = sprintf('%s/myrms %s %s/btable.txt %s -time %u -particle %u -space 3 -permeability %.5f -dintra %.2f -dextra %.2f -mspoints 50\n',...
            rootrms,...                             % directory to myrms CUDA code
            fullfile(target,filename),...           % input file name of geometry
            target,...                              % input file name of b-table
            fullfile(target,'output.bin'),...       % outoput file name
            100,...                                 % total time, ms
            2e6,...                                 % number of random walker
            perm,...                                % permeability, um/ms
            Din,...                                 % ICS intrinsic diffusivity, um2/ms
            Dex);                                   % ECS intrinsic diffusivity, um2/ms
        fprintf(fileID,command);
    end
end
fclose(fileID);
% Please run the shell script in terminal: sh job.sh

%% Read simulation results
cv = [0 0.1 0.2];          % coefficient of variation of radius = std(r)/mean(r)
f =  [0.4 0.5 0.6];        % targeted intra-cellular volume fraction
MD = zeros(5000,3,3);
MK = zeros(5000,3,3);
tex = zeros(3,3);
fgt = zeros(3,3);
for i = 1:numel(cv)
    cvi = cv(i);
    for j = 1:numel(f)
        fj = f(j);
        % Project name
        proj = sprintf('cv%u_f%u', cvi*100, fj*100);
        target = fullfile(root,proj);
        
        % Read geometry
        RMS = rmsobj();
        [BW, vs] = RMS.readSubstrate(fullfile(target,'fiber.bin'));
        vs = vs*1e-3;                   % side length of each pixel, nm

        % Setup permeability
        fi = nnz(BW==1)/numel(BW);      % actual ICS volume fraction
        Sx = nnz(diff(BW==1,1,1))*vs^2; % surface area in x-direction
        Sy = nnz(diff(BW==1,1,2))*vs^2; % surface area in y-direction
        Sz = nnz(diff(BW==1,1,3))*vs^2; % surface area in z-direction
        S = Sx+Sy+Sz;                   % surface area, um^2
        V = nnz(BW==1)*vs^3;            % volume, um^3
        
        sim = rmsobj(fullfile(target,'output.bin'),fullfile(target,'btable.txt'));
        perm = sim.kappa;               % membrane permeability, um/ms
        fgt(i,j) = fi;                  % actual volume fraction
        tex(i,j) = (1-fi)/perm/(S/V);   % actual exchange time, ms
        
        [Ki,Di] = sim.akc_mom(normr(randn(100,3)));
        MDi = mean(Di,2);
        MKi = mean(Ki,2);
        MD(:,i,j)  = MDi;
        MK(:,i,j)  = MKi;
    end
end
t   = sim.TD;
save(fullfile(root,'simulation_results.mat'),'cv','f','fgt','t','MD','MK','tex');


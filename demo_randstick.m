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
seed = 1;       % seed of random number generation
abar = 5;       % mean distance between beadings, micrometer
astd = 2.5;     % std of distance between beadings, micrometer
lbar = 5;       % bead width, micrometer
rcsa = 0.5;     % mean radius of cross-section, micrometer
cv = 0.2;       % coefficient of variation of radius = std(r)/mean(r)
% cv = 0;
Lvox = 30;      % length of the field of view, micrometer
Nvox = 600;     % matrix size of the geometry

f = 0.5;        % intra-cellular volume fraction

% The estimation of number of randomly positioned, randomly oriented fibers
% 0.8964 corrects the mean length, 1.4 corrects the overlapping
N = f*Lvox^3/(pi*rcsa^2*Lvox*0.8965)*1.4;
N = round(N);

% Initialize random position and random orientation
rng(seed);                  % use the same seed for reproducibility
pos = rand(N,3);            % cylinder position
dir = normr(randn(N,3));    % cylinder direction

% Generate the geometry of randomly positioned, randomly oriented fibers
% 1 = ICS, 0 = ECS
rs = randstick();
tic;
BW = rs.multistick(Nvox,Lvox,pos,dir,seed,abar,astd,lbar,rcsa,cv);
toc;

% Output the actual volume fraction
fprintf('ICS volume fraction = %.2f.\n', nnz(BW==1)/numel(BW));

% In simulations, 1 = ICS, 2 = ECS
BW = 2-uint8(BW);

% Save geometry in the bin file
RMS = rmsobj();
RMS.saveSubstrate(fullfile(root,'cv20','fiber.bin'), BW, round(Lvox/Nvox*1e3));

%% Visualize geometry in 3-dimension

% Read geometry
RMS = rmsobj();
[BW, vs] = RMS.readSubstrate(fullfile(root,'cv20','fiber.bin'));

% Undersample for fast visualization (optional)
BW_undersample = imresize3(BW, 0.5, 'nearest');
% BW_undersample = BW;

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

%% Check the geometry in slices
RMS = rmsobj();
BW00 = RMS.readSubstrate(fullfile(root,'cv00','fiber.bin'));
BW20 = RMS.readSubstrate(fullfile(root,'cv20','fiber.bin'));

figure; 
s_ = 200;
subplot(121); imagesc(BW00(:,:,s_)); pbaspect([1 1 1]); axis off
subplot(122); imagesc(BW20(:,:,s_)); pbaspect([1 1 1]); axis off

%% Run simulations
% Project name
proj = 'cv20';
target = fullfile(root,proj);

% Read geometry
RMS = rmsobj();
[BW, vs] = RMS.readSubstrate(fullfile(root,proj,'fiber.bin'));
vs = vs*1e-3;                   % side length of each pixel, nm

% Setup permeability
f = nnz(BW==1)/numel(BW);       % ICS volume fraction
Sx = nnz(diff(BW,1,1))*vs^2;    % surface area in x-direction
Sy = nnz(diff(BW,1,2))*vs^2;    % surface area in y-direction
Sz = nnz(diff(BW,1,3))*vs^2;    % surface area in z-direction
S = Sx+Sy+Sz;                   % surface area, um^2
V = nnz(BW==1)*vs^3;            % volume, um^3
tex = 25;                       % exchange time, ms
perm = (1-f)/(tex*S/V);         % permeability, um/ms
Din = 1;                        % ICS intrinsic diffusivity, um2/ms
Dex = 2;                        % ECS intrinsic diffusivity, um2/ms

% Design b-value for PGSE, not really used in this simulation
gdir = eye(3);
bval = 1;
DEL = 2;
del = 1;
TE = DEL + del;
btab = btabgen(gdir,bval,[DEL del TE],fullfile(target,'btable.txt'));

% Create shell script for simulations
fileID = fopen(fullfile(target,'job.sh'),'w');
rootrms = fullfile(root0,'lib','rms');
fprintf(fileID,'#!/bin/bash\n');
filename = 'fiber.bin';
command = sprintf('%s/myrms %s %s/btable.txt %s -time %u -particle %u -space 3 -permeability %.5f -dintra %.2f -dextra %.2f -mspoints 50',...
    rootrms,...                             % directory to myrms CUDA code
    fullfile(target,filename),...           % input file name of geometry
    target,...                              % input file name of b-table
    fullfile(target,'output.bin'),...       % outoput file name
    20,...                                  % total time, ms
    1e7,...                                 % number of random walker
    perm,...                                % permeability, um/ms
    Din,...                                 % ICS intrinsic diffusivity, um2/ms
    Dex);                                   % ECS intrinsic diffusivity, um2/ms
fprintf(fileID,command);
fclose(fileID);

% Please run the shell script in terminal: sh job.sh

%% Analysis
projs = {'cv00','cv20'};
MD = [];
MK = [];
for i = 1:numel(projs)
    target = fullfile(root,projs{i});
    sim = rmsobj(fullfile(target,'output.bin'),fullfile(target,'btable.txt'));
    [Ki,Di] = sim.akc_mom(normr(randn(100,3)));
    MDi = mean(Di,2);
    MKi = mean(Ki,2);
    MD  = cat(2,MD,MDi);
    MK  = cat(2,MK,MKi);
end
t   = sim.TD;

tex = 25;   % exchange time, ms
f = 0.5;    % ICS volume fraction

figure('unit','inch','position',[0 0 10 5]);
clear hd hk
lgtxt = {'no beading','beading'};
cmap = colormap('lines');
for i = 1:2
    subplot(121);
    hold on;
    hd(i) = plot(t,MD(:,i),'color',cmap(i,:));
    xlabel('$t$, ms','interpreter','latex','fontsize',20);
    ylabel('MD, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
    box on; grid on;
    pbaspect([1 1 1]);
    xlim([0 max(t)]);
    ylim([0 2]);
    yticks(0:0.4:2);
    title('$t_{\rm ex}$ = 25 ms, $f$ = 50\%','interpreter','latex',...
        'fontsize',20,'FontWeight','normal');
    
    subplot(122);
    hold on;
    hk(i) = plot(t,MK(:,i),'color',cmap(i,:));
    xlabel('$t$, ms','interpreter','latex','fontsize',20);
    ylabel('MK','interpreter','latex','fontsize',20);
    box on; grid on;
    pbaspect([1 1 1]);
    xlim([0 max(t)]);
    ylim([0 1]);
    yticks(0:0.2:1);
    [~, I] = max(MK(:,i));
%     plot(t(I),MK(I,i),'v','color',cmap(i,:));
    title('$t_{\rm ex}$ = 25 ms, $f$ = 50\%','interpreter','latex',...
        'fontsize',20,'FontWeight','normal');
end

subplot(121);
legend(hd,lgtxt,'interpreter','latex','fontsize',20);

subplot(122);
legend(hk,lgtxt,'interpreter','latex','fontsize',20);



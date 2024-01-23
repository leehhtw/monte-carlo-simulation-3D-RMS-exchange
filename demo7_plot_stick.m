%% Set up and perform Monte Carlo simulation of diffusion in randomly positioned, randomly oriented sticks
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');

%% Generate geometry
rs = randstick();

seed = 1;       % seed of random number generation
abar = 5;       % mean distance between beadings, micrometer
astd = 2.5;     % std of distance between beadings, micrometer
lbar = 5;       % bead width, micrometer
rcsa = 0.5;     % mean radius of cross-section, micrometer
cv = [0 0.2];   % coefficient of variation of radius = std(r)/mean(r)
Lvox = 30;      % length of the field of view, micrometer
Nvox = 600;     % matrix size of the geometry

for i = 1:numel(cv)
    cvi = cv(i);
    N = 15;
    % Initialize random position and random orientation
    rng(seed);                  % use the same seed for reproducibility
    pos = [linspace(0.1,0.9,N).' 0.5*ones(N,1) 0.5*ones(N,1)];  % cylinder position
    dir = repmat([0 0 1],N,1);  % cylinder direction

    % Generate the geometry of randomly positioned, randomly oriented fibers
    % 1 = ICS, 0 = ECS
    tic;
    BW = rs.multistick(Nvox,Lvox,pos,dir,seed,abar,astd,lbar,rcsa,cvi);
    toc;

    % In simulations, 1 = ICS, 2 = ECS
    BW = 2-uint8(BW);

    % Save geometry in the bin file
    RMS = rmsobj();
    filename = sprintf('cv%02u_parallel', cvi*100);
    mkdir(fullfile(root,filename));
    RMS.saveSubstrate(fullfile(root,filename,'fiber.bin'), BW, round(Lvox/Nvox*1e3));
end

%% Visualize geometry in 3-dimension
cv = [0 0.2];         % coefficient of variation of radius = std(r)/mean(r)
for i = 1:numel(cv)
    proj = sprintf('cv%02u_parallel', cv(i)*100);

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
    figure('unit','inch','position',[0 0 8 8]);
    tic;
    trisurf(TR,'edgealpha',0,'facecolor',0.7*[1 1 1]); xlim([0 n1]); ylim([0 n2]); zlim([0 n3]);
    pbaspect([1 1 1]);
    toc;
    material dull
    camlight(75,0);
    axis off
    view(90,25);
end


%% Analyze simulation results
clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root0 = fileparts(filePath);
addpath(genpath(fullfile(root0,'lib')));
root = fullfile(root0,'data');

%% Karger model fitting

load(fullfile(root,'simulation_results_fix_tex.mat'));
tex_fit = zeros(2,1);
K0_fit  = zeros(2,1);
for i = 1:2
    for j = 1
        [~, I1] = min(abs(t-50));
        [~, I2] = min(abs(t-100));
        list = I1:I2;
        ti = t(list);
        Ki = MK(list,i,j);
        init = [10, 2];
        X = kargerfit(ti,Ki,init);
        K0_fit(i,j)  = X(1);
        tex_fit(i,j) = X(2);
    end
end

%% Plot exchange time

figure('unit','inch','position',[0 0 5 5]);
hold on;
k = 0;
mk = {'v','o','x'};
cmap = colormap('lines');
clear h lgtxt
for j = 1:1
    for i = 1:2
        k = k+1;
        h(k) = plot(tex(i,j), tex_fit(i,j), mk{i}, 'linewidth', 1, 'color', cmap(j,:));
        lgtxt{k} = sprintf('CV($r$)=%.2f, $f$=%.2f',cv(i), fgt(i,j));
    end
end
xlim([0 30]);
ylim([0 30]);
xticks(0:5:30);
yticks(0:5:30);
pbaspect([1 1 1]);
xlabel('exchange time (theory), ms','fontsize',14);
ylabel(['exchange time (K' char(228) 'rger model), ms'],'fontsize',14);
hr = refline(1,0); set(hr,'color',[0.5 0.5 0.5]);
box on; grid on;
legend(h,lgtxt,'interpreter','latex','fontsize',12,'box','off','location','northwest');

%% Plot figure

figure('unit','inch','position',[0 0 8 4]);
clear hd hk lgtxt
lgtxt = {'no beading','beading'};
cmap = colormap('lines');
for i = 1:numel(cv)
    cvi = cv(i);
    subplot(1,2,1);
    hold on;
    
    for j = 1:numel(f)
        fj = fgt(i,j);
        MDi = MD(:,i,j);
        hd(i) = plot(t, MDi, 'linewidth', 2, 'color', cmap(i,:));
%         lgtxt{i} = sprintf('CV($r$)=%.2f',cvi);
        xlabel('$t$, ms','interpreter','latex','fontsize',20);
        ylabel('MD, $\mu$m$^2$/ms','interpreter','latex','fontsize',20);
        box on; grid on;
        pbaspect([1 1 1]);
        xlim([0 max(t)]);
        ylim([0 2]);
        yticks(0:0.4:2);
    end
    legend(hd,lgtxt,'interpreter','latex','fontsize',20,'box','off');
%     title(sprintf('CV($r$)=%.2f',cvi),'interpreter','latex','fontsize',20);
    
    subplot(1,2,2);
    hold on;
    for j = 1:numel(f)
        fj = fgt(i,j);
        MKi = MK(:,i,j);
        
        % theory
        t0 = t/tex_fit(i,j);
        Kt = K0_fit(i,j)*2./t0.*(1-1./t0.*(1-exp(-t0)));
        
        hk(i) = plot(t, MKi, 'linewidth', 2, 'color', cmap(i,:));
        plot(t, Kt, '--', 'linewidth', 2, 'color', cmap(i,:));
        xlabel('$t$, ms','interpreter','latex','fontsize',20);
        ylabel('MK','interpreter','latex','fontsize',20);
        box on; grid on;
        pbaspect([1 1 1]);
        xlim([0 max(t)]);
        box on; grid on;
        pbaspect([1 1 1]);
        xlim([0 max(t)]);
        ylim([0 1]);
        yticks(0:0.2:1);
    end
    legend(hk,lgtxt,'interpreter','latex','fontsize',20,'box','off');
%     title(sprintf('CV($r$)=%.2f',cvi),'interpreter','latex','fontsize',20);
    
end






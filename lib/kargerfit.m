function X = kargerfit(t,K,init)

karger = @(t,tex) 2./(t./tex).*(1-1./(t./tex).*(1-exp(-t./tex)));
options = fitoptions('Method','NonlinearLeastSquares',...
    'Algorithm','Trust-Region',...
    'Upper',[Inf Inf],'Lower',[0 0],'StartPoint',init);
ft = fittype(@(K0,tex,t) K0*karger(t,tex),...
    'independent','t');
[kfit,gof] = fit(t,K,ft,options);
ci = confint(kfit);
X(1) = kfit.K0;
X(2) = kfit.tex;
X(3) = range(ci(:,1))/4;
X(4) = range(ci(:,2))/4;
X(5) = gof.rsquare;

end
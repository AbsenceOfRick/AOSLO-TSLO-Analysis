function coefEsts=GammaFit(xdata,ydata)

modelFun =  @(p,x) p(3) .* (x ./ p(1)).^(p(2)-1) .* exp(-(x ./ p(1)).^p(2))+p(4);
startingVals = [10 2 5 0];
coefEsts = lsqcurvefit(modelFun,startingVals,xdata',ydata');

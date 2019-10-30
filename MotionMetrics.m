%This function will perform various analyses on the eye motion to pull out
%metrics such as drift span, saccade speed, etc. This is a flexible
%function where new analyses can be constantly built in. Any additional
%outputs should be formatted into the EM structure. 

%Inputs
%X Stream
%Y Stream
%Sampling Rate
%Drift Start
%Drift End
%Saccade Start
%Saccade End
%Blink Start
%Blink End

%Outputs
%EM structure containing analysis results

%Note: This function is still very much a work in progress.

%Norick Bowers, Spring 2017. 

function [EM] = MotionMetrics(xx,yy,SampRate,DriftS,DriftE,SaccS,SaccE,DropS,DropE,TCorrX,TCorrY)

%%%%%%%%%%%%%%%%% Basic Info %%%%%%%%%%%%%%%%%
EM.Total.xx = xx; EM.Total.yy = yy;
EM.Total.Dur = length(xx) * 1000/SampRate;

%Start/End times of events
EM.Drift.Start = DriftS; EM.Drift.End = DriftE;
EM.Drop.Start = DropS; EM.Drop.End = DropE;
EM.Sacc.Start = SaccS; EM.Sacc.End = SaccE;

%%%%%%%%%%%%%%%%%% Drifts %%%%%%%%%%%%%%%%%%
for aa = 1:length(DriftS)
    
    try
        xxtmp = xx(DriftS(aa):DriftE(aa));
        yytmp = yy(DriftS(aa):DriftE(aa));
    catch
        xxtmp = [];
        yytmp = [];
    end
    
    if ~isempty(xxtmp)
    
    %Error where 1 NaN value is in front or end, eliminate.
    xxtmp(find(isnan(xxtmp))) = [];
    yytmp(find(isnan(yytmp))) = [];
    
    %Amplitude
    EM.Drift.Amp(aa) = round(sqrt( ((xxtmp(1)-xxtmp(end))^2) + ((yytmp(1)-yytmp(end))^2) ),4);
    
    %Speed
    EM.Drift.Spd{aa} =  (abs(diff(sqrt( (xxtmp).^2 + (yytmp.^2) ) * SampRate)));
    EM.Drift.Spd_Mean(aa) = mean(EM.Drift.Spd{aa});
    
    %Span
    mx = nanmean(xxtmp);
    my = nanmean(yytmp);
    
    EM.Drift.Span(aa) = max(sqrt( ((mx - xxtmp).^2) + ((my-yytmp).^2) ));
    
    %Duration
    EM.Drift.Dur(aa) = round(length(xxtmp) * (1000/SampRate)); 
    
    %torsion? onset/offset? goodness of fit? adjusted for framerate?
    
    %Positional Heatmap
    EM.Drift.XPos{aa} = (xxtmp - xxtmp(1))';
    EM.Drift.YPos{aa} = (yytmp - yytmp(1))';
    
    else
        EM.Drift.XPos{aa}= NaN;
        EM.Drift.YPos{aa} = NaN;
        EM.Drift.Dur(aa) = NaN;
        EM.Drift.Span(aa) = NaN;
        EM.Drift.Spd{aa} = NaN;
        EM.Drift.Spd_Mean(aa) = NaN;
        EM.Drift.Amp(aa) = NaN;
    end
    
end

%Heatmap
% Bins = 2;
% maxc = .25;
% Xtmp = cell2mat(EM.Drift.XPos);
% Ytmp = cell2mat(EM.Drift.YPos);
% Xtmp(find(isnan(Xtmp))) = [];
% Ytmp(find(isnan(Ytmp))) = [];
% ndhist(Xtmp,Ytmp,'bins',Bins,'filt',6);
% axis square
% colormap hot


%BCEA & ISOA

%Total X and Y traces
XTMP = xx; XTMP(find(isnan(XTMP))) = []; %remove NaN's X
YTMP = yy; YTMP(find(isnan(YTMP))) = []; %remove NaN's y

[isoa, bcea, PRL, PRL2, density, xGrid, yGrid, fh] = ...
    ComputeFixationStability(XTMP, YTMP, 0.68, 0);
EM.Stability.isoa = isoa;
EM.Stability.bcea = bcea;
EM.Stability.PRL = PRL;
EM.Stability.PRL2 = PRL2;
EM.Stability.density = density;
EM.Stability.xGrid = xGrid;
EM.Stability.yGrid = yGrid;
EM.Stability.fh = fh;

%Create vector with only real drift positions (not stitched together)
tmpX = []; tmpY = [];
for aa = 1:length(DriftS)
    tmpX = [tmpX ; xx(DriftS(aa):DriftE(aa))-xx(DriftS(aa))];
    tmpY = [tmpY ; yy(DriftS(aa):DriftE(aa))-yy(DriftS(aa))];
end
tmpX(find(isnan(tmpX))) = [];
tmpY(find(isnan(tmpY))) = [];

%BCEA Drift Only
[isoa, bcea, PRL, PRL2, density, xGrid, yGrid, fh] = ...
    ComputeFixationStability(tmpX, tmpY, 0.68, 0);

EM.Stability.DriftOnly.isoa = isoa;
EM.Stability.DriftOnly.bcea = bcea;
EM.Stability.DriftOnly.PRL = PRL;
EM.Stability.DriftOnly.PRL2 = PRL2;
EM.Stability.DriftOnly.density = density;
EM.Stability.DriftOnly.xGrid = xGrid;
EM.Stability.DriftOnly.yGrid = yGrid;
EM.Stability.DriftOnly.fh = fh;


%%%%%%%%%%%%%% Saccades %%%%%%%%%%%%%%%%%

%Saccade Rate
EM.Sacc.SaccRate = length(SaccS) / ((length(xx)*(1000/SampRate))/1000);

for aa = 1:length(SaccS)
    
    %Amplitude
    x_s = xx(SaccS(aa)); x_e = xx(SaccE(aa));
    y_s = yy(SaccS(aa)); y_e = yy(SaccE(aa));
    EM.Sacc.Amp(aa) =  sqrt( (x_s-x_e)^2 + (y_s-y_e)^2 );
    
    %Direction
    Deg = rad2deg(atan2((y_e-y_s),(x_e-x_s))); %saccade direction 180 to -180
    EM.Sacc.Angle(aa) = Deg;
    
    %Speed
    EM.Sacc.Spd{aa} = abs(diff(sqrt( (xx(SaccS(aa):SaccE(aa)).^2 + yy(SaccS(aa):SaccE(aa)).^2) ) .* SampRate));
    
    EM.Sacc.Spd_Mean = mean(EM.Sacc.Spd{aa});
    
    %Duration
    EM.Sacc.Dur(aa) = length(xx(SaccS(aa):SaccE(aa))) * 1000/SampRate;
    
end

% Software by Norick Bowers and Josselin Gautier
% School of Optometry, UC Berkeley, California, USA

clc; clear all; close all;
%% Initialize
PPD_H = 560;
PPD_V = 560;

ManualWin = 0.25; %Window for manual checking function (as a percentage of the total trace)
FiltWindow = 51; %Window for loess filter (Lower value = less smoothing/more false positives for saccades)

%Turn on/off to enable/disable certain functions
ShowSep = 1; %Plot Saccade,Drift,Blink seperation
Manual_Check = 1; %Manually check abnormal drift traces
AnalyzeMetrics = 1; %Analyze metrics of eye motion and generate plots
CorrectTorsion = 0; %Correct torsion
Analyze_Fourier = 0; %Analyze spectral properties of drifts
Save_On = 0; %Save workspace in directory
Load_Demarcation = 0; %Load eye trace demarcation from processed file

if ispc
    [Curr_File Directory] = uigetfile('./',...
        'select eye trace as MAT file');
else
    [Curr_File Directory] = uigetfile('./',...
        'select eye trace as MAT file');
end

%load(sprintf('%s/%s.mat',Directory,Curr_File)); %Load
load(sprintf('%s%s',Directory,Curr_File));

if Load_Demarcation %Load previous demarcations
    File_Name = sprintf('%s/%s.mat',Directory,sprintf('%s_Processed',Curr_File));
    load(File_Name,'DropS','DropE','DriftS','DriftE','SaccS','SaccE','RejectedS','RejectedE', ...
        'XFilt','YFilt','IVals','xx','yy','xx_original','yy_original');
end


%Initialize variables
PxArcminH = ( (512/PPD_H) * 60 ) / 512; %Pixel to Arcmin Conversion (arcmin in 1 pixel)
PxArcminV = ( (512/PPD_V) * 60 ) / 512;
ChooseSamp = 1:length(frameshifts_strips_spline);  %Which samples to choose (must be from beginning, default to all)
FrameHeight = round(framewidth * PxArcminH);
FrameRate = videoframerate; %Framerate
SampRate = samplerate;  %Samplerate (Equal to FrameRate * SPF)
SPF = SampRate/FrameRate; %Samples per frame
if rem(max(ChooseSamp)/SPF,1) ~= 0 %Make sure samples and frames align (only for choosesamp variable)
    ChooseSamp = [1: ((floor(max(ChooseSamp)/SPF) +1)*SPF)];
end
SampFrame = 1:SPF:length(ChooseSamp) ; %Starts of each frame
UsedFrame = floor(ChooseSamp(SampFrame)/SPF); %Which frames of the video are in samples chosen
BlinkFrames = blinkframes(find(blinkframes>min(min(UsedFrame)) & blinkframes<max(max(UsedFrame)))); %Get dropped frames
Samp_per_Sec = min(min(diff(timeaxis_secs))); %Time seperation between samples (sec)
DropE = []; DropS = []; SaccS = 1; SaccE = [];


%% Arrange

%X and Y positions of entire movie
if ~Load_Demarcation
    xx = frameshifts_strips_spline(:,1);
    yy = frameshifts_strips_spline(:,2);
    IVals = zeros(length(maxvals_strips),1); %interpolated values
    IVals(find(maxvals_strips==0)) = 1;
    
    %Chosen samples
    xx = xx(ChooseSamp);
    yy = yy(ChooseSamp);
    
    %Convert to arcmin
    xx = xx.*PxArcminH;
    yy = yy.*PxArcminV;
end

%Rebuild and obtain start/end of dropped frames
if ~Load_Demarcation
    [xx,yy,DropS,DropE,xx_Original,yy_Original,IVals] = SepBlink(xx,yy,Samp_per_Sec,timeaxis_secs,SPF,IVals);
end

TimeAx = linspace(0,timeaxis_secs(end),length(xx));

%Pull out repeated frames
[RepeatedFrames,RepeatedSamples] = SepRepeat(xx,yy,SPF);

%Get start and end times for saccades
if ~Load_Demarcation
    [SaccS,SaccE,XFilt,YFilt,xx,yy] = SepSaccade(xx,yy,SampRate,FiltWindow,SPF);
end

%Get start and end times for drifts
if ~Load_Demarcation
    if ~isempty(SaccS)
        [DriftS,DriftE] = SepDrift(SaccS,SaccE,DropS,DropE,xx);
    else
        DriftS = 1; DriftE = length(xx);
        warning(sprintf('\nSaccades not separated, using entire trace\n'));
    end
end

%Find incorrectly labeled blinks/saccades
[SaccS,SaccE,autoRejS,autoRejE,DropS,DropE] = IncorrectBlinks(SaccS,SaccE,DropS,DropE,xx,yy);

%Manually Check for missed saccades
if ~Load_Demarcation
    if Manual_Check
        [SaccS,SaccE,RejectedS,RejectedE,DriftS,DriftE] = ManualCheck(xx,yy,SPF,SaccS,SaccE,DriftS,DriftE,DropS,DropE, autoRejS, autoRejE,ManualWin,XFilt,YFilt);%
    else
        RejectedS = [];
        RejectedE = [];
    end
    
    %drop negative values and NaN values.
    Dstmp = DriftS; Detmp = DriftE;
    Dstmp(find((DriftE-DriftS)<=1)) = [];
    Detmp(find((DriftE-DriftS)<=1)) = [];
    tmp = find(isnan(Dstmp) | isnan(Detmp) );
    Dstmp(tmp) = [];
    Detmp(tmp) = [];
    DriftS = Dstmp;
    DriftE = Detmp;
    
    Sstmp = SaccS; Setmp = SaccE;
    Sstmp(find((SaccE-SaccS)<=1)) = [];
    Setmp(find((SaccE-SaccS)<=1)) = [];
    tmp = find(isnan(Sstmp) | isnan(Setmp) );
    Sstmp(tmp) = [];
    Setmp(tmp) = [];
    SaccS = Sstmp;
    SaccE = Setmp;
end

%Correct Torsional Sawtooth
if CorrectTorsion
    [DriftSegX,DriftSegY] = ReferenceFrameCorrection(xx,yy,DriftS,DriftE,SPF); %Get stitched together drifts (ignore reference correction for now)
    [TCorrX,TCorrY] =  TorsionCorrection(DriftSegX,DriftSegY,SPF);
else
    DriftSegX = xx'; DriftSegY = yy';
    TCorrX = xx'; TCorrY = yy';
end

%Eye Motion Analysis
if AnalyzeMetrics
    [EM] = MotionMetrics(xx,yy,SampRate,DriftS,DriftE,SaccS,SaccE,DropS,DropE,TCorrX,TCorrY);
end

%Run Spectral Analysis
if  Analyze_Fourier
    PmtmWin = 1000;
    THB = 2.5;
    Overlap = 50;
    
    [F PxxM PyyM Lambda DPS_Seq] = MultiTaper_PowSpec(DriftSegX,DriftSegY,PmtmWin,SampRate,THB,Overlap);
    
    [F PxxMT PyyMT Lambda DPS_Seq] = MultiTaper_PowSpec(TCorrX,TCorrY,PmtmWin,SampRate,THB,Overlap);
    
end


%Save Entire Workspace
if Save_On == 1
    save(sprintf('%s/%s.mat',Directory,sprintf('%s_Processed',Curr_File)));
    fprintf(sprintf('\n\nSaving...\n\n'))
end

%% Figures
if ShowSep == 1
    
    %Plot X and Y motion
    figure;
    plot(xx-xx(1),'b'); hold on;
    plot(yy-yy(1),'r'); hold on;
    legend('X','Y','Location','Best');
    xlabel('Frames');
    ylabel('Position (arcmin)');
    ylim( [min([min(yy)-10,min(xx)-10]), max(max(yy)+10,max(xx)+10) ]);
    Yh = get(gca,'ylim');
    
    %Highlight Saccades w/ magenta
    for aa = 1:length(SaccS)
        H = fill([SaccS(aa) SaccE(aa) SaccE(aa) SaccS(aa)],[max(Yh) max(Yh) min(Yh) min(Yh)],'m');
        set(H, 'FaceAlpha', 0.4,'EdgeAlpha',0); hold on;
    end
    
    
    %Highlight Blinks/Dropped Frames w/ black
    for aa = 1:length(DropS)
        H = fill([DropS(aa) DropE(aa) DropE(aa) DropS(aa)],[max(Yh) max(Yh) min(Yh) min(Yh)],'k');
        set(H, 'FaceAlpha', 0.4,'EdgeAlpha',0); hold on;
    end
    
    %Mark everything else as drift
    for aa = 1:length(DriftS)
        H = fill([DriftS(aa) DriftE(aa) DriftE(aa) DriftS(aa)],[max(Yh) max(Yh) min(Yh) min(Yh)],'c');
        set(H, 'FaceAlpha', 0.2,'EdgeAlpha',0); hold on;
    end
    
    
    if exist('RejectedS','var')%if Manual_Check
        
        for aa = 1:length(RejectedS)
            H = fill([RejectedS(aa) RejectedE(aa) RejectedE(aa) RejectedS(aa)],...
                [max(Yh) max(Yh) min(Yh) min(Yh)],'k');
            set(H, 'FaceAlpha', 0.6,'EdgeAlpha',0); hold on;
        end
    end
    
    if exist('autoRejS','var')%if Manual_Check
        
        for aa = 1:length(autoRejS)
            H = fill([autoRejS(aa) autoRejE(aa) autoRejE(aa) autoRejS(aa)],...
                [max(Yh) max(Yh) min(Yh) min(Yh)],'k');
            set(H, 'FaceAlpha', 0.6,'EdgeAlpha',0); hold on;
        end
    end
    
    IPercent = length(find(maxvals_strips==0))/length(maxvals_strips);
    title(sprintf('Entire Trace, %dHz, %.2f Interpolated, %d Frames',SampRate,IPercent,[length(xx)/SPF]));
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    xlim([1 length(xx)]);
    set(gca,'xgrid','on','XTick',[0:SPF:length(xx)],'XTickLabel', ...
        0:max([0:SPF:length(xx)]),'FontSize',10);
    legend('off');
    
    %Mark interpolated values
    hold on;
    scatter(find(IVals==1),xx(find(IVals==1))-xx(1),'k.')
    hold on;
    scatter(find(IVals==1),yy(find(IVals==1))-yy(1),'k.')
    
end
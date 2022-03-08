% software by Norick Bowers and Josselin Gautier
% School of Optometry, UC Berkeley, California, USA
% -------------------------------------------------

clc; clear all; close all;
%% Initialize

%Adjust parameters below by hand as desired before running

if ispc
    [Curr_File Directory] = uigetfile('G:\My Drive\PRL_project\aoslo_data\2nd_Exp\20109R\*.mat',...
        'select eye trace as MAT file');    %2nd_Exp\*.mat',...
else
    [Curr_File Directory] = uigetfile(['/Volumes/GoogleDrive/My Drive/PRL_project/aoslo_data/',...
        '2nd_Exp/20196L/*.mat'],...%20196L/10_4_2019_18_28_57/*.mat',...
        'select eye trace as MAT file');
end

%Px Arcmin Calculation:  512/PPD = FieldSize(Deg).  FieldSize(Deg)*60 = FieldSize(Arc). FieldSize(Arc)/512 = PxArcmin.
%PPD_H = 556; %556 % Pixels per degree X
%PPD_V = 554; %564; %Pixels per degree Y
% Request the exact calibrated Pixel Per Degree %JGadd
prompt = {'Enter Horizontal Pixel per Degree:','Vertical Pixel per Degree:'};
titleTemp = 'Input';%JGc: title is already a function, do not named it as a var
dims = [1 35];
definput = {'558','554'};
answer = inputdlg(prompt,titleTemp,dims,definput);
% Use str2num to convert the cell array to numerals.
PPD_H = str2num(answer{1});
PPD_V = str2num(answer{2});

% -------------------------------------------------
%
PxArcminH = ( (512/PPD_H) * 60 ) / 512; %Pixel to Arcmin Conversion (arcmin in 1 pixel)
PxArcminV = ( (512/PPD_V) * 60 ) / 512;
ManualWin = 0.5; %Window for manual checking function (as a percentage of the total trace)

FiltWindow = 51; %Window for loess filter (Lower value = less smoothing/more false positives for saccades)

%Turn on/off to enable/disable certain functions
ReVAS = 1; %If using ReVAS output, set to 1, if Stevenson set to 0
ShowSep = 1; %Plot Saccade,Drift,Blink seperation
Manual_Check = 1; %Manually check abnormal drift traces
CorrectTorsion = 0; %Correct torsion
Analyze_Fourier = 0; %Analyze spectral properties of drifts (Unavailable, need to update to pmtm method from Welsch)
AnalyzeMetrics = 1; %Analyze metrics of eye motion and generate plots
Save_On = 1; %Save workspace in directory
Load_Demarcation = 0; %Load eye trace demarcation from processed file

%load(sprintf('%s/%s.mat',Directory,Curr_File)); %Load
load(sprintf('%s%s',Directory,Curr_File));
EM=[];% for saving, even if EM not written later
if Load_Demarcation %Load previous demarcations
    File_Name = sprintf('%s/%s.mat',Directory,sprintf('%s_Processed',Curr_File));
    load(File_Name,'DropS','DropE','DriftS','DriftE','SaccS','SaccE','RejectedS','RejectedE', ...
        'XFilt','YFilt','IVals','xx','yy','xx_original','yy_original');
end


%Initialize variables
if ReVAS %ReVAS output
    FrameWidth = 512;
    FrameHeight = 512;
    SampRate = params.samplingRate;
    SPF = FrameHeight/params.stripHeight;
    ChooseSamp = 1:length(position(:,1));
    SampFrame = 1:SPF:length(ChooseSamp) ; %Starts of each frame
    BlinkFrames = params.badFrames;
    Samp_per_Sec = min(min(diff(timeSec))); %JGadd Time seperation between samples (sec)
else %Stevenson output
    
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
end

%% Arrange

%X and Y positions of entire movie
if ~Load_Demarcation
    
    if ReVAS %ReVAS output
        xx = position(:,1);
        yy = position(:,2);
    else %Stevenson output
        xx = frameshifts_strips_spline(:,1);
        yy = frameshifts_strips_spline(:,2);
    end
    
    %Chosen samples
    xx = xx(ChooseSamp);
    yy = yy(ChooseSamp);
    
    %Convert to arcmin
    xx = xx.*PxArcminH;
    yy = yy.*PxArcminV;
    
    xx_UnInt = xx;
    yy_UnInt = yy;
end


%Interpolation
if ~Load_Demarcation
    
    %%%% Entire Trace
    if ReVAS %ReVAS output
        
        %JG mod save nan, as hampel tend to filter a big part of them
        xxnan=isnan(xx);
        yynan=isnan(yy);
        %JGmod : detect outliers BEFORE interpolation
        [xx,Ix] = hampel(xx,30,2);%3,1);
        [yy,Iy] = hampel(yy,30,2);%3,1);
        xx(xxnan)=nan;
        yy(yynan)=nan;
        
        %X Trace interpolation
        SI = 1:length(xx); %sample index
        Samps = xx; %Samples
        Query = find(isnan(xx)); %NaN Values
        ISamps = pchip(SI,Samps,Query); %Do interpolation
        xx(Query) = ISamps; %Replace Values
        
        %Y Trace interpolation
        SI = 1:length(yy); %sample index
        Samps = yy; %Samples
        Query = find(isnan(yy)); %NaN Values
        ISamps = pchip(SI,Samps,Query); %Do interpolation
        yy(Query) = ISamps; %Replace Values
        
        IVals = zeros(1,length(xx));
        IVals(Query) = 1;
        
        %[xx,Ix] = hampel(xx,5,5);
        %[yy,Iy] = hampel(yy,5,5);
        IVals(find(Ix)==1) = 2;
        IVals(find(Iy)==1) = 2;
        IVals=IVals';
        
        %     %%%% Per Frame
        %     for aa = 1:length(SampFrame)-1
        %
        %         IVals(SampFrame(aa):SampFrame(aa+1)-1) = 0;
        %
        %         if BlinkFrames(aa) ~= 1
        %
        %             %Sample for this frame
        %             ThisX = xx(SampFrame(aa):SampFrame(aa+1)-1);
        %             ThisY = yy(SampFrame(aa):SampFrame(aa+1)-1);
        %
        %             if length(find(~isnan(ThisX))) > 2
        %
        %                 %PChip for X
        %                 SI = 1:length(ThisX); %sample index
        %                 Samps = ThisX; %Samples
        %                 Query = find(isnan(ThisX)); %NaN Values
        %                 ISamps = pchip(SI,Samps,Query); %Do interpolation
        %                 ThisX(Query) = ISamps; %Replace Values
        %
        %                 %PChip for Y
        %                 SI = 1:length(ThisY); %sample index
        %                 Samps = ThisY; %Samples
        %                 Query = find(isnan(ThisY)); %NaN Values
        %                 ISamps = pchip(SI,Samps,Query); %Do interpolation
        %                 ThisY(Query) = ISamps; %Replace Values
        %
        %                 %Run hampel filter to remove outliers
        %                 ThisX = hampel(ThisX,3);
        %                 ThisY = hampel(ThisY,3);
        %
        %                 %Replace NaN's with interpolated values
        %                 xx(SampFrame(aa):SampFrame(aa+1)-1) = ThisX;
        %                 yy(SampFrame(aa):SampFrame(aa+1)-1) = ThisY;
        %
        %                 %Turn off Warning
        %                 w = warning('query','last');
        %                 warning('off',w.identifier)
        %
        %                 %Keep track of which samples were interpolated
        %                 ThisInterp = zeros(1,length(ThisX));
        %                 ThisInterp(Query) = 1;
        %                 IVals(SampFrame(aa):SampFrame(aa+1)-1) = ThisInterp;
        %
        %             else
        %                 BlinkFrames(aa) = 1;
        %                 xx(SampFrame(aa):SampFrame(aa+1)-1) = NaN;
        %                 yy(SampFrame(aa):SampFrame(aa+1)-1) = NaN;
        %             end
        %
        %         end
        %     end
        %
        % end
        
        %Mark Bad Frames as NaN;
        Blinks = find(BlinkFrames==1); %Frames w/ Blinks
        if ~isempty(Blinks)
            NonConsec = find(diff(Blinks)~=1); %Non consecutive blinks
            Starts = [1 NonConsec+1];
            Ends = [NonConsec length(Blinks)];
            DropS = Blinks(Starts)*SPF;
            DropE = Blinks(Ends)*SPF;
            DropS = DropS - 31; %This technique scoots this forward 1 frame?
        else
            DropS = [];
            DropE = [];
        end
        
        %Eliminate
        for aa = 1:length(DropS)
            IVals(DropS(aa):DropE(aa)) = 0;
            xx(DropS(aa):DropE(aa)) = NaN;
            yy(DropS(aa):DropE(aa)) = NaN;
        end
        
        TimeAx = linspace(0,timeSec(end),length(xx));%JGadd
        
    else %Stevenson output
        
        IVals = zeros(length(maxvals_strips),1); %interpolated values
        IVals(find(maxvals_strips==0)) = 1;
        
        %Rebuild and obtain start/end of dropped frames
        if ~Load_Demarcation
            [xx,yy,DropS,DropE,xx_Original,yy_Original,IVals] = SepBlink(xx,yy,Samp_per_Sec,timeaxis_secs,SPF,IVals);
        end
        
        TimeAx = linspace(0,timeaxis_secs(end),length(xx));
    end
end%additional 'end' here? JG

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
        [SaccS,SaccE,RejectedS,RejectedE,DriftS,DriftE] = ManualCheckPlus(xx,yy,SPF,SaccS,SaccE,DriftS,DriftE,DropS,DropE, autoRejS, autoRejE,ManualWin,XFilt,YFilt);
    else
        RejectedS = [];
        RejectedE = [];
    end
    
    if ~iscolumn(RejectedS)%if line vector instead of column vector % -------------------------------------JG to commit
        RejectedS=RejectedS';
        RejectedE=RejectedE';
    end
    RejectedS = sort([RejectedS;autoRejS]);
    RejectedE = sort([RejectedE;autoRejE]);
    
    %%
    %drop negative values and NaN values (some sort of bug?).
    Dstmp = DriftS; Detmp = DriftE;
    if length(DriftS)~=length(DriftE)% -----------JG to commit
        fprintf('\nError after ManualCheckPlus, different sizes of vector\n');
    end
    
    
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
    
    [F PxxMT PyyMT Lambda DPS_Seq] = MultiTaper_PowSpec(TCorrX,TCorrY,PmtmWin,SampRate,THB,Overlap)
    
    %Temporary figure to compare torsion correction for thesis proposal.
    %Add below or delete later.
    
    %     figure;
    %     plot(log10(F),log10(PxxM),'Color','b','LineWidth',1.5);
    %     xlim([0 2]);
    %     hold on;
    %     plot(log10(F),log10(PxxMT),'Color','r','LineWidth',1.5);
    %     title('X Trace, Corrected vs Uncorrected');
    %     xlabel('Frequency (Hz)');
    %     ylabel('Log Amplitude (arcmin)');
    %     title('Horizontal Amplitude');
    %     ylim([-2.5 0.5]);
    %     xlim([0 2]);
    %     set(gca,'XTick',[0 1 2],'XTickLabel',{'0','10','100'});
    %     set(gca,'ytick',[-2:1:1])
    %     axis square
    %     Yh = get(gca,'YLim');
    %     Xh = gca;
    %     Xh.XAxis.MinorTick = 'on';
    %     Xh.XAxis.MinorTickValues = double([log10(1:9) log10(10:10:90)]);
    %     legend('Corrected','Uncorrected','Location','Best');
    %     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
end


%Save Entire Workspace
if Save_On == 1
    %{
    save(sprintf('%s/%s.mat',Directory,sprintf('%s_Processed',Curr_File)));
    fprintf(sprintf('\n\nSaving...\n\n'))
    %}
    %JG or save just usefull variables for seek of space when streaming
    %%{
    
    save(sprintf('%s/%s.mat',Directory,sprintf('%s_Processed',Curr_File)),...
        'PPD_H','PPD_V','PxArcminH','PxArcminV',...
        'TimeAx','xx','yy','XFilt','YFilt',...
        'SaccS','SaccE',...
        'DriftS','DriftE','EM','IVals',...
        'autoRejS','autoRejE','RejectedS','RejectedE','DropS','DropE');
    fprintf(['\n file saved:\n']);
    fprintf([Curr_File '_Processed.mat\n']);
    fprintf(['correctly saved in \n' Directory '\n']);
    %}
else
    fprintf('Saving was set to OFF, mat file not saved');
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
    Limits = [prctile(xx,95),prctile(yy,95)];
    ylim([-max(Limits)-10,max(Limits)+10]);
    %     ylim( [min([min(yy)-10,min(xx)-10]), max(max(yy)+10,max(xx)+10) ]);
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
    
    
    if exist('RejectedS')%if Manual_Check
        
        for aa = 1:length(RejectedS)
            H = fill([RejectedS(aa) RejectedE(aa) RejectedE(aa) RejectedS(aa)],...
                [max(Yh) max(Yh) min(Yh) min(Yh)],'k');
            set(H, 'FaceAlpha', 0.6,'EdgeAlpha',0); hold on;
        end
    end
    
    IPercent = length(find(IVals==1))/length(IVals);
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


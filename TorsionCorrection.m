%This function takes the stitched together drift segments and calculates
%the slope of the sawtooth arising at the framerate of the system due to
%torsion. First the real motion is estimate for each slope by taking a
%linear fit of the current frame and the two surrounding frames. Then the
%remaining slope is calculated and removed from the trace

% Input 1: Stitched Together Drift Segment (X trace)
% Input 2: Stitched Together Drift Segment (Y trace)
% Input 3: # of samples per frame

% Output 1: Torsion corrected X trace
% Output 2: Torsion corrected Y trace

%Note: Ctrl-F "Figure" to find figures. There is a figure (and instructions) in the
%middle of the code to visualize correction of an individual frame.

% Norick Bowers, 12/2017.

function [TCorrX,TCorrY] = TorsionCorrection(DriftSegX,DriftSegY,SPF)
% Look for consistent slopes indicating torsional movements
FrameAvg = 10; % # of frames to average torsion over

FWin = round(FrameAvg/2); % Window of frames for averaging
SegFramesX = reshape(DriftSegX,SPF,[]); %Matrix where C=Frames & R=Strips
SegFramesY = reshape(DriftSegY,SPF,[]); %Matrix where C=Frames & R=Strips
TorCorrX = zeros(size(SegFramesX)); %Initialize torsion line samples (built from measured slopes)
TorCorrY = zeros(size(SegFramesY)); %Initialize torsion line samples (built from measured slopes)
AvgSlopeX = zeros(1,length(SegFramesX));
AvgSlopeY = zeros(1,length(SegFramesY));

%Loop through only the frames designated by TorFrames
for aa = 2:length(SegFramesX(1,:))-1
    
    IDX = aa-1:aa+1; %Last, current, and next frame indices
    
    %% X Torsion
    
    %Approximate real motion (X Trace)
    LMean = mean(SegFramesX(:,IDX(1))); %Last frame mean
    CMean = mean(SegFramesX(:,IDX(2))); %Current frame mean
    NMean = mean(SegFramesX(:,IDX(3))); %Next frame mean
    
    LSlope = polyfit(1:2,[LMean,CMean],1); LSlope=LSlope(1); %Slope from last frame
    NSlope = polyfit(1:2,[CMean,NMean],1); NSlope=NSlope(1); %Slope to next frame
    CSlope = mean([LSlope,NSlope]); %Current frame's slope
    
    RealMotion = linspace(-CSlope/2,CSlope/2,SPF); % Real Motion Slope
    RealMotion = RealMotion(:); %Columnize
    
    %Remove real motion and calculate torsion
    ThisFrame = SegFramesX(:,aa); %This frame's motion
    TFC = ThisFrame-RealMotion; %This frame (corrected for real motion)
    
    XTFC{aa} = TFC-TFC(1); %Save copy of frame motions (plotting purposes only)
    
    FitX = polyfit( (1:length(TFC))' , TFC , 1); %Fit torsion slope (Fittmp(1) = slope, Fittmp(2) = y-intercept
    TSlopeX(aa) = FitX(1); %This frame's slope
    
    
    %Apply Smooth Torsion Correction (Indexed FWin frames back)
    if aa > FWin*2 && aa <= length(SegFramesX(1,:))-FWin-1;
        CurrFrame = aa - FWin; %Apply T correction retrospectively
        CurrWin = (aa - (2*FWin)) : (aa); %Window +/- FWin
        CurrWin(FWin+1)= []; %Exclude local slope
        
        AvgSlopeX(CurrFrame) = mean(TSlopeX(CurrWin)); %Avg T Slope of +/- FWin
        TorCorrX(:,CurrFrame) = (linspace(-SPF/2,SPF/2,SPF) .* AvgSlopeX(CurrFrame))' ; %Linear T correction for Frame (aa-FWin)
    end
    
    %     %FIGURE
    %     %In order to visualize how real motion is estimated and removed from
    %     %the trace, break above (at TorCorrX) and then uncomment and plot the
    %     %code below. FigA shows how the real eye motion of a frame is
    %     %calculated as the slope of the nearby frames. FigB shows how this real
    %     %motion is the applied to the current frame's samples so the effects of
    %     %real eye motion does not influence our estimate of torsion. Torsion is
    %     %calculated as the slope of the real-motion-corrected frame (the
    %     %magenta line in FigB). Note that FigA and FigB are not
    %     %necessarily to scale.
    
%     %     %Current Slope (Fanagled for visualization purposes)
%     CSlopetmp = [RealMotion(1)*2 RealMotion(length(RealMotion)/2)*2 RealMotion(end)*2];
%     
%     subplot(1,2,1);
%     plot([LMean CMean],'b','LineWidth',1)
%     hold on; plot([2 3],[CMean NMean],'r','LineWidth',1)
%     hold on; plot(CSlopetmp+CMean,'k','LineWidth',1);
%     hold on; scatter(1,LMean,'b','filled');
%     hold on; scatter(2,CMean,'k','filled');
%     hold on; scatter(3,NMean,'r','filled');
%     legend('Last Frame','Next Frame','Current Frame','Location','Best');
%     ylabel('Arcmin');
%     set(gca,'xticklabel',{'Last Mean','','Current Mean','','Next Mean'});
%     title('Estimated real motion from 3 frame mean positions');
%     
%     subplot(1,2,2);
%     plot(ThisFrame-CMean,'k:','LineWidth',1.5);
%     hold on; plot(RealMotion,'k','LineWidth',1.5);
%     hold on; plot(TFC-CMean,'m:','LineWidth',1.5);
%     legend('Current Frame','Estimated Real Motion','Corrected Frame','Location','Best');
%     xlim([0 SPF]);
%     xlabel('Samples'); ylabel('Normalized Arcmin');
%     title('Real motion correction of current frame');
%     set(gcf,'unit','normalized','outerposition',[0 0 1 1]);
      
    %% Y Torsion
    
    %Clear redundant variables (shouldn't matter, but just in case)
    clear LMean CMean NMean LSlope NSlope CSlope RealMotion ThisFrame TFC
    
    %Approximate real motion (Y Trace)
    LMean = mean(SegFramesY(:,IDX(1))); %Last frame mean
    CMean = mean(SegFramesY(:,IDX(2))); %Current frame mean
    NMean = mean(SegFramesY(:,IDX(3))); %Next frame mean
    
    LSlope = polyfit(1:2,[LMean,CMean],1); LSlope=LSlope(1); %Slope from last frame
    NSlope = polyfit(1:2,[CMean,NMean],1); NSlope=NSlope(1); %Slope to next frame
    CSlope = mean([LSlope,NSlope]); %Current frame's slope
    
    RealMotion = linspace(-CSlope/2,CSlope/2,SPF); % Real Motion Slope
    RealMotion = RealMotion(:); %Columnize
    
    %Remove real motion and calculate torsion
    ThisFrame = SegFramesY(:,aa); %This frame's motion
    TFC = ThisFrame-RealMotion; %This frame (corrected for real motion)
    
    YTFC{aa} = TFC-TFC(1); %Save copy of frame motions (plotting purposes only)
    
    FitY = polyfit( (1:length(TFC))' , TFC , 1); %Fit torsion slope (Fittmp(1) = slope, Fittmp(2) = y-intercept
    TSlopeY(aa) = FitY(1); %This frame's slope
    
    
    %Apply Smooth Torsion Correction (Indexed FWin frames back)
    if aa > FWin*2 && aa <= length(SegFramesY(1,:))-FWin-1;
        AvgSlopeY(CurrFrame) = mean(TSlopeY(CurrWin)); %Avg T Slope of +/- FWin
        TorCorrY(:,CurrFrame) = (linspace(-SPF/2,SPF/2,SPF) .* AvgSlopeY(CurrFrame))' ; %Linear T correction for Frame (aa-FWin)
    end
end

%% Remove Torsion
TCorrX = SegFramesX - TorCorrX; %Apply torsion correction
TCorrX = reshape(TCorrX,[],1); %Reshape to single vector
TCorrX = TCorrX(:);

TCorrY = SegFramesY - TorCorrY; %Apply torsion correction
TCorrY = reshape(TCorrY,[],1); %Reshape to single vector
TCorrY = TCorrY(:);

% %% Figures
% figure;
% plot(DriftSegX,'k');
% hold on;
% H = plot(TCorrX,'m');
% H.Color(4) = 0.75;
% hold on;
% plot(TCorrX+10,'c')
% xlabel('Samples');
% ylabel('Position');
% title('Torsion Correction Vs Original');
% legend('Original','Torsion Corrected','Offset T-corrected','Location','Best');
% xlim([0 length(DriftSegX)]);
% set(gcf,'Units','Normalized','Outerposition',[0 0 1 1]);
% set(gca,'xgrid','on','XTick',[0:SPF:length(DriftSegX)],'XTickLabel', ...
%     0:max([0:SPF:length(DriftSegX)]),'FontSize',10);
% 
% TortmpX = reshape(TorCorrX,[],1);
% TortmpY = reshape(TorCorrY,[],1);
% 
% figure;
% subplot(2,1,1);
% plot(TortmpX,'b');
% title('X Torsion Correction per Frame');
% xlim([0 length(TortmpX)]);
% ylim([min(min([TortmpX TortmpY])) max(max([TortmpX TortmpY]))]);
% subplot(2,1,2);
% plot(AvgSlopeX,'b');
% title('X Average Slope per Frame');
% xlim([0 length(AvgSlopeX)]);
% ylim([min(min([AvgSlopeX AvgSlopeY])) max(max([AvgSlopeX AvgSlopeY]))]);
% 
% figure;
% subplot(2,1,1);
% plot(TortmpY,'r');
% title('Y Torsion Correction per Frame');
% xlim([0 length(TortmpY)]);
% ylim([min(min([TortmpX TortmpY])) max(max([TortmpX TortmpY]))]);
% subplot(2,1,2);
% plot(AvgSlopeY,'r');
% title('Y Average Slope per Frame');
% xlim([0 length(AvgSlopeY)]);
% ylim([min(min([AvgSlopeX AvgSlopeY])) max(max([AvgSlopeX AvgSlopeY]))]);
% 
% 
% 
% 


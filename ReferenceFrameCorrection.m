%This function estimates the reference frame artifact and removes it from
%subsequent drift frames. First, average eye motion is calculated and
%removed from the estimation of the reference frame artifact. The reference
%frame artifact estimation is calculated as the mean movement per strip
%across a number of frames (default 10, ~330ms). After this estimation is
%obtained it is subsequently removed frame-by-frame for all drift frames.

%Input 1: X stream
%Input 2: Y stream
%Input 3: Drift start (indexed to samples)
%Input 4: Drift end (indexed to samples)
%Input 5: Number of samples per frame

%Output 1: Stitched together drift segments (X stream)
%Output 1: Stitched together drift segments (Y stream)
%Output 3: Reference-corrected stitched together drift segments (X stream)
%Output 3: Reference-corrected stitched together drift segments (Y stream)
%Output 4: Reference Artifact X
%Output 5: Reference Artifact Y

%Norick Bowers, 12/2017.
%Adapated from Austin Roorda's "ReferenceFrame_Torsion_Corrector_April2017"

function [DriftSegX,DriftSegY,CorrSegX,CorrSegY,RX,RY] = ReferenceFrameCorrection(SegX,SegY,DriftS,DriftE,SPF)
%% Organize Drift Frames

%Initial Variables:
FrameAvg = 1:10; % # of frames to average for reference frame
MinLngth = 50; %Minimum length of drifts to be stitched together

%Initial drift frames
DframeS = ceil(DriftS/SPF)+1;
DframeE = floor(DriftE/SPF);
DriftSegX = []; DriftSegY = [];

%Add first frame
if DriftS(1) == 1
    DframeS(1) = 1;
end

%If first drift starts and ends in first frame, eliminate first frame (unlikely bug)
if DframeE(1) < DframeS(1)
    
    DframeS(1) = 1;
    DframeE(1) = 1;
    
end

%Stitch all drift frames together
for aa = 1:length(DframeS); %Loop through all frames
    
    tmpS = (DframeS(aa)*SPF)-SPF+1;
    tmpE = (DframeE(aa)*SPF);
    
    %Current Drift Samples
    CurrDX = SegX(tmpS:tmpE)' ;
    CurrDY = SegY(tmpS:tmpE)' ;
    
    %Stitch Together
    if aa ~= 1 && ~isempty(CurrDX) && length(CurrDX)>MinLngth;
        DriftSegX = [DriftSegX , CurrDX - (CurrDX(1)-DriftSegX(end))];
        DriftSegY = [DriftSegY , CurrDY - (CurrDY(1)-DriftSegY(end))];
    elseif aa == 1 && ~isempty(CurrDX);
        DriftSegX = CurrDX;
        DriftSegY = CurrDY;
    end
    
end


SegFramesX = reshape(DriftSegX,SPF,[]);
SegFramesY = reshape(DriftSegY,SPF,[]);
XFrame_M = mean(SegFramesX); %Mean frame motion X
YFrame_M = mean(SegFramesY); %Mean frame motion Y


%% Average Motion Per Frame (Real Movement)

AvgFitX = polyfit((FrameAvg)',XFrame_M(FrameAvg)',1); %Fit Slope
AvgFitY = polyfit((FrameAvg)',YFrame_M(FrameAvg)',1); %Fit Slope



if AvgFitX(1) ~= 0 && AvgFitY(1) ~=0
    
    %Create one frame's worth of samples for fit
    tmp = linspace(FrameAvg(1),FrameAvg(end),SPF);
    CorrectX = (tmp*(AvgFitX(1)/(length(FrameAvg)-1)));
    CorrectX = CorrectX - mean(CorrectX); %Normalize around 0
    CorrectY = (tmp*(AvgFitY(1)/(length(FrameAvg)-1)));
    CorrectY = CorrectY - mean(CorrectY); %Normalize around 0
    
else
    
    %No reference artifact (likely a bug)
    CorrectX = zeros(1,SPF);
    CorrectY = zeros(1,SPF);
    warning(sprintf('Mean motion over %d frames = 0! Likely a bug.',FrameAvg));
    
end

%% Average Motion Per Strip (Reference Artifact)

%Use first drift for strip averaging
XFrametmp = DriftSegX( 1 : (FrameAvg(end)*SPF) ) ;
YFrametmp = DriftSegY( 1 : (FrameAvg(end)*SPF) ) ;

%Initialize Strips Indices
CurrStrip = (0 : SPF : (length(FrameAvg)*SPF)-1) + 1;
for aa = 1:SPF
    
    %Equal Linear Correction
    StripMeanX(aa) = mean( XFrametmp(CurrStrip)-XFrame_M(FrameAvg) );
    StripMeanY(aa) = mean( YFrametmp(CurrStrip)-YFrame_M(FrameAvg) );
      
    CurrStrip = CurrStrip + 1; %Next Strip
    
end


%Remove real motion from estimated reference artifact
ReferenceCorrX = StripMeanX - CorrectX;
ReferenceCorrY = StripMeanY - CorrectY;

%Start from 0
ReferenceCorrX = ReferenceCorrX - ReferenceCorrX(1);
ReferenceCorrY = ReferenceCorrY - ReferenceCorrY(1);

%Apply reference correction (using matrices where C = frames and R = samples)
CorrSegX = reshape(DriftSegX,SPF,[]) - (repmat(ReferenceCorrX,length(XFrame_M),1))' ;
CorrSegY = reshape(DriftSegY,SPF,[]) - (repmat(ReferenceCorrY,length(YFrame_M),1))' ;
CorrSegX = reshape(CorrSegX,[],1)' ; %reshape
CorrSegY = reshape(CorrSegY,[],1)' ; %reshape

RX = ReferenceCorrX; RY = ReferenceCorrY; %Simplify for return

end
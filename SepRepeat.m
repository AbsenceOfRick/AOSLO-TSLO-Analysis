%This function takes the entire X and entire Y stream and highlights
%repeated frames for removal then returns the corrected traces and indices
%of removed samples and frames. Repeated frames are artifacts from
%processfromraw.m.

%Input 1: X stream
%Input 2: Y stream
%Input 3: # of samples per frame

%Output 1: Bad frames
%Output 2: Bad samples

%Norick Bowers, Fall 2016

function [RepeatedFrames,RepeatedSamples] = SepRepeat(xx,yy,SPF)

% Get all samples from all frames
FrameCnt = 1:SPF:length(xx) ;

for aa = 1:length(FrameCnt)-1; %Loop through all frames
    
    %Excise (# of samples per frame) of X stream
    tmpX = xx(FrameCnt(aa):FrameCnt(aa+1)-1);
    tmpY = yy(FrameCnt(aa):FrameCnt(aa+1)-1);
    
    %Ignore repeated NaN values from blink frames
    if aa~=1
        
        %If (samples in this frame) == (samples in last frame) & (not a blink)
        if sum(tmpX - tmpX_Prev) == 0 & sum(tmpY - tmpY_Prev) == 0 & ...
                isempty(find(isnan(tmpX) | isnan(tmpX_Prev))) ;
            Repeat(aa) = 1; %Mark repeats
        else
            Repeat(aa) = 0; 
        end
    else
        Repeat(aa) = 0;
    end
    
    %Save (current frame) to be (last frame) on next loop iteration
    tmpX_Prev = tmpX;
    tmpY_Prev = tmpY;
    
end

RepeatedFrames = find(Repeat==1); %pull out repeated frames
RepeatedFrames = RepeatedFrames-1; %offet 1
RepeatedSamples = [];

%Get sample index for repeated frames
for aa = 1:length(RepeatedFrames);
    
    tmp = (RepeatedFrames(aa) * SPF); %Current frame indexed for samples
    RepeatedSamples = [RepeatedSamples (tmp:tmp+(SPF-1))] ; 
    
end

%Figures? Maybe a subplot with all original frames and all repeated frames?
%How to visualize this?

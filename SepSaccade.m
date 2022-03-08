%This function takes X and Y positions as well as the sampling rate in
%order to seperate saccades into start and end times (samples). Spatial
%positions are first Loess filtered to eliminate abberations arising from
%nystagmus or noise. Saccade onset is defined as when the eye moves beyond
%a speed threshold (set below, default to 90 arcmin/sec) and saccade
%offset is defined as when the eye falls below this threshold. Saccades
%must have a minimum duration. Saccades that happen too close together are
%merged into one event. Initial variables can be adjusted below.

%This function will also interpote single jumps in the data greater than a
%threshold (set below). This is to remove obvious artifacts.
% (See "Spike Interpolation" below section for details.)

%Input:
%Unfiltered x position
%Unfiltered y position
%Sampling rate (in Hertz)
%Window size for loess filter

%Output:
%Saccade start times
%Saccade end times
%Filtered X positions
%Filtered Y positions
%Adjusted X stream (no spikes)
%Adjusted Y stream (no spikes)

%Norick Bowers, Spring 2017

function [Sstart,Send,XFilt,YFilt,xx,yy] = SepSaccades(xx,yy,SampRate,FiltWindow,SPF)
%% Initialize

% Initial Variables

SpdThresh = 90;%65;%90;%70;%JG 90; %speed threshold for saccades (arcmin/sec, default 90) 90<=>1.5
MinDur = 5; %Minimum duration of saccade (in ms, default 5)
MaxDur = 100;%JG150; %Maximum duration of saccade (in ms, default 100)
MinGap = 50; %Minimum gap between saccades (in ms, default 50)
MinAmp = 1.5;%1.5;%1.2;%1%1.8;%JG 3; %Minimum amplitude of saccade (in arcmin, default 3)

%Adjust Filter Window to sampling rate (51 worked at 32SPF, adjust accordingly)
%FiltWindow = 51;%51;
%% Low-pass filtered x position (for speed analysis only)
Xtmp = xx;
Ytmp = yy;
Xtmp(find(isnan(xx))) = 0;
Ytmp(find(isnan(yy))) = 0;

XFilt = smooth(Xtmp,FiltWindow,'loess');
YFilt = smooth(Ytmp,FiltWindow,'loess');

clear Xtmp; clear Ytmp;

XFilt(find(isnan(xx))) = NaN;
YFilt(find(isnan(xx))) = NaN;

InstVel = abs(diff(sqrt(XFilt.^ 2 + YFilt.^ 2) .* SampRate)); %speed in Arcmin/Sec

%% Saccade Seperation

%Pull out start/end times
SaccVals = find(InstVel > SpdThresh); %index for values > threshold (indexed into InstVel)
if ~isempty(SaccVals)
    SaccVals = SaccVals(:) ; %Make sure it's a column or else they don't concatenate properly
    tmp = diff(SaccVals); tmp(end+1) = NaN; %Keep indices the same length (diff removes 1)
    NC = find(tmp~=1); %Nonconsecutive SaccVals
    Send = [SaccVals(NC)+1;SaccVals(end)] ; %End Values %JGtemp
    Sstart = [SaccVals(1);SaccVals(NC(1:end-1)+1)] ; %Start Values
    Sstart(end+1) = SaccVals(NC(end)); %Add last saccade start
    
    %eliminate saccades that are too short or too long in duration
    tmp = find( ((Send - Sstart) * (1000/SampRate))<MinDur | ((Send - Sstart) * (1000/SampRate))>MaxDur );
    Sstart(tmp) = [];
    Send(tmp) = [];
    
    %Combine saccades that happen too close together
    Sstarttmp = [Sstart ; NaN]; %offset starts
    Sendtmp = [NaN ; Send]; %offset ends
    tmp = find((Sstarttmp - Sendtmp) * (1000/SampRate) < MinGap) ; %difference in offset starts/ends
    Sstart(tmp) = [];
    Send(tmp-1) = [];
    
    %eliminate saccades that are too small in amplitude
    tmp = find(sqrt((xx(Sstart)-xx(Send)).^2 + (yy(Sstart)-yy(Send)).^2) < MinAmp);
    Sstart(tmp) = [];
    Send(tmp) = [];
    
    
else
    Sstart = []; Send = [];
end

%% Figures
% %
% % Instantaneous speeds w/ saccade samples highlighted in red
% SaccSpeeds = zeros(length(InstVel),1);
% SaccSpeeds(SaccVals) = InstVel(SaccVals);
% SaccSpeeds(find(SaccSpeeds==0)) = NaN;
% 
% figure;
% plot(InstVel,'b');
% hold on;
% plot(SaccSpeeds,'r');
% text(0.5,0.8,sprintf('Red are traces which\nexceed speed threshold'),'units','normalized');
% ylabel('Instantaneous Speed');
% xlabel('Samples');
% 
% %Smoothed / Unsmoothed comparison
% figure;
% subplot(2,1,1);
% plot(xx,'c'); hold on;
% plot(XFilt,'k'); hold on;
% title('X Position');
% xlim([0 length(XFilt)]);
% ylabel('H Position');
% legend('Unfiltered','Filtered','Location','Best');
% subplot(2,1,2);
% plot(yy,'c'); hold on;
% plot(YFilt,'k'); hold on;
% title('Y Position');
% ylabel('V Position');
% xlim([0 length(YFilt)]);
% legend('Unfiltered','Filtered','Location','Best');
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);
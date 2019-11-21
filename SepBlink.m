%This function takes x and y data streams and fills in missing data with
%NaN values corrosponding to the length of the blink/dropped frame.

%Note that this function is assuming blink frames are present as
%discontinuities in the time axis. This may not be the case for all
%versions of stabilizefromraw.m. In this case blinks will likely be marked
%as saccades.

%Input 1: X stream.
%Input 2: Y stream.
%Input 3: Time gap between samples in seconds.
%Input 4: Number of samples in a single frame.
%Input 5: Time Axis (with dropped frames representing discrete gaps).

%Output 1: X stream with NaN's replacing missing data.
%Output 2: Y stream with NaN's replacing missing data.
%Output 3: Start of dropped data (indexed for output 1-2).
%Output 4: End of dropped data (indexed for output 1-2).
%Output 5: Original x stream.
%Output 6: Original y stream.

%Norick Bowers, Fall 2016

function [xx,yy,BadS,BadE,xx1,yy1,IVals] = SepBlink(xx,yy,Samp_per_Sec,TimeAx,SPF,IVals)

xx1 = xx; yy1 = yy; %save original data stream

%Find big gaps in time axis
TimeGap = diff(TimeAx);

%Round (so things match)
TimeGap = round(TimeGap,10);
Samp_per_Sec = round(Samp_per_Sec,10);

BadS = find(TimeGap ~= Samp_per_Sec); % Gaps larger than gap between samples @ given sample rate
BadS_Sec = TimeAx(BadS); %Bad start in seconds (for figure only)
BlinkSec = TimeGap(BadS); %time of blink length in seconds
BadE_Sec = BadS_Sec+BlinkSec; %Bad end in seconds (for figure only)
BlinkSamp = floor(BlinkSec/Samp_per_Sec); %time of blink length in samples
BadE = BadS + BlinkSamp; %Blink ends
% BadS = BadS+1; %Adjust
BadE = BadE -1;
cnt = 0; %Counter to adjust original xx index

for aa = 1:length(BadS)
    
    if BadS(aa) ~= BadE(aa) %check if dropped data is 1 sample (error in artificial data?)
        %NaN's to replace missing data
        DropVals = (1:length(BadS(aa):BadE(aa))) .* NaN;
        DropVals = DropVals' ; %Flip vertically
        
        if mod(length(DropVals),SPF) ~= 0 %Sometimes adds 1 extra sample? Just drop it.
            DropVals(end) = [];
        end
        
        %Adjust Indices
        BadS(aa) = BadS(aa)+cnt;
        BadE(aa) = BadE(aa)+cnt;
        
        %Insert rows (http://www.mathworks.com/matlabcentral/fileexchange/9984-insertrows-a-b-ind-)
        [xx] = insertrows(xx,DropVals, BadS(aa) ); %insertrows("into this matrix", "these values", "at this location")
        [yy] = insertrows(yy,DropVals, BadS(aa) );
        [IVals] = insertrows(IVals,DropVals, BadS(aa) );
        
        cnt = cnt+length(DropVals); %Adjust counter
    end
end

%Check if dropped data is 1 sample (error from artificial traces?)
if BadS == BadE;
    BadS =[];
    BadE = [];
end

%Drop blinks that are too short 
%(error from occasional skips in continuous time axis)
MinLength = 30; %thirty samples (it's usually just 1 or 2 samples long)
tmp = find(BadE - BadS < MinLength);
BadS(tmp) = []; BadE(tmp) = [];

% % Uncomment below in order to check original x stream vs adjusted x stream
% figure;
% subplot(2,1,1);
% plot(xx,'b');
% Yh = get(gca,'ylim');
% hold on
% for aa = 1:length(BadS);
%     H = area([BadS(aa) BadE(aa)],[max(Yh) max(Yh)],'FaceColor','k','EdgeAlpha',0);
%     alpha(H,.4); hold on;
%     H = area([BadS(aa) BadE(aa)],[min(Yh) min(Yh)],'FaceColor','k','EdgeAlpha',0);
%     alpha(H,.4); hold on;
% end
% hold off
% xlim([0 length(xx)]);
% title('NaN Values Inserted');
% 
% subplot(2,1,2);
% plot(TimeAx,xx1,'r');
% Yh = get(gca,'ylim');
% hold on;
% for aa = 1:length(BadS_Sec);
%     H = area([BadS_Sec(aa) BadE_Sec(aa)],[max(Yh) max(Yh)],'FaceColor','k','EdgeAlpha',0);
%     alpha(H,.4); hold on;
%     H = area([BadS_Sec(aa) BadE_Sec(aa)],[min(Yh) min(Yh)],'FaceColor','k','EdgeAlpha',0);
%     alpha(H,.4); hold on;
% end
% hold off
% xlim([0 max(TimeAx)]);
% title('Time Axis')
% Yh = get(gca,'ylim');
% set(gcf,'units','normalized','outerposition',[0 0 1 1]);

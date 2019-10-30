%Appendix

function [F PxxM PyyM Lambda DPS_Seq] = MultiTaper_PowSpec(xx,yy,WinSize,SampRate,THB,Overlap)

Overlap = Overlap/100; %Overlap percentage
LastWin = 1; %Initialize
cnt = 0; %Initialize

% Run DPSS once in order to obtain windows and Lambda values
SeqCnt = (2*THB-1); % # of windows with spectral concentration near 1;
[DPS_Seq,Lambda] = dpss(WinSize,THB,SeqCnt); %Slepian Sequences

%Run Windowed Multitaper Analysis
while LastWin + WinSize <= length(xx)
    
    cnt = cnt+1; %Iterate
    
    %Current Samples
    if cnt == 1
        WinStart = 1;
        WinEnd = WinSize;
    else
        WinStart = round( LastWin - (WinSize * Overlap));
        WinEnd = WinStart + WinSize ;
        WinStart = WinStart + 1;
    end
    
    CurrWin = [WinStart : WinEnd];
    Xtmp = xx(CurrWin);
    Ytmp = yy(CurrWin);
    Xtmp = Xtmp - Xtmp(1);
    Ytmp = Ytmp - Ytmp(1);
    
    %Run pmtm
    [Pxx(cnt,:) F] = pmtm(Xtmp,{THB,SeqCnt},WinSize,SampRate);
    [Pyy(cnt,:) F] = pmtm(Ytmp,{THB,SeqCnt},WinSize,SampRate);
    
    %Take amplitude from power
    Pxx(cnt,:) = sqrt(Pxx(cnt,:));
    Pyy(cnt,:) = sqrt(Pyy(cnt,:));
    
    LastWin = CurrWin(end); %Last window end
    
end


%Take Average
SP = size(Pxx);
if SP(1) ~= 1 %Added for long drift segments
    PxxM = mean(Pxx);
    PyyM = mean(Pyy);
else
    PxxM = Pxx;
    PyyM = Pyy;
end

%% Figures
% figure;
% for aa = 1:length(Pxx(:,1))
%     hold on; plot(log10(F),log10(Pxx(aa,:)));
% end
% hold on;
% plot(log10(F),log10(PxxM),'k','LineWidth',2);
% ylabel('Log Amplitude');
% xlabel('Log Frequency');
% title('Amplitude: Each Window and Mean');
% 
% 
% figure;
% for aa = 1:length(DPS_Seq(1,:))
%     subplot(2,1,1);
%     hold on;
%     plot(DPS_Seq(:,aa),'LineWidth',2);
%     title('Slepian Sequences');
%     xlim([1 WinSize]);
%     
%     subplot(2,1,2);
%     hold on;
%     bar(aa,Lambda(aa));
%     ylabel('Spectral Concentration');
%     title('Lambda Values');
%     set(gca,'XTickLabel',[]);
% end
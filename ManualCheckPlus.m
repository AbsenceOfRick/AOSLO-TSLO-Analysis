% This function will manually loop through the entire eye trace and allow
% you to manually select bad data as well as identify missed saccades. The
% GUI will have 4 buttons (defined below) as well as two plots for the X and
% Y eye trace. The GUI will loop through the whole eye trace in a series of
% windows (last input) to allow you to examine the trace closer as it's
% rather hard to see things in the entire trace. The function will output new
% drift and saccade end points based on your selections, as well as next
% indices indicating data that has been manually rejected.
%
% Buttons:
% 1. Make Selection: This calls ginput to make 4 selections that must be made
% in the correct order (beginning (X), end (X), beginning (Y) end (Y). to
% select a portion of the trace.
% 2. Mark Saccade: This button will identify the selected data to be
% classified as a saccade.
% 3. Mark Rejected: This button will identify the selected data to be
% classified as rejected. This will automatically select the nearest frames
% to classify as reject (unless you bisect a saccade with the start/end, then
% the respective start/end will be that saccade, should be a rare occurence).
% 2. Next: This button will open the next chunk of the trace to manually
% examine.
%
% Input:
% X trace
% Y trace
% # of samples per frame
% Saccade starts
% Saccade ends
% Drift starts
% Drift ends
% Blink/dropped starts
% Blink/dropped ends
% Auto-rejected starts
% Auto-rejected ends
% Window size (as percentage of total trace)
%
% Output:
% New saccade starts
% New saccade ends
% Rejected starts
% Rejected ends
% New drift starts
% New drift ends
% Old saccade starts
% Old saccade ends
% Old drift starts
% Old drift ends

function [AllSs,AllSe,Rs,Re,AllDs,AllDe,Ss,Se,Ds,De] = ManualCheckPlus(xx,yy,SPF,Ss,Se,Ds,De,Bs,Be,ARs,ARe,MW,XFilt,YFilt)
%% Initialize

%Initialize new indices (Just pad with zerosfor now)
NewSs = zeros(1,length(Ss));
NewSe = zeros(1,length(Se));
Rs = zeros(1,length(Ss));
Re = zeros(1,length(Se));
AdjVals = zeros(1,length(Rs));

%Windows to the nearest frame
MW = round(length(xx) * MW);
FDem = [0:SPF:length(xx)];
MW = find(abs(FDem-MW) == abs((min(abs(FDem - MW))))); %Frame demarcations
MW = MW(1) * SPF;
Ws = 1:MW:length(xx);
We = Ws - 1;
We(1) = [];
We = [We length(xx)];

SampRate=960;
vel=[diff(XFilt).*SampRate./60, diff(YFilt).*SampRate./60];% for deg/s speed
pythVel=abs(diff(sqrt(XFilt.^ 2 + YFilt.^ 2) .* SampRate))/60;

vel=[vel(1,:);vel];
pythVel=[pythVel(1); pythVel];

%% Display GUI for Windows over entire trace
aa=0; RIdx = []; Cnt = 0; MarkEnd = 1;
while aa < length(Ws)
    
    InputS = []; InputE = [];
    NextPlot = 0; ResetCurr = 0;
    aa = aa+1;
    WinAdj = (MW * aa) - MW;
    
    %Find current saccade/blink/drift values for this window
    Sstmp = Ss(find(Ss<We(aa) & Ss>=Ws(aa))); Setmp = Se(find(Se<We(aa) & Se>=Ws(aa))); %Curr saccades
    Dstmp = Ds(find(Ds<We(aa) & Ds>=Ws(aa))); Detmp = De(find(De<We(aa) & De>=Ws(aa))); %Curr drifts
    Bstmp = Bs(find(Bs<We(aa) & Bs>=Ws(aa))); Betmp = Be(find(Be<We(aa) & Be>=Ws(aa))); %Curr blinks
    ARstmp = ARs(find(ARs<We(aa) & ARs>=Ws(aa))); ARetmp = ARe(find(ARe<We(aa) & ARe>=Ws(aa))); %Curr auto-rejected
    
    %Force it to be vertical arrays for later concatentations
    Sstmp = Sstmp(:)-WinAdj; Setmp = Setmp(:)-WinAdj; Dstmp = Dstmp(:)-WinAdj;
    Detmp = Detmp(:)-WinAdj; Bstmp = Bstmp(:)-WinAdj; Betmp = Betmp(:)-WinAdj;
    ARstmp= ARstmp(:)-WinAdj;ARetmp= ARetmp(:)-WinAdj;
    
    %% JG add
    if length(Dstmp)~=length(Detmp)
        fprintf([ num2str(aa) '/' num2str(length(Ws)) ': # Drift starts and endings different BEFORE processing of first and last elements\n']);
    end
    
    if ~isempty(Sstmp) && ~isempty(Setmp)%JG bug #6
        if  Sstmp(1)>Setmp(1) %Saccades
            Sstmp = [1;Sstmp];
        end
        if Setmp(end)<Sstmp(end)
            Setmp(end+1) = length(xx(Ws(aa):We(aa)));
        end
    end
    
    if ~isempty(Dstmp) && ~isempty(Detmp)
        if  Dstmp(1)>Detmp(1) %| (Dstmp(1)~=1 & Dstmp(1)>Detmp(1))%JG mod %Drift
            Dstmp = [1;Dstmp];
        end
        if Detmp(end)<Dstmp(end)
            Detmp(end+1) = length(xx(Ws(aa):We(aa)));
        end
    end
    
    if ~isempty(Bstmp)
        if  Bstmp(1)>Betmp(1) %Blinks
            Bstmp = [1;Bstmp];
        end
        if Betmp(end)<Bstmp(end)
            Betmp(end+1) = length(xx(Ws(aa):We(aa)));
        end
    end
    
    if ~isempty(ARstmp) && isempty(ARetmp)%JG bug#3
        ARetmp= length(xx(Ws(aa):We(aa)));
    end
    
    if ~isempty(ARstmp)
        
        if  ARstmp(1)>ARetmp(1) %Previously auto-rejected
            ARstmp = [1;ARstmp];
        end
        if ARetmp(end)<ARstmp(end)
            ARetmp(end+1) = length(xx(Ws(aa):We(aa)));
        end
        
    end
    
    
    if length(Dstmp)~=length(Detmp)%JG bug#4
        fprintf([ num2str(aa) '/' num2str(length(Ws)) ': # Drift starts and endings different AFTER processing of first and last elements\n']);
        
        if length(Dstmp)<length(Detmp)
            % we might have an extra event at the end
            if Detmp(end-1)>Dstmp(end)
                Detmp=Detmp(1:end-1);
            end
        else %length(Dstmp)>length(Detmp), then trim the first elements
            if Dstmp(2)<Detmp(1)
                Dstmp=Dstmp(2:end);
            end
            if Dstmp(3)<Detmp(2)
                Dstmp(3)=[];
            end
        end
    end
    
    if length(Dstmp)~=length(Detmp)
        fprintf([ num2str(aa) '/' num2str(length(Ws)) ': # Drift starts and endings different still AFTER processing \nof first and last elements\n']);
        beep;
    end
    
    %if any( (Detmp-Dstmp)==length(xx(Ws(aa):We(aa))) )% if %JG bug#5
    %    fprintf('\ Error: one drift segment the size of the full window\n');
    %end
    
    
    %     %Event starting in one window and ending in another
    %     if length(Bstmp) ~= length(Betmp) && Bstmp(end)>Betmp(end) %Blinks
    %         Betmp(end+1) = length(xx(Ws(aa):We(aa)));
    %     elseif length(Bstmp) ~= length(Betmp) && Bstmp(1)>Betmp(1)
    %         Bstmp = [1 ; Bstmp];
    %     end
    %
    %     if length(Sstmp) ~= length(Setmp) && Sstmp(end)>Setmp(end) %Blinks
    %         Setmp(end+1) = length(xx(Ws(aa):We(aa)));
    %     elseif length(Sstmp) ~= length(Setmp) && Sstmp(1)>Setmp(1)
    %         Sstmp = [1 ; Sstmp];
    %     end
    %
    %       if length(Dstmp) ~= length(Detmp) && Dstmp(end)>Detmp(end) %Blinks
    %         Detmp(end+1) = length(xx(Ws(aa):We(aa)));
    %     elseif length(Dstmp) ~= length(Detmp) && Dstmp(1)>Detmp(1)
    %         Dstmp = [1 ; Bstmp];
    %     end
    
    %Make UI Panels and Buttons
    figure('units','normalized','position',[.1 .1 .7 .7]); %open figure
    UIPH_Fig = uipanel('BackGroundColor','white','Position',[0 0.0 .85 1], ...
        'Title',sprintf('Window %d/%d, %d Samples',aa,length(Ws),MW));
    
    UIPH = uipanel('BackgroundColor','white','Position',[.85 0.35 .15 .3]); %Panel for buttons & 2D plot.
    
    UIPH_2D =    uipanel('BackgroundColor','white','Position',[.85 0.65 .15 .35]);
    UIPH_2Dspeed=uipanel('BackgroundColor','white','Position',[.85 0.0 .15 .35]);
    
    %Make Selection Button
    UIH_Select = uicontrol('parent',UIPH,'style','pushbutton','string','Make Selection','BackgroundColor',[.1,.1,.1],...
        'units','normalized','position',[.01 .61 .49 .3],'Fontsize',12,'ForeGroundColor','w','CallBack',{@UserR});
    %Make Next Button
    UIH_Next = uicontrol('parent',UIPH,'style','pushbutton','string','Next','BackgroundColor',[.8,.8,.8],...
        'units','normalized','position',[.49 .61 .49 .3],'Fontsize',12,'ForeGroundColor','k','CallBack',{@UserR});
    %Accept Button
    UIH_Sacc = uicontrol('parent',UIPH,'style','pushbutton','string','Mark Saccade','BackgroundColor',[.2,.7,.2],...
        'units','normalized','position',[.01 .31 .49 .3],'Fontsize',12,'CallBack',{@UserR});
    %Reject Button
    UIH_Reject = uicontrol('parent',UIPH,'style','pushbutton','string','Mark Rejected','BackgroundColor',[.7,.2,.2], ...
        'units','normalized','position',[.49 .31 .49 .3],'Fontsize',12,'CallBack',{@UserR});
    %Reset Button
    UIH_Reset = uicontrol('parent',UIPH,'style','pushbutton','string','RESET','BackgroundColor',[.1,.7,.7], ...
        'units','normalized','position',[.01 .01 .97 .3],'Fontsize',16,'CallBack',{@UserR});
    
    
    
    %Plot X/Y (parent is first uipanel)
    SP1 = subplot(3,1,1,'parent',UIPH_Fig);
    plot(xx(Ws(aa):We(aa)),'b','parent',SP1);
    title('X Trace');
    ylabel('Arcmin'); %xlabel('Frames');
    xlim(SP1,[1 length(xx(Ws(aa):We(aa)))]);
    pos = get(gca, 'Position');%JGadd
    set(gca,'xgrid','on','XTick',[0:SPF:length(xx)],'XTickLabel', ...
        0:max([0:SPF:length(xx)]),'FontSize',10,'Position',[0.02 pos(2) .96 pos(4)]);
    Yh1 = get(gca,'ylim');
    Ax1 = gca;
    
    %Labels
    hold on
    for bb = 1:length(Sstmp) %Saccades
        try
            H = fill([Sstmp(bb) Setmp(bb) Setmp(bb) Sstmp(bb)],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'m');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    
    for bb = 1:length(Bstmp) %Blinks
        try
            H = fill([Bstmp(bb) Betmp(bb) Betmp(bb) Bstmp(bb)],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'k');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(Dstmp) %Drifts
        try
            H = fill([Dstmp(bb) Detmp(bb) Detmp(bb) Dstmp(bb)],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'c');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(ARstmp) %Previously Auto-rejected
        try
            H = fill([ARstmp(bb) ARetmp(bb) ARetmp(bb) ARstmp(bb)],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'k');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0); hold on;
        end
    end
    
    hold off
    
    SP2 = subplot(3,1,2);
    plot(yy(Ws(aa):We(aa)),'r','parent',SP2);
    title('Y Trace');
    ylabel('Arcmin'); %xlabel('Frames');
    xlim(SP2,[1 length(yy(Ws(aa):We(aa)))]);
    pos = get(gca, 'Position');%JGadd
    set(gca,'xgrid','on','XTick',[0:SPF:length(xx)],'XTickLabel', ...
        0:max([0:SPF:length(xx)]),'FontSize',10,'Position',[0.02 pos(2) .96 pos(4)]);
    Yh2 = get(gca,'ylim');
    Ax2 = gca;
    
    %Labels
    hold on
    for bb = 1:length(Sstmp) %Saccades
        try
            H = fill([Sstmp(bb) Setmp(bb) Setmp(bb) Sstmp(bb)],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'m');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(Bstmp) %Blinks
        try
            H = fill([Bstmp(bb) Betmp(bb) Betmp(bb) Bstmp(bb)],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'k');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(Dstmp) %Drifts
        try
            H = fill([Dstmp(bb) Detmp(bb) Detmp(bb) Dstmp(bb)],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'c');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(ARstmp) %Auto-rejected
        try
            H = fill([ARstmp(bb) ARetmp(bb) ARetmp(bb) ARstmp(bb)],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'k');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0); hold on;
        end
    end
    hold off
    
    SP3 = subplot(3,1,3);
    plot(pythVel(Ws(aa):We(aa)),'Color',[0.6 0 0.6],'parent',SP3);
    title('Pythagorian Velocity Trace');
    ylabel('Velocity (Deg/sec)'); %xlabel('Frames');
    xlim(SP3,[1 length(yy(Ws(aa):We(aa)))]);
    pos = get(gca, 'Position');%JGadd
    set(gca,'xgrid','on','XTick',[0:SPF:length(xx)],'XTickLabel', ...
        0:max([0:SPF:length(xx)]),'FontSize',10,'Position',[0.02 pos(2) .96 pos(4)]);
    Yh3 = get(gca,'ylim');
    Ax3 = gca;
    linkaxes([Ax1,Ax2,Ax3],'x');
    %Labels
    hold on
    for bb = 1:length(Sstmp) %Saccades
        try
            H = fill([Sstmp(bb) Setmp(bb) Setmp(bb) Sstmp(bb)],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'m');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(Bstmp) %Blinks
        try
            H = fill([Bstmp(bb) Betmp(bb) Betmp(bb) Bstmp(bb)],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'k');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(Dstmp) %Drifts
        try
            H = fill([Dstmp(bb) Detmp(bb) Detmp(bb) Dstmp(bb)],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'c');
            set(H, 'FaceAlpha', 0.1,'EdgeAlpha',0); hold on;
        end
    end
    
    for bb = 1:length(ARstmp) %Auto-rejected
        try
            H = fill([ARstmp(bb) ARetmp(bb) ARetmp(bb) ARstmp(bb)],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'k');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0); hold on;
        end
    end
    hold off
    
    
    %JG add------ add the
    %Plot X/Y (parent is first uipanel)
    SP4 = axes(UIPH_2D);%SP1 = plot(3,1,1,'parent',UIPH_Fig);
    plot(SP4,xx(Ws(aa):We(aa)),yy(Ws(aa):We(aa)),'Color',[0.6 0.6 0.6]);%SP1);
    hold on
    for bb = 1:length(Sstmp)
        plot(SP4,xx(Ws(aa)-1+Sstmp(bb):Ws(aa)-1+Setmp(bb)),yy(Ws(aa)-1+Sstmp(bb):Ws(aa)-1+Setmp(bb)),'m')
    end
    for bb = 1:length(Bstmp)
        plot(SP4,xx(Ws(aa)-1+Bstmp(bb):Ws(aa)-1+Betmp(bb)),yy(Ws(aa)-1+Bstmp(bb):Ws(aa)-1+Betmp(bb)),'k')
    end
    for bb = 1:length(Dstmp)
        plot(SP4,xx(Ws(aa)-1+Dstmp(bb):Ws(aa)-1+Detmp(bb)),yy(Ws(aa)-1+Dstmp(bb):Ws(aa)-1+Detmp(bb)),'Color',[0.4 0.4 0.4])
    end
    for bb = 1:length(ARstmp)
        plot(SP4,xx(Ws(aa)-1+ARstmp(bb):Ws(aa)-1+ARetmp(bb)),yy(Ws(aa)-1+ARstmp(bb):Ws(aa)-1+ARetmp(bb)),'k')
    end
    hold off
    title('2D plan');
    %ylabel('Arcmin'); %xlabel('Frames');
    pos = get(gca, 'Position');%JGadd
    set(gca,'xgrid','on','XColor',[0 0 1],'YColor',[1 0 0]);
    %'XTick',,[0:SPF:length(xx)],'XTickLabel', ...
    %   0:max([0:SPF:length(xx)]),'FontSize',10,'Position',[0.02 pos(2) .96 pos(4)]);
    
    
    %Plot X/Y (parent is first uipanel)
    SP6 = axes(UIPH_2Dspeed);%SP1 = subplot(3,1,1,'parent',UIPH_Fig);
    plot(SP6,vel(Ws(aa):We(aa),1),vel(Ws(aa):We(aa),2),'Color',[0.6 0 0.6]);%SP1);
    hold on
    for bb = 1:length(Sstmp)
        plot(SP6,vel(Ws(aa)-1+Sstmp(bb):Ws(aa)-1+Setmp(bb),1),vel(Ws(aa)-1+Sstmp(bb):Ws(aa)-1+Setmp(bb),2),'m')
    end
    for bb = 1:length(Bstmp)
        plot(SP6,vel(Ws(aa)-1+Bstmp(bb):Ws(aa)-1+Betmp(bb),1),vel(Ws(aa)-1+Bstmp(bb):Ws(aa)-1+Betmp(bb),2),'k')
    end
    for bb = 1:length(Dstmp)
        plot(SP6,vel(Ws(aa)-1+Dstmp(bb):Ws(aa)-1+Detmp(bb),1),vel(Ws(aa)-1+Dstmp(bb):Ws(aa)-1+Detmp(bb),2),'Color',[0.4 0.4 0.4])
    end
    for bb = 1:length(ARstmp)
        plot(SP6,vel(Ws(aa)-1+ARstmp(bb):Ws(aa)-1+ARetmp(bb),1),vel(Ws(aa)-1+ARstmp(bb):Ws(aa)-1+ARetmp(bb),2),'k')
    end
    hold off
    title('Velocity');
    %ylabel('Vert. Speed (deg/s)'); %xlabel('Frames');
    pos = get(gca, 'Position');%JGadd
    axis([-2 2 -2 2]);
    set(gca,'xgrid','on','XColor',[0 0 1],'YColor',[1 0 0]);%,'XTick',[0:SPF:length(xx)],'XTickLabel', ...
    %    0:max([0:SPF:length(xx)]),'FontSize',10,'Position',[0.02 pos(2) .96 pos(4)]);
    
    % Put variables needed for nested function in UI button handles
    while NextPlot == 0 && ResetCurr == 0  % Update until next plot/reset is called
        
        UIH_Select.UserData.Idx = RIdx;  UIH_Select.UserData.CurrIt = aa; %Select
        UIH_Select.UserData.Ax1 = Ax1;   UIH_Select.UserData.Ax2 = Ax2; UIH_Select.UserData.Ax3 = Ax3;
        UIH_Select.UserData.Cnt = Cnt;
        
        UIH_Sacc.UserData.Ax1 = Ax1;   UIH_Sacc.UserData.Ax2 = Ax2; UIH_Sacc.UserData.Ax3 = Ax3;%Saccade
        UIH_Sacc.UserData.Idx = RIdx;  UIH_Sacc.UserData.CurrIt = aa;
        UIH_Sacc.UserData.NewSs = NewSs;  UIH_Sacc.UserData.NewSe = NewSe;
        UIH_Sacc.UserData.InputS = InputS;  UIH_Sacc.UserData.InputE = InputE;
        UIH_Sacc.UserData.Cnt = Cnt;
        
        UIH_Reject.UserData.Idx = RIdx;  UIH_Reject.UserData.CurrIt = aa; %Reject
        UIH_Reject.UserData.Ax1 = Ax1;   UIH_Reject.UserData.Ax2 = Ax2; UIH_Reject.UserData.Ax3 = Ax3;
        UIH_Reject.UserData.Rs = Rs;   UIH_Reject.UserData.Re = Re;
        UIH_Reject.UserData.InputS = InputS;  UIH_Reject.UserData.InputE = InputE;
        UIH_Reject.UserData.Cnt = Cnt;
        
        UIH_Reset.UserData.Idx = RIdx;   UIH_Reset.UserData.CurrIt = aa; %Reset
        UIH_Reset.UserData.ResetCurr = ResetCurr;
        
        UIH_Next.UserData.NextPlot = NextPlot; UIH_Next.UserData.CurrIt = aa; %Next
        
        uiwait %Wait for UI
    end
    
    %Adjustments, resets, etc
    
    if ResetCurr == 0 %If next plot
        
        MarkEnd = max(max( [find(Rs ~= 0) find(NewSs ~= 0)] ));
        
        if isempty(MarkEnd)
            MarkEnd = 0;
        end
        
        if aa == 1 %Adjustments for new window positions
            AdjVals(1:MarkEnd) = 1;
        else
            AdjVals(LastEnd : MarkEnd) = WinAdj;
        end
        
    else %If reset
        
        if aa == 0 %Reset Current Values Only
            
            Rs(1:MarkEnd) = 0; NewSs(1:MarkEnd) = 0;
            
        elseif LastEnd ~= MarkEnd
            
            Rs(LastEnd:MarkEnd) = 0;
            NewSs(LastEnd:MarkEnd ) = 0;
            
        end
    end
    LastEnd = MarkEnd + 1; %Previous ending for manually marked data
    
    close all %Close current plot
    
end

%% Modify Indices

%Reindex adjustment for windows & remove zeros
Rs(find(Rs == 0)) = NaN; Rs = Rs + AdjVals;
NewSs(find(NewSs == 0)) = NaN; NewSs = NewSs + AdjVals;
Rs(find(isnan(Rs))) = []; NewSs(find(isnan(NewSs))) = [];
Re(find(Re == 0)) = NaN; Re = Re + AdjVals;
NewSe(find(NewSe == 0)) = NaN; NewSe = NewSe + AdjVals;
Re(find(isnan(Re))) = []; NewSe(find(isnan(NewSe))) = [];

AllDs = Ds(:); AllDe = De(:); AllSs = Ss(:); AllSe = Se(:);

%Eliminate any auto-identified saccades that overlap with manual ones
if ~isempty(NewSs)
    
    for ii = 1:length(NewSs)
        
        if NewSs(ii) > NewSe(ii) %If starts before end, user probably clicked backwards, just adjust for it.
            tmpS = NewSs; tmpE = NewSe;
            NewSs = tmpE; NewSe = tmpS;
        end
        
        DiffS = abs(NewSs(ii)-AllSs);
        
        Nearest = find(DiffS == min(DiffS)); %Nearest saccade
        
        CurrS = NewSs(ii):NewSe(ii);
        NearS = AllSs(Nearest):AllSe(Nearest);
        
        ISS = intersect(CurrS,NearS); %Intersecting values
        
        if ~isempty(ISS)
            
            %Fix Drift first
            DIdxs = find(AllSe(Nearest)+1 == AllDs); %Nearest drift start
            DIdxe = find(AllSs(Nearest)-1 == AllDe); %Nearest drift end
            
            if ~isempty(DIdxs) && ~isempty(DIdxe)
                AllDs(DIdxs) = NewSe(ii)+1;
                AllDe(DIdxe) = NewSs(ii)-1;
                
                %Fix Saccade
                AllSs(Nearest) = NewSs(ii);
                AllSe(Nearest) = NewSe(ii);
                
            else %Saccade improperly place (on blink, two saccades, etc)
                NewSs(ii) = [];
                NewSe(ii) = [];
            end
            
        end
        ISS = []; %Black intersecting values again
        
    end
end
%End NRB New Stuff

for ii = 1:length(NewSs)
    
    %Adjust drifts for new saccade inclusion
    AllDs = unique(sort([AllDs ; NewSe(ii)+1]));
    AllDe = unique(sort([AllDe ; NewSs(ii)-1]));
    
    %JG bug#2 different sizes for Ds and De
    if length(AllDs)~=length(AllDe)
        fprintf('ManualCheckPlus(): different length of Drifts starts and ends\n');
    end
    
    %Adjust saccades for new saccade inclusion
    AllSs = unique(sort([AllSs ; NewSs(ii)]));
    AllSe = unique(sort([AllSe ; NewSe(ii)]));
    
    %JS EDITS: Get all the values withing these newly identified
    %saccades
    %NewSse=[NewSse;AllSs(aa):AllSe(aa)];%JGadd TO DO
    
end

%     AddedD = find(AllDs==NewSe+1);
%     AddedS = find(AllSs==NewSs);

%     AddedD = [];
%     AddedS = [];



%Mark rejected data to the nearest frame, redefine drifts and saccades accordingly
FDem(1) = []; %Remove the 0 value
for ii = 1:length(Rs)
    
    %Reindex to the nearest frames
    NFS = find(abs(FDem-Rs(ii)) == (min(abs(FDem-Rs(ii))))); %Nearest Frame Start
    NFS = NFS(1);
    if FDem(NFS)>Rs(ii) %Nearest PREVIOUS frame only
        NFS = NFS-1;
    end
    Rs(ii) = FDem(NFS);
    
    NFE = find(abs(FDem-Re(ii)) == (min(abs(FDem-Re(ii)))));  %Nearest Frame End
    NFE = NFE(end); %break tie to furthest frame
    if FDem(NFE)<Re(ii) %Nearest FOLLOWING frame only
        NFE = NFE+1;
    end
    Re(ii) = FDem(NFE);
    
    RSamps = Rs(ii):Re(ii); %samples of rejected data
    
    
    %Eliminate Saccades within rejected range
    Iss = intersect(AllSs,RSamps); %intersecting Ss
    if ~isempty(Iss)
        for tt = 1:length(Iss)
            DropSs(tt) = find(AllSs==Iss(tt));
        end
    else
        DropSs = [];
    end
    
    Ise = intersect(AllSe,RSamps); %intersecting Se
    if ~isempty(Ise)
        for tt = 1:length(Ise)
            DropSe(tt) = find(AllSe==Ise(tt));
        end
    else
        DropSe = [];
    end
    
    Sdrop = intersect(DropSs,DropSe);
    AllSs(Sdrop) = NaN; AllSe(Sdrop) = NaN; %drop saccades
    
    %If a saccade start/end is without it's corrosponding end/start
    if length(DropSs)>length(DropSe) %starts inside, ends outside
        Re(ii) = AllSe(DropSs(end)); %Move rejected area end to saccades end
        AllSe(DropSs(end)) = NaN; %Drop whole saccade
    end
    if length(DropSe)>length(DropSs) %starts outside, ends inside
        Rs(ii) = AllSs(DropSe(1)); %Move rejected area start to saccades start
        AllSs(DropSe(end)) = NaN; %Drop whole saccade
    end
    
    
    
    
    %Adjust drift starts/ends withing selected range
    Ids = intersect(AllDs,RSamps); %intersecting Ds
    if ~isempty(Ids)
        for tt = 1:length(Ids)
            DropDs(tt) = find(AllDs==Ids(tt));
        end
    else
        DropDs = [];
    end
    
    
    Ide = intersect(AllDe,RSamps); %intersecting Se
    if ~isempty(Ide)
        for tt = 1:length(Ide)
            DropDe(tt) = find(AllDe==Ide(tt));
        end
    else
        DropDe = [];
    end
    
    
    Ddrop = intersect(DropDs,DropDe); %Drifts entirely encapsulated by rejected area
    AllDs(Ddrop) = NaN; AllDs(Ddrop) = NaN;
    
    %Starts inside, ends outside
    for tt = 1:length(DropDs)
        if isempty(find(DropDs(tt) == Ddrop))
            AllDs(DropDs(tt)) = Re(ii)+1;
        end
    end
    
    %Starts outside, ends inside
    for tt = 1:length(DropDe)
        if isempty(find(DropDe(tt) == Ddrop))
            AllDe(DropDe(tt)) = Rs(ii)-1;
        end
    end
    
    %Rejected region wholly within one drift
    if isempty(Ide) && isempty(Ids)
        
        NDI = find(abs(AllDs-Rs(ii)) == (min(abs(AllDs-Rs(ii)))));
        DSamps = AllDs(NDI):AllDe(NDI);
        
        if ~isempty(intersect(DSamps,RSamps)) %Some part of nearest drift intersects
            
            AllDs = unique(sort([AllDs ; Re(ii)+1]));
            AllDe = unique(sort([AllDe ; Rs(ii)-1]));
            
        end
    end
    
    
    DropDe=[];DropDs=[];DropSs=[];DropSe=[]; %Clear placeholder vars
    
end

%% UserResponse nested function
    function UserR(Source,~) %Can only take handle & event data inputs?
        aa = Source.UserData.CurrIt;
        
        
        if strcmp(Source.String,'Make Selection') %If make selection button
            [Xinput,YInput] = ginput(4);
            Ax1 = Source.UserData.Ax1;
            Ax2 = Source.UserData.Ax2;
            Ax3 = Source.UserData.Ax3;
            Cnt = Source.UserData.Cnt +1;
            InputS(Cnt) = nanmean([Xinput(1),Xinput(3)]); %Marked Start
            InputE(Cnt) = nanmean([Xinput(2),Xinput(4)]); %Marked End
            
            %Fill in GUI for visualization of selected points
            axes(Ax1);
            hold on;
            H = fill([InputS(Cnt) InputE(Cnt) InputE(Cnt) InputS(Cnt)],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'k');
            set(H, 'FaceAlpha', 0.2,'EdgeAlpha',0);
            
            axes(Ax2);
            hold on
            H = fill([InputS(Cnt) InputE(Cnt) InputE(Cnt) InputS(Cnt)],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'k');
            set(H, 'FaceAlpha', 0.2,'EdgeAlpha',0);
            
            axes(Ax3);
            hold on
            H = fill([InputS(Cnt) InputE(Cnt) InputE(Cnt) InputS(Cnt)],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'k');
            set(H, 'FaceAlpha', 0.2,'EdgeAlpha',0);
            
            %Assign variables back into common workspace
            assignin('base','InputS',InputS);
            assignin('base','InputE',InputE);
            assignin('base','Cnt',Cnt);
            uiresume
            
            
        elseif strcmp(Source.String,'Mark Rejected') %If rejected button...
            
            %Inputs
            InputStmp = Source.UserData.InputS(end);
            InputEtmp = Source.UserData.InputE(end);
            Cnt = Source.UserData.Cnt;
            Rs = Source.UserData.Rs;
            Re = Source.UserData.Re;
            
            %Save Rejected Indices
            Rs(Cnt) = round(InputStmp); Re(Cnt) = round(InputEtmp);
            
            %Fill
            axes(Ax1);
            hold on;
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'r');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            axes(Ax2);
            hold on
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'r');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            axes(Ax3);
            hold on
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'r');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            %Reassign to original workspace
            assignin('base','Rs',Rs);
            assignin('base','Re',Re);
            uiresume
            
        elseif strcmp(Source.String,'Mark Saccade') %If saccade button...
            
            %Inputs
            InputStmp = Source.UserData.InputS(end);
            InputEtmp = Source.UserData.InputE(end);
            Cnt = Source.UserData.Cnt;
            NewSs = Source.UserData.NewSs;
            NewSe = Source.UserData.NewSe;
            
            
            %Save New Saccade Indices
            NewSs(Cnt) = round(InputStmp); NewSe(Cnt) = round(InputEtmp);
            
            %Fill
            axes(Ax1);
            hold on;
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh1) max(Yh1) min(Yh1) min(Yh1)],'m');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            axes(Ax2);
            hold on
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh2) max(Yh2) min(Yh2) min(Yh2)],'m');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            axes(Ax3);
            hold on
            H = fill([InputStmp InputEtmp InputEtmp InputStmp],[max(Yh3) max(Yh3) min(Yh3) min(Yh3)],'m');
            set(H, 'FaceAlpha', 0.3,'EdgeAlpha',0);
            
            %Reassign to original workspace
            assignin('base','NewSs',NewSs);
            assignin('base','NewSe',NewSe);
            uiresume
            
            
        elseif strcmp(Source.String,'Next') %If next button
            NextPlot = Source.UserData.NextPlot;
            NextPlot = 1;
            assignin('base','NextPlot',NextPlot);
            uiresume
            
            
        elseif strcmp(Source.String,'RESET') %if Reset button...
            ResetCurr = Source.UserData.ResetCurr;
            ResetCurr = 1;
            aa = aa-1; %Redo current iteration
            assignin('base','aa',aa);
            assignin('base','ResetCurr',ResetCurr);
            uiresume
            
        end
        
    end

end
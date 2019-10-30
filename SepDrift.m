%This function takes the start and ends  of all other events (Saccades and dropped
%frames as of the original version) and defines drifts as anything not
%defined as a different event.

%Input:
%Pairs of event-starts and event-ends seperated as different variables.
%Last input MUST BE one of the traces.
%First input will be Event1-start and next will be Event1-end, then third
%input will be Event2-start and next will be Event2-end, etc etc.

%Output:
%Drift starts
%Drift ends

%Note: Thus far this will only work if the trace both starts and ends with
%a drift. If it begins with a saccade or blink this code needs to be
%modified.

%Norick Bowers, Fall 2016.

function [DriftS DriftE] = SepDrift(varargin) %variable input

VarS = []; VarE = [];
DriftS = []; DriftE = [];
TriggerWarn = 0;


for aa = 1:length(varargin)-1; %Loop through all but the last
    
    if isempty(find(isnan(varargin{aa}))) && ~isempty(varargin{aa});
        
        varargin{aa} = varargin{aa}(:) ;
        
        if mod(aa,2)~= 0
            VarS = [VarS ; varargin{aa}]; %Event(aa) start
        else
            VarE = [VarE ; varargin{aa}];%Event(aa) end
        end
        
    end
    
end

Sample = varargin{end} ;

%Merge consecutive events into one
tmpS = sort([VarS']); tmpE = sort([VarE']); %sort

tmpS = [tmpS NaN]; tmpE = [NaN tmpE]; %offset

ConsecEvents = find((tmpS-tmpE)==1 | (tmpS-tmpE==2)); %find 1 | 2 sample difference (why 2?)

tmpS(find(isnan(tmpS))) = []; tmpE(find(isnan(tmpE))) = []; %undo offset

tmpS(ConsecEvents) = []; tmpE(ConsecEvents-1) = []; %combine consecutive events


%Loop through all events
for aa = 1:length(tmpS);
    
    if isempty(find(1==tmpS)) && aa==1; %Trace does not start with dropped frame/saccade
        
        %1:first event
        DriftS(aa) = 1;
        DriftE(aa) = tmpS(aa)-1;
        
    elseif aa~=1 && (tmpE(aa-1) ~= tmpS(aa)-1) ;
        
        %Event(aa) end : Event(aa) start
        DriftS(aa) = tmpE(aa-1)+1;
        DriftE(aa) = tmpS(aa)-1;
        
    else
        DriftS(aa) = NaN;
        DriftE(aa) = NaN;
        warning('Drift marked as NaN, likely error, removing...')
        TriggerWarn = 1;
    end
    
end

%Mark end as drift (if not marked as something else)
if tmpE(end) ~= length(Sample);
    DriftS(end+1) = tmpE(end)+1;
    DriftE(end+1) = length(Sample);
end

%Remove NaN's (usually because trace started with blink/saccade)
if TriggerWarn
    warning(sprintf('Removing %d drifts due to NaN values',length(find(isnan(DriftS)) )));
end
DriftS(find(isnan(DriftS))) = [];
DriftE(find(isnan(DriftE))) = [];

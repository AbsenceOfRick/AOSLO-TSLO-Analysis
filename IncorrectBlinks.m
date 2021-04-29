%This function will take blinks that are incorrectly identified and fix
%them. Often subjects will blink while also making a saccade and this will
%cause the trace to be highlighted as 3 seperate events (Saccade->Blink->Saccade)
%when it's really just one event. This is a relatively rare occurence and
%not something that can be analyzed since it occurs with the eye closed.
%This function will simply mark the entire thing as rejected while
%eliminating them from blink and saccade indices. 

%Input:
%saccade starts
%saccade ends
%blink starts
%blink ends
%X trace
%Y traxe

%Output:
%saccade starts
%saccade ends
%rejected starts
%rejected ends
%blink starts
%blink ends


function [SaccS,SaccE,NewRs,NewRe,DropS,DropE] = IncorrectBlinks(SaccS,SaccE,DropS,DropE,xx,yy)

%Find "Blinks" inside saccades and mark whole thing as rejected

ToDropS = zeros(1,length(SaccS)); %Initialize
ToDropB = zeros(1,length(DropS));
NewRs = []; NewRe = [];

for aa = 1:length(SaccS)
    
    CurrS = SaccS(aa):SaccE(aa); %Current saccade
    CurrS = CurrS(:); %Columnate
    
    if ~isempty(intersect(DropS,CurrS)) %Blink within saccade

        BadBlinks = intersect(DropS,CurrS); %All blinks within current S
        ToDropS(aa) = 1; %Delete Marker

        %Drop all blinks within saccade
        for bb = 1:length(BadBlinks)
            CurrB = find(BadBlinks(bb) == DropS); %Index for DropS
            ToDropB(CurrB) = 1; %Drop blink
        end
    end
    
    %Find Saccades near or within Blinks and expand blink
    if ~isempty(find( (SaccE(aa) + 1) == DropS )) && ...
            isempty(find( (SaccS(aa) - 2) == DropE )) %Saccade ending on a blink but not starting
        ToDropS(aa) = 2; %Move blink start marker
    end
    
    if ~isempty(find( (SaccS(aa) - 2) == DropE )) && ...
            isempty(find( (SaccE(aa) + 1) == DropS )) %Saccade immediately coming out of a blink

        ToDropS(aa) = 3; %Move blink end marker

    end
    
    if ~isempty(find( (SaccE(aa) + 1) == DropS )) && ...
            ~isempty(find( (SaccS(aa) - 2) == DropE ))  %Seperate saccade fully encapsulated by blink
        ToDropS(aa) = 1; %Delete Marker
    end

    

end


%Reindex rejected values
for aa = 1:length(ToDropS)
    
    if ToDropS(aa) == 1 %B in S or S in B, simply mark as deleted
        NewRs = [NewRs;SaccS(aa)]; 
        NewRe = [NewRe;SaccE(aa)];
    end
    
    if ToDropS(aa) == 2 %Saccade ending on a blink but not starting
        
        i = find((SaccE(aa) + 1) == DropS);
        DropS(i) = SaccS(aa);
        
    end
    
    if ToDropS(aa) == 3 %Saccade immediately coming out of a blink
        
        i = find((SaccS(aa) - 2) == DropE);
        DropE(i) = SaccE(aa);
        
    end
end



SaccS(find(ToDropS~=0)) = [];
SaccE(find(ToDropS~=0)) = [];

DropS(find(ToDropB~=0))=[];
DropE(find(ToDropB~=0))=[];


    DropRejS = intersect(DropS-1,NewRe);
    DropRejE = intersect(DropE+2,NewRs);
    ToDropBs = zeros(1,length(DropS));
    ToDropBe = zeros(1,length(DropE));
    ToDropRej = zeros(1,length(DropRejS));


 for aa = 1:length(DropRejS)
     ToDropRej(aa) = find(NewRe == DropRejS(aa));

     ToDropBs(find(DropRejS(aa) == DropS-1)) = 1;
     ToDropBe(find(DropRejE(aa) == DropE+2)) = 1;


 end

 if ~isempty(ToDropRej)
    NewRs(ToDropRej) = [];
    NewRe(ToDropRej) = [];
 end

DropS(find(ToDropBs)) = [];
DropE(find(ToDropBe)) = [];







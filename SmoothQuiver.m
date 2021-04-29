%This function is not currently implemented as it requires offsets to align
%to a master retinal image. 

%This function creates smoothed quiver plots. It takes the X Y U V inputs
%that you want to pass into the default quiver function and creates a
%SizeXSize grid of smoothed vectors. A point is included in the grid so
%long as (MinVals) real data points fall within 1std of the X and 1std of the Y.
%The vector in the smoothed data point is simply the average of the
%(MinVals) real data points.


function [XYUV_Smooth] = SmoothQuiver(XS,YS,XE,YE,Size,MinVals,Spacing,Window)

XYUVtmp = []; %Init
It = 1;

%For recentering everything
CX = round(nanmedian(XS));
CY = round(nanmedian(YS));
%XY values for smoothed grid
X = (CX-Size):(CX+Size);
Y = (CY-Size):(CY+Size); %Y is reversed in images
Adj = 10; %For readjusting sizes, used below.


if nanstd(XS) > Size || nanstd(YS) > Size
    Size = round( max(max( [std(YS),std(XS)] )) );
    fprintf(sprintf('Smooted Quiver Plot Area Extended to %d',Size));
end

%Search window (1std)
if exist('Window') ~= 1
    WinX = nanstd(XS);
    WinY = nanstd(YS);
else
    WinX = Window;
    WinY = Window;
end

for xi = 1:Spacing:length(X)
    for yi = 1:Spacing:length(Y)
        
        %SS within 1std of current point on grid
        XVals = find( ((XS > X(xi)-WinX) & XS < (X(xi)+WinX)) );
        YVals = find( ((YS > Y(yi)-WinY) & YS < (Y(yi)+WinY)) );
        UseVals = intersect(XVals,YVals);
        
        if numel(UseVals) >= MinVals
            U(It) = nanmean(XE(UseVals)-XS(UseVals));
            V(It) = nanmean(YE(UseVals)-YS(UseVals));
            
            %Find (MinVals) SaccS within +/- 1STD
            XYUVtmp = [ XYUVtmp ; X(xi), Y(yi), U(It), V(It) ];
        else
            XYUVtmp = [XYUVtmp ; X(xi), Y(yi), NaN, NaN ];
        end
        
        It = It+1;
        
    end
    
end

%Eliminate NaN's
XYUV_Smooth(:,1) = XYUVtmp(find(~isnan(XYUVtmp(:,3))),1);
XYUV_Smooth(:,2) = XYUVtmp(find(~isnan(XYUVtmp(:,3))),2);
XYUV_Smooth(:,4) = XYUVtmp(find(~isnan(XYUVtmp(:,3))),4);
XYUV_Smooth(:,3) = XYUVtmp(find(~isnan(XYUVtmp(:,3))),3);



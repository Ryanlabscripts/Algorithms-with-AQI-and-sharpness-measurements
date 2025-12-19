function [lhSharp]=imSharp10(imgResult,sharp_line,sharp_start,sharp_RHe)

% Image sharpness for inclusion in batch processing

% -9993 value for lhSharp removed.
% If LOWER threshold transition is not found, then assume LH edge of image
% (i.e. 1) as position for threshold crossing.
% Additional check also made about location for max value; warning issued
% if the max value is found at sharp_RHe.
% Parameter added for toggle of 'Section' plot; plotON='True' or 'False'
% plotON='False';
plotON=0;

% 'Section' plot reinstated. New figure window for each call to imSharp

% ERROR checks:
% Check made for max value occurring at edge of image; indicated by lhSharp
% value of -9991
% Failure to find UPPER threshold transition indicated by lhSharp value of -9992
% Failure to find LOWER threshold transition indicated by lhSharp value of -9993


% 'sharp_line,sharp_RHe' added as input variables

%-----------------------------------------------------
% C is the y-value for the 'section line'
C=sharp_line;
m=0; % Gradient of 'section line' if needed in the future
%-----------------------------------------------------

% Set upper and lower threshold values
tU=0.9; % 90%
tL=0.1; % 10%

[nX,~]=size(imgResult);

% Pre-define for 'section line'
x=zeros(1,nX);
y=zeros(1,nX);
% Pre-define for 'extracted line from image'
tmp=zeros(1,nX);

for i=1:nX
    % 'section line'
    x(i)=i;
    y(i)=round(m*x(i)+C);
    % Extracted line from image
    tmp(i)=imgResult(x(i),y(i));
end

if plotON
% Plot 'section'
    figure; % Open a new default figure window
    plot(1:nX,tmp); hold on;
end

%-----------------------------------------------------
% Request input of search extent for LH edge of object
% disp(' ');
% sE=input('Enter x-axis search extent: ');
sE=sharp_RHe;
%-----------------------------------------------------
sS=sharp_start;

% Find max value in the section of data up to sE
[lhUv,lhI]=max(tmp(sS:sE));
lhI=lhI+sS;
if plotON
    plot(lhI,lhUv,'r*'); % Show max value on the plot
end

if lhI == 1
    fprintf('\n\nERROR: Max value found is at LH edge of image.\n\n')
    lhSharp=-9991;
    hold off
    return
elseif lhI == sharp_RHe
    fprintf('\nWARNING: Max value found is at edge of x-axis search extent.')
    fprintf('\nConsider increasing value of sharp_RHe: current value is %d.\n',sharp_RHe)
end

% Search back from LH max until value falls below upper threshold of max value
test=lhUv*tU;
i=lhI-1;
while tmp(i)>test
    i=i-1;
    if i==0
        fprintf('\n\nERROR: LH value does NOT fall below UPPER threshold of max value.\n\n')
        lhSharp=-9992;
        hold off
        return
    end
end
lhU=i+1;
if plotON
    plot(lhU,test,'g*')
end

% Search back from LH max until value falls below lower threshold of max value
test=lhUv*tL;
i=lhU-1;
while tmp(i)>test
    i=i-1;
    if i==0
        % If LOWER threshold transition is not found, then assume LH edge of image
        % (i.e. lhL=1) as position for threshold crossing.
        fprintf('\nWARNING: LH value does NOT fall below LOWER threshold of max value.')
        fprintf('\nLH edge of image assumed as position for threshold crossing.\n')
        break
    end
end
lhL=i+1;

if plotON
    plot(lhL,test,'c*')
    hold off
end

lhSharp=lhU-lhL;

end


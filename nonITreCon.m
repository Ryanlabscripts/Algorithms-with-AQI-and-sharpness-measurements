function [imgResult,tTaken,SharpVal,AqiVal]=nonITreCon(proj,geo,angles,method,lH,result_image,sharp_line,sharp_start,sharp_RHe)
%
% Sharpness & AQI values calculated for selected image

%--------------------------------------------------------------------------
%   'method' will be used to select the specific reconstruction script
%   Possible 'method' options are:
opts=     {'FDK','DUMMY'};

% Remove any leading or trailing spaces and convert to uppercase
method=upper(strtrim(method));

% Check to see if 'method' is an allowable option
if ~contains(opts,method)
    error(['Non-iterative reconstruction method "' method '" not available' ]);
end

%--------------------------------------------------------------------------

%   lH = level of beam hardening to be applied (integer 1 to 6)
if abs(lH) > 6
    error('ERROR: Invalid input for beam hardening level (<=6).')
end

% Apply beam hardening
fprintf('\n ===== \n');
fprintf('Applying beam hardening. \n');
% Polynomial coefficients
a=[0, 0, 0, 0, 0, 0];
b=[1, 0.75, 0.5, 0.2, 0, 0];
c=[0, 0.25, 0.5, 0.8, 0.8, 0.2];
d=[0, 0, 0, 0, 0.2, 0.8];
scale=[1, 1.32, 1.94, 4.4, 10, 139];
if lH==1 % No beam hardening
    projH=proj;
else
    [projH] = beamH(proj,a(lH),b(lH),c(lH),d(lH),scale(lH));
end

%--------------------------------------------------------------------------

% Check all angles <= 0
if (max(angles)>0) && min(angles)>=0
    % All angles >= 0, so this correction is REQUIRED
    disp('Correcting sign of "angles" vector; all values now negative.')
    angles=-angles;
elseif (max(angles)<=0) && min(angles)<0
    % All angles <= 0, so no action required
else
    disp('Positive & negative values found in "angles"; CHECK.')
    return
end

close all

%--------------------------------------------------------------
tic;
switch method
    case 'FDK'
        fprintf('\nNon-iterative reconstruction method selected:  %s\n\n',method);
        % Reconstruct image using FDK
        imgResult=FDK(single(projH),geo,angles,'filter','hann');
%         imgResult=FDK((projH),geo,angles,'filter','hann');

        % Extract image result_image
        result_Req=imgResult(:,:,result_image);

        % Calculate Sharpness Quality Index
        [SharpVal]=imSharp10(result_Req,sharp_line,sharp_start,sharp_RHe);

        % Prepare & calculate AQI Quality Index
        % Values selected for running aqindex (change as required)
        % Current selection recommended by 'Anisotropic blind image quality assessment' paper.
        N=8;
        nod=6;
        firstangle=0;
        angleunits='degree';
        mode='color';
        average='common';
        [AqiVal,~,~]=aqindex(result_Req,N,nod,firstangle,angleunits,mode,average);

    case 'DUMMY'
        fprintf('\nNon-iterative reconstruction method selected:  %s\n\n',method);
        return
        
    otherwise
        error(['Invalid input name: ', method,'. No such option available.']);
end

tTaken=toc; % Time taken for the non-iterative reconstruction

fprintf('===== \n');
fprintf('Time taken for the %s non-iterative reconstruction: %d seconds.\n',...
    method,round(tTaken));
% Alternatively:
disp(['(Time taken :    ',secs2hms(round(tTaken)),')']);
fprintf(' \n');

end


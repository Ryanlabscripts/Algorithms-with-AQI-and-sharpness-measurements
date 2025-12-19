function [xGM,QQ,lhSharp,resultImageOut]= CGLS_Sharpness_mk9_final(proj,geo,angles,lH,niter,Freewheel,sharp_line,sharp_start,sharp_RHe,varargin)
%
% Multiple slice output
% A text file named 'sliceData.txt' contains the numbers of the slices
% required.
% This file must contain at least 1 number.
% e.g.
%      sliceData=[123,456,789]
% The first number is the number of the slice that will be used for
% stopping assessment. 
% resultImageOut now contains all the requested slice images for all of the
% iterations.
% slice data
imgSlices=load('sliceData.txt','-ascii');
[~,nSlices]=size(imgSlices);
if nSlices<1
    return
end
% Option to 'freewheel' iterations and ignore n2Stop added.
% itFreewheel needs to be set in script call; place value (separated by
% comma) just after max number of iterations
% if Freewheel = 1, 'freewheel' and ignore n2Stop
% if Freewheel = 0, do not 'freewheel' and use n2Stop as previously
Freewheel=logical(Freewheel); % Force itFreewheel input to be logical
%
% Initially, assume NO freewheel
itFreewheel=false;
% CGLS_CBCT solves the CBCT problem using the conjugate gradient least
% squares
% 
%  CGLS_CBCT(PROJ,GEO,ANGLES,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
% 
%  CGLS_CBCT(PROJ,GEO,ANGLES,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
% 
% 
%  'Init'    Describes diferent initialization techniques.
%             * 'none'     : Initializes the image to zeros (default)
%             * 'FDK'      : intializes image to FDK reconstrucition
%             * 'multigrid': Initializes image by solving the problem in
%                            small scale and increasing it when relative
%                            convergence is reached.
%             * 'image'    : Initialization using a user specified
%                            image. Not recomended unless you really
%                            know what you are doing.
%  'InitImg'    an image for the 'image' initialization. Avoid.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri 
%--------------------------------------------------------------------------

% Preallocate for AQI assessment
QQ=zeros(niter,1);
Xval=(1:niter)';

% Preallocate for sharpness assessment
lhSharp=zeros(niter,1);

% Values selected for running aqindex (change as required)
% Current selection based upon 'Anisotropic blind image quality assessment' paper.
N=8;
nod=6;
firstangle=0;
angleunits='degree';
mode='color';
average='common';

%--------------------------------------------------------------------------

[ip,iq,ir]=size (proj);
disp(['size(proj) is:  ',num2str(ip),' by ',num2str(iq),' by ',num2str(ir)]);
clear ip iq ir;

%--------------------------------------------------------------------------

%%

%--------------------------------------------------------------------------
%   lH = level of beam hardening to be applied (integer 1 to 6)
if abs(lH) > 6
    error('ERROR: Invalid input for beam hardening level (<=6).')
end


tic;

% Apply beam hardening
fprintf('\n===== \n');
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
clear proj
% Now reset proj to the beam hardened values
proj=projH;

tTaken=toc; % Time at completion of beam hardening
fprintf('===== \n');
disp(['Time at completion of beam hardening :    ',secs2hms(round(tTaken))]);
fprintf(' \n');

%--------------------------------------------------------------------------

[verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,varargin);

measurequality=~isempty(QualMeasOpts);

qualMeasOut=zeros(length(QualMeasOpts),niter);

% //doi: 10.1088/0031-9155/56/13/004

r=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids);
p=Atb(r,geo,angles,'matched','gpuids',gpuids);
gamma=norm(p(:),2)^2;

errorL2=zeros(1,niter);

% Preallocate for frame results
[ip,iq,ir]=size(x);
disp(['size(x) is:  ',num2str(ip),' by ',num2str(iq),' by ',num2str(ir)]);
resultImageOut=zeros(ip,iq,niter,nSlices); % Base size of (resultImageOut) on size of (x)
clear ip iq ir;
[ip,iq,ir,is]=size (resultImageOut);
disp(['size(resultImageOut) is:  ',num2str(ip),' by ',num2str(iq),' by ',num2str(ir),' by ',num2str(is)]);
clear ip iq ir is;
tIterationResults=zeros(niter,1);
% ii4MinResults=zeros(niter,1);
n2StopResults=zeros(niter,1);
itStop=zeros(niter,1);

for ii=1:niter
    x0 = x;
    if (ii==1 && verbose);tic;end
    
    q=Ax(p,geo,angles,'Siddon','gpuids',gpuids);
    alpha=gamma/norm(q(:),2)^2;
    x=x+alpha*p;

    %     ===============================================================
    %       Investigate 'stopping' based on Sharpness
    %
    % NB: The stopping is based upon a specific image; 'result_image' in this case.
    % Extract image result_image

    result_Req=x(:,:,imgSlices(1));
    resultImageOut(:,:,ii,1)=result_Req(:,:); % Save result image (first slice) at each iteration

    % Calculate Sharpness Quality Index
    [lhSharp(ii)]=imSharp10(result_Req,sharp_line,sharp_start,sharp_RHe);

    % Prepare & calculate AQI Quality Index
    [QQ(ii),~,~]=aqindex(result_Req,N,nod,firstangle,angleunits,mode,average);

    tTaken=toc; % Time at completion of iteration

    % If sharpness=0 or NaN this probably means that imSharp has encountered an
    % ERROR, so stop here.
    if (lhSharp(ii)==0) || (isnan(lhSharp(ii)))
        disp('Sharpness=0 or NaN, probably indicative of error encountered in imSharp9.');
        disp('STOP');
        fprintf(' \n');
        disp(['Iteration :    ',num2str(ii)]);
        disp(['AQI value :    ',num2str(QQ(ii))]);
        disp(['Time at completion of iteration :    ',secs2hms(round(tTaken))]);
        fprintf(' \n');
        return
    end

    disp(['Iteration :    ',num2str(ii)]);
    disp(['AQI value :    ',num2str(QQ(ii))]);
    disp(['Sharpness value :    ',num2str(lhSharp(ii))]);
    disp(['Time at completion of iteration :    ',secs2hms(round(tTaken))]);

% ===========================Revision of stopping rules====================================================================
    if ii==1 % Set-up
        xGM=x;
        n2Stop=niter;

        if itFreewheel
            disp('N.B. "Freewheel" iteration to max number of iterations selected.');
        end
        %
        tIterationResults(ii)=round(tTaken);
        tOld=tTaken;
    else 
        xGM=x;

        % Start consideration of whether to stop further iteration
        if (lhSharp(ii) < 40)&&(ii>2)
            if max(lhSharp(ii-2:ii))==min(lhSharp(ii-2:ii))
                n2Stop=ii;
                disp(['STOP iteration no. reset to: ',num2str(n2Stop)]);
            end
        end
        tIteration=tTaken-tOld;
        disp(['Time for iteration :    ',secs2hms(round(tIteration))]);
        tIterationResults(ii)=round(tIteration);
        tOld=tTaken;
    end
 % ==================================================================================================================

    n2StopResults(ii)=n2Stop;

    fprintf(' \n'); % Print line space

    % Save result image (for remainder of slices) for this iteration
    for iSlice=2:nSlices
        result_Req=x(:,:,imgSlices(iSlice));
        resultImageOut(:,:,ii,iSlice)=result_Req(:,:);
    end

    % Check whether to stop iteration process
    if (not(itFreewheel) && ii==n2Stop)
        if(Freewheel)
            itStop(ii)=ii;
            fprintf('\n\nQuality achieved. Iteration would stop at this point if "Freewheel" had not been requested. \n');
            % Update status of itFreewheel and continue through to max number of
            % iterations
            itFreewheel=true;
        end
    end

    % Check whether to stop iteration process
    if (not(itFreewheel) && ii==n2Stop) || (itFreewheel && ii==niter)
        if (not(itFreewheel) && ii==n2Stop)
            fprintf('\n\nQuality achieved. Stopping further iterations. \n');
        else
            fprintf('\n\nFull set of iterations completed. \n');
        end
        disp(' ');
        % disp(['Minimum sharpness value achieved at iteration: ',num2str(ii4Min)]);
        % Now plot the progression of AQI throughout the iterations.
        plot(Xval(1:ii),QQ(1:ii),'-r');
        title('AQI Plot');
        xlabel('No. of iterations');
        ylabel('AQI value');

        % Indicate location of maximum AQI value
        [ymax,imax]=max(QQ);
        hold on
        plot(Xval(imax),ymax,'ob');
        hold off

        ylim([0 1.1*ymax])
        shg;
        % Save the plot as a .fig file with a default filename for the time being.
        savefig('CGLSs7_AQI')
        clf('reset')


        % Now plot the progression of sharpness throughout the iterations.
        figure
        plot(Xval(1:ii),lhSharp(1:ii),'-r');
        title('Sharpness Plot');
        xlabel('No. of iterations');
        ylabel('Sharpness value');

        % % Indicate location of minimum sharpness value
        % hold on
        % plot(Xval(ii4Min),lhSharpmin,'ob');
        hold off

        ymax=max(lhSharp);
        ylim([0 1.1*ymax])
        shg;
        % Save the plot as a .fig file with a default filename for the time being.
        savefig('CGLSs7_Sharpness')
        clf('reset')

        % Write iteration timing data to Excel file
        % Times given to nearest second for each iteration
        % Form Excel filename (which includes date & time stamp)
        filename=strcat('tIterationResults',string(datetime('now','Format','dd-MM-yy_HH-mm-SS')),'.xlsx');
        % Create table
        T = table((1:niter)',tIterationResults,n2StopResults,itStop,...
            lhSharp,QQ,'VariableNames',...
            {'Iteration','Iteration time','n2Stop','itStop','Sharpness','AQI'});
        % Write table to Excel spreadsheet
        writetable(T,filename);

        break
    end

    %     ===============================================================

    aux=proj-Ax(x,geo,angles,'Siddon','gpuids',gpuids); %expensive, is there any way to check this better?
    errorL2(ii)=im3Dnorm(aux,'L2');
    
    if measurequality
        qualMeasOut(:,ii)=Measure_Quality(x0,x,QualMeasOpts);
    end

    if ii>1 && errorL2(ii)>errorL2(ii-1)
        % OUT!
       x=x-alpha*p;
       if verbose
          disp(['CGLS stoped in iteration N', num2str(ii),' due to divergence.']) 
       end
       return; 
    end
    % If step is adecuatem, then continue withg CGLS
    r=r-alpha*q;
    
    s=Atb(r,geo,angles,'matched','gpuids',gpuids);
    gamma1=norm(s(:),2)^2;
    beta=gamma1/gamma;
    gamma=gamma1;
    p=s+beta*p;
    
   
     if (ii==1 && verbose)
        expected_time=toc*niter;   
        disp('CGLS');
        disp(['Expected duration   :    ',secs2hms(expected_time)]);
        disp(['Expected finish time:    ',string(datetime('now')+seconds(expected_time))]);   
        disp('');
     end
end

end


%% parse inputs'
function [verbose,x,QualMeasOpts,gpuids]=parse_inputs(proj,geo,angles,argin)
opts=     {'init','initimg','verbose','qualmeas','gpuids'};
defaults=ones(length(opts),1);

% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:CGLS:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
       error('TIGRE:CGLS:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]); 
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
   if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('TIGRE:CGLS:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        case 'init'
            x=[];
            if default || strcmp(val,'none')
                x=zeros(geo.nVoxel','single');
                continue;
            end
            if strcmp(val,'FDK')
                x=FDK(proj,geo,angles);
                continue;
            end
            if strcmp(val,'multigrid')
                x=init_multigrid(proj,geo,angles);
                continue;
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue;
            end
            if isempty(x)
               error('TIGRE:CGLS:InvalidInput','Invalid Init option') 
            end
            % % % % % % % ERROR
        case 'initimg'
            if default
                continue;
            end
            if exist('initwithimage','var')
                if isequal(size(val),geo.nVoxel')
                    x=single(val);
                else
                    error('TIGRE:CGLS:InvalidInput','Invalid image for initialization');
                end
            end
        %  =========================================================================
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:CGLS:InvalidInput','Invalid quality measurement parameters');
                end
            end
         case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE:Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end
        otherwise 
            error('TIGRE:CGLS:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in CGLS()']);
    end
end


end


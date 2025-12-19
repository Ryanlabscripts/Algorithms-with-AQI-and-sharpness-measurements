function [resGM,QQ,lhSharp,resultImageOut]=OS_SART_Sharpness_mk9_final(proj,geo,angles,lH,niter,Freewheel,sharp_line,sharp_start,sharp_RHe,varargin)
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
%   lH = level of beam hardening to be applied (integer 1 to 6)
%

% OS_SART solves Cone Beam CT image reconstruction using Oriented Subsets
%              Simultaneous Algebraic Reconxtruction Techique algorithm
%
%   OS_SART(PROJ,GEO,ALPHA,NITER) solves the reconstruction problem
%   using the projection data PROJ taken over ALPHA angles, corresponding
%   to the geometry descrived in GEO, using NITER iterations.
%
%   OS_SART(PROJ,GEO,ALPHA,NITER,OPT,VAL,...) uses options and values for solving. The
%   possible options in OPT are:
%
%   'BlockSize':   Sets the projection block size used simultaneously. If
%                  BlockSize = 1 OS-SART becomes SART and if  BlockSize = length(alpha)
%                  then OS-SART becomes SIRT. Default is 20.
%
%   'lambda':      Sets the value of the hyperparameter. Default is 1
%
%   'lambda_red':   Reduction of lambda.Every iteration
%                  lambda=lambdared*lambda. Default is 0.95
%
%   'Init':        Describes diferent initialization techniques.
%                  'none'     : Initializes the image to zeros (default)
%                  'FDK'      : intializes image to FDK reconstrucition
%                  'multigrid': Initializes image by solving the problem in
%                               small scale and increasing it when relative
%                               convergence is reached.
%                  'image'    : Initialization using a user specified
%                               image. Not recomended unless you really
%                               know what you are doing.
%   'InitImg'      an image for the 'image' initialization. Aviod.
%
%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
%   'QualMeas'     Asks the algorithm for a set of quality measurement
%                  parameters. Input should contain a cell array of desired
%                  quality measurement names. Example: {'CC','RMSE','MSSIM'}
%                  These will be computed in each iteration.
% 'OrderStrategy'  Chooses the subset ordering strategy. Options are
%                  'ordered' :uses them in the input order, but divided
%                  'random'  : orders them randomply
%                  'angularDistance': chooses the next subset with the
%                                     biggest angular distance with the ones used.
%
% OUTPUTS:
%
%    [img]                       will output the reconstructed image
%    [img,errorL2]               will output the L2 norm of the residual
%                                (the function being minimized)
%    [img,errorL2,qualMeas]      will output the quality measurements asked
%                                by the input 'QualMeas'
%
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


%% Deal with input parameters

[blocksize,lambda,res,lambdared,verbose,QualMeasOpts,OrderStrategy,nonneg, gpuids]=parse_inputs(proj,geo,angles,varargin);
measurequality=~isempty(QualMeasOpts);

qualMeasOut=zeros(length(QualMeasOpts),niter);

if nargout>1
    computeL2=true;
else
    computeL2=false;
end
% does detector rotation exists?
if ~isfield(geo,'rotDetector')
    geo.rotDetector=[0;0;0];
end
%% weigth matrices
% first order the projection angles
[alphablocks,orig_index]=order_subsets(angles,blocksize,OrderStrategy);


% Projection weigth, W
geoaux=geo;
geoaux.sVoxel([1 2])=geo.sVoxel([1 2])*1.1; % a Bit bigger, to avoid numerical division by zero (small number)
geoaux.sVoxel(3)=max(geo.sDetector(2),geo.sVoxel(3)); % make sure lines are not cropped. One is for when image is bigger than detector and viceversa
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'Siddon', 'gpuids', gpuids);  %
W(W<min(geo.dVoxel)/2)=Inf;
W=1./W;


% Back-Projection weigth, V
V=computeV(geo,angles,alphablocks,orig_index, gpuids);

clear A x y dx dz;


%% hyperparameter stuff
nesterov=false;
if ischar(lambda)&&strcmp(lambda,'nesterov')
    nesterov=true;
    lambda=(1+sqrt(1+4))/2;
    gamma=0;
    ynesterov=zeros(size(res),'single');
    ynesterov_prev=ynesterov;
end
%% Iterate
errorL2=[];
offOrigin=geo.offOrigin;
offDetector=geo.offDetector;
rotDetector=geo.rotDetector;
DSD=geo.DSD;
DSO=geo.DSO;


% TODO : Add options for Stopping criteria
for ii=1:niter

    % If verbose, time the algorithm
    if (ii==1 && verbose==1);tic;end
    % If quality is going to be measured, then we need to save previous image
    % THIS TAKES MEMORY!
    if measurequality
        res_prev=res;
    end


    for jj=1:length(alphablocks)
        % Get offsets
        if size(offOrigin,2)==size(angles,2)
            geo.offOrigin=offOrigin(:,orig_index{jj});
        end
        if size(offDetector,2)==size(angles,2)
            geo.offDetector=offDetector(:,orig_index{jj});
        end
        if size(rotDetector,2)==size(angles,2)
            geo.rotDetector=rotDetector(:,orig_index{jj});
        end
        if size(DSO,2)==size(angles,2)
            geo.DSO=DSO(:,orig_index{jj});
        end
        if size(DSD,2)==size(angles,2)
            geo.DSD=DSD(:,orig_index{jj});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Slower and memory-eating code (but clearer)%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %proj_err=proj(:,:,orig_index{jj})-Ax(res,geo,alphablocks{:,jj},'interpolated'); %                                 (b-Ax)
        %weighted_err=W(:,:,orig_index{jj}).*proj_err;                                 %                          W^-1 * (b-Ax)
        %backprj=Atb(weighted_err,geo,alphablocks{:,jj},'FDK');                          %                     At * W^-1 * (b-Ax)
        %weigth_backprj=bsxfun(@times,1./sum(V(:,:,orig_index{jj}),3),backprj);        %                 V * At * W^-1 * (b-Ax)
        %res=res+lambda*weigth_backprj;                                                % x= x + lambda * V * At * W^-1 * (b-Ax)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nesterov
            % The nesterov update is quite similar to the normal update, it
            % just uses this update, plus part of the last one.
            ynesterov=res +bsxfun(@times,1./sum(V(:,:,jj),3),Atb(W(:,:,orig_index{jj}).*(proj(:,:,orig_index{jj})-Ax(res,geo,alphablocks{:,jj}, 'gpuids', gpuids)),geo,alphablocks{:,jj}, 'gpuids', gpuids));
            res=(1-gamma)*ynesterov+gamma*ynesterov_prev;
        else
            res=res+lambda* bsxfun(@times,1./sum(V(:,:,jj),3),Atb(W(:,:,orig_index{jj}).*(proj(:,:,orig_index{jj})-Ax(res,geo,alphablocks{:,jj}, 'gpuids', gpuids)),geo,alphablocks{:,jj}, 'gpuids', gpuids));
        end

        % Non-negativity constrain
        if nonneg
            res=max(res,0);
        end

    end

    if ii==1
        % Preallocate for frame results
        [ip,iq,ir]=size(res);
        disp(['size(res) is:  ',num2str(ip),' by ',num2str(iq),' by ',num2str(ir)]);
        resultImageOut=zeros(ip,iq,niter,nSlices); % Base size of (resultImageOut) on size of (res)
        clear ip iq ir;
        [ip,iq,ir,is]=size (resultImageOut);
        disp(['size(resultImageOut) is:  ',num2str(ip),' by ',num2str(iq),' by ',num2str(ir),' by ',num2str(is)]);
        clear ip iq ir is;
        tIterationResults=zeros(niter,1);
        n2StopResults=zeros(niter,1);
        itStop=zeros(niter,1);
    end

    %     ===============================================================
    %       Investigate 'stopping' based on Sharpness
    %
    % NB: The stopping is based upon a specific image; 'result_image' in this case.
    % Extract image result_image

    result_Req=res(:,:,imgSlices(1));
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
        resGM=res;
        n2Stop=niter;

        if itFreewheel
            disp('N.B. "Freewheel" iteration to max number of iterations selected.');
        end
        %
        tIterationResults(ii)=round(tTaken);
        tOld=tTaken;
    else 
        resGM=res;

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
        result_Req=res(:,:,imgSlices(iSlice));
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


        % Now plot the progression of AQI throughout the iterations.
        plot(Xval(1:ii),QQ(1:ii),'-r');
        title('AQI Plot');
        xlabel('No. of iterations');
        ylabel('AQI value');

        [ymax,imax]=max(QQ);
        hold on
        plot(Xval(imax),ymax,'ob');
        hold off

        ylim([0 1.1*ymax])
        shg;
        % Save the plot as a .fig file with a default filename for the time being.
        savefig('OS_SARTs7_AQI')
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
        savefig('OS_SARTs7_Sharpness')
        clf('reset')

        % Write iteration timing data to Excel file
        % Times given to nearest second for each iteration
        % Form Excel filename (which includes date & time stamp)
        filename=strcat('tIterationResults',datestr(now,'dd-mm-yy_HH-MM-SS'),'.xlsx');
        % Create table
        T = table((1:niter)',tIterationResults,n2StopResults,itStop,...
            lhSharp,QQ,'VariableNames',...
            {'Iteration','Iteration time','n2Stop','itStop','Sharpness','AQI'});
        % Write table to Excel spreadsheet
        writetable(T,filename);

        break
    end

    %     ===============================================================





    % If quality is being measured
    if measurequality

        %Can save quality measure for every iteration here
        %See if some image quality measure should be used for every
        %iteration?
        qualMeasOut(:,ii)=Measure_Quality(res_prev,res,QualMeasOpts);
    end

    % reduce hyperparameter
    if nesterov
        gamma=(1-lambda);
        lambda=(1+sqrt(1+4*lambda^2))/2;
        gamma=gamma/lambda;
    else
        lambda=lambda*lambdared;
    end

    if computeL2 || nesterov
        % Compute error norm2 of b-Ax
        geo.offOrigin=offOrigin;
        geo.offDetector=offDetector;
        geo.DSD=DSD;
        geo.rotDetector=rotDetector;
        errornow=im3Dnorm(proj-Ax(res,geo,angles,'Siddon', 'gpuids', gpuids),'L2');
        %     If the error is not minimized
        if ii~=1 && errornow>errorL2(end) % This 1.1 is for multigrid, we need to focus to only that case
            if verbose
                disp(['Convergence criteria met, exiting on iteration number:', num2str(ii)]);
            end
            return;
        end
        %     Store Error
        errorL2=[errorL2 errornow];
    end
    % If timing was asked
    if ii==1 && verbose==1
        expected_time=toc*(niter-1);
        expected_duration=toc*(niter);
        disp('OS-SART');
        disp(['Expected duration  :    ',secs2hms(expected_duration)]);
        disp(['Expected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end

end
end



function initres=init_multigrid(proj,geo,alpha)

finalsize=geo.nVoxel;
% start with 64
geo.nVoxel=[64;64;64];
geo.dVoxel=geo.sVoxel./geo.nVoxel;
if any(finalsize<geo.nVoxel)
    initres=zeros(finalsize');
    return;
end
niter=100;
nblock=20;
initres=zeros(geo.nVoxel');
while ~isequal(geo.nVoxel,finalsize)


    % solve subsampled grid
    initres=OS_SART(proj,geo,alpha,niter,'BlockSize',nblock,'Init','image','InitImg',initres,'Verbose',0);

    % Get new dims.
    geo.nVoxel=geo.nVoxel*2;
    geo.nVoxel(geo.nVoxel>finalsize)=finalsize(geo.nVoxel>finalsize);
    geo.dVoxel=geo.sVoxel./geo.nVoxel;
    % Upsample!
    % (hopefully computer has enough memory............)
    [y, x, z]=ndgrid(linspace(1,size(initres,1),geo.nVoxel(1)),...
        linspace(1,size(initres,2),geo.nVoxel(2)),...
        linspace(1,size(initres,3),geo.nVoxel(3)));
    initres=interp3(initres,x,y,z);
    clear x y z
end

end

%% Parse inputs
function [block_size,lambda,res,lambdared,verbose,QualMeasOpts,OrderStrategy,nonneg, gpuids]=parse_inputs(proj,geo,alpha,argin)
opts=     {'blocksize','lambda','init','initimg','verbose','lambda_red','qualmeas','orderstrategy','nonneg', 'gpuids'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('TIGRE:OS_SART:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('TIGRE:OS_SART:InvalidInput',['Optional parameter "' argin{ii} '" does not exist' ]);
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
            error('TIGRE:OS_SART:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]);
        end
        val=argin{jj};
    end

    switch opt
        % % % % % % % Verbose
        case 'verbose'
            if default
                verbose=1;
            else
                verbose=val;
            end
            if ~is2014bOrNewer
                warning('TIGRE: Verbose mode not available for older versions than MATLAB R2014b');
                verbose=false;
            end
            % % % % % % % hyperparameter, LAMBDA
        case 'lambda'
            if default
                lambda=1;
            elseif ischar(val)&&strcmpi(val,'nesterov')
                lambda='nesterov'; %just for lowercase/upercase
            elseif length(val)>1 || ~isnumeric( val)
                error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
            else
                lambda=val;
            end
        case 'lambda_red'
            if default
                lambdared=1;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:OS_SART:InvalidInput','Invalid lambda')
                end
                lambdared=val;
            end
        case 'blocksize'
            if default
                block_size=20;
            else
                if length(val)>1 || ~isnumeric( val)
                    error('TIGRE:OS_SART:InvalidInput','Invalid BlockSize')
                end
                block_size=val;
            end

        case 'init'
            res=[];
            if default || strcmp(val,'none')
                res=zeros(geo.nVoxel','single');
                continue;
            end
            if strcmp(val,'FDK')
                res=FDK(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'multigrid')
                res=init_multigrid(proj,geo,alpha);
                continue;
            end
            if strcmp(val,'image')
                initwithimage=1;
                continue;
            end
            if isempty(res)
                error('TIGRE:OS_SART:InvalidInput','Invalid Init option')
            end
            % % % % % % % ERROR
        case 'initimg'
            if default
                continue;
            end
            if exist('initwithimage','var')
                if isequal(size(val),geo.nVoxel')
                    res=single(val);
                else
                    error('TIGRE:OS_SART:InvalidInput','Invalid image for initialization');
                end
            end
        case 'qualmeas'
            if default
                QualMeasOpts={};
            else
                if iscellstr(val)
                    QualMeasOpts=val;
                else
                    error('TIGRE:OS_SART:InvalidInput','Invalid quality measurement parameters');
                end
            end
        case 'orderstrategy'
            if default
                OrderStrategy='random';
            else
                OrderStrategy=val;
            end
        case 'nonneg'
            if default
                nonneg=true;
            else
                nonneg=val;
            end
        case 'gpuids'
            if default
                gpuids = GpuIds();
            else
                gpuids = val;
            end

        otherwise
            error('TIGRE:OS_SART:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in OS_SART_CBCT()']);
    end
end

end

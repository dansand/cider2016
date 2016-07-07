% Script to calculate G matrix for a surface wave dataset using ray theory

%% BLOCK ONE: Set some basic variables for use in the script
clear all;

%add m_map tools to path
addpath(fullfile(cd,'m_map')); %Make sure to remove this at the end if you don't want to leave this in your path

%model and data setup
blocksize = input('Set an approximate block dimension in degrees = '); 
%blocksize=10.0; %the approximate block dimension in degrees

dataperiod=input('Set a central period of the data in seconds (050, 100, or 150) = ','s'); 

%dataperiod='050'; %the central period of the data in s (050,100, or 150)
%dataperiod='100';
%dataperiod='150';
datafile=['R' dataperiod '.raw1.sel']; %the phase measurement file
cref=input('Set the reference velocity (km/s) appropriate for your choice of central period = '); 

%cref=3.952; %the reference velocity for 50 s data
%cref=4.080; %the reference velocity for 100 s data
%cref=4.280; %the reference velocity for 150 s data

%data and model covariance setup
emult=input('Set a multiplier on data error estimates = '); 
%emult=1; %multiplier for data error estimates 

% Output filenames and force recalculation flags
%Matrix calculation
outputmatrix=['G' dataperiod '.' num2str(blocksize) '.mat']; %output file for G matrix
force_G_calculate=0; %set to 1 if recalculation of matrix is required even if outputmatrix exists, otherwise 0

%some constants
rad=pi/180.0;
fac=2*pi*6371.0/360.0;


% Set up basic data and model geometry
[ nblk,bsize,nlat,mlat,hsize ] = blks2d(blocksize);
fprintf('There are %d blocks in the model\n',nblk);

% get max number of data (some may be skipped)
fid=fopen(datafile);
ndatamax=linecount(fid);
fprintf('There are %d measurements in file %s\n',ndatamax,datafile);

% Check for existence of outputmatrix
if(exist(fullfile(cd,outputmatrix),'file')==2&&force_G_calculate==0)
    fprintf('Reading saved G matrix\n');
    load(outputmatrix);

    %Make a plot of the sampling of the model space
    sampling=1.0;
    mhit=zeros(nblk,1);
    for i=1:nblk
        mhit(i)=numel(find(G_sparse(:,i)));
    end
    [modlat,modlon,hz]=blks_resample(nblk,bsize,nlat,mlat,hsize,mhit,sampling);%sample model on regular grid
    figure;%Pacific centered plot
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    m_pcolor(modlon,modlat,log10(hz)); shading flat;
    plotplates(360);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title({['log_{10} of hit count for period ' dataperiod]});

else
% Set up space for some matrices and vectors
    fprintf('Calculating G matrix\n');
    %G_temp=spalloc(ndatamax,nblk,round(0.05*ndatamax*nblk));  %assumes max 5% nonzero in G (may need to tune for very coarse parameterization or finite frequency)
    G_temp=zeros(ndatamax,nblk); %Can be huge. Use spalloc above if memory is tight.  Will run slightly slower
    d_obs_temp=zeros(ndatamax,1);
    d_err_temp=zeros(ndatamax,1);

% Loop over data file

    fgetl(fid); %first line is junk
    nproc=0;
    nskip=0;

    tic;
    tline=fgetl(fid);
    while ischar(tline)
        %parse out relevant data from line
        A=sscanf(tline,'%7c %8c %f %f %f %f %f %f %f %f %f %f %f %c %f %f %f %f %f %f');
        cmtname=char(A(1:7))';
        stat=char(A(8:12))';
        cmtlat=A(16);
        cmtlon=A(17);
        delc=A(21);
        distance=A(22);
        dphase=A(23);
        period=A(28);
        weight=A(32);
        weight2=A(33);
    
        if(weight<=0||weight2<=0) %skip zero weighted data
            nskip=nskip+1;
            fprintf('skip 0\n');
            tline=fgetl(fid);
            continue
        end
    
        om=2.0*pi/period;
        arg=om*distance*delc/dphase;
        if(arg>=0) %inconsistent data
            nskip=nskip+1;
            fprintf('skip inconsistent\n');
            tline=fgetl(fid);
            continue
        end
        cvel=sqrt(-arg);
        dc=delc/cref;
        errc=cref/(om*distance*weight2);
        norb=1;
        t0=(90.0-cmtlat)*rad;
        p0=cmtlon*rad;
        t0=geocen(t0);
    
        %fix some station name issues
        if(stat(4)=='-'||strcmp(stat(1:3),'PFO'))
            stat(4:5)='  ';
        end
        if(strcmp(stat(1:4),'IPAS'))
            stat='PAS  ';
        end
        if(stat(5)=='-')
            stat(5)=' ';
        end
        [ts,ps,found]=getstninfo(stat);
        if(found==0) %skip if station not found in stations file
            nskip=nskip+1;
            fprintf('skip found\n');
            tline=fgetl(fid);
            continue
        end
    
        %calculate row of G matrix according to distances from ray theory
        [row,delt]=sray(t0,p0,ts,ps,norb,nblk,nlat,bsize,mlat,hsize);
        rowsum=sum(row/delt);
        if(abs(rowsum-1.0)>0.005) %throw out data if total sum of row ~= delt
            nskip=nskip+1;
            fprintf('skip sum\n');
            tline=fgetl(fid);
            continue
        end
        dist=distance/fac;
        nproc=nproc+1;
        if(mod(nproc,100)==0) %output status every 100 measurements
            fmtstring=repmat('%14.6g ',1,7);
            fprintf(['%9d ' fmtstring '\n'],nproc,rowsum,delt,delt*rowsum,dist,dc,errc,cvel);
        end
    
        %Add in row to sparse G matrix
        G_temp(nproc,:)=0.01*row/delt; %conversion for model in percent perturbation
        d_obs_temp(nproc)=dc;
        d_err_temp(nproc)=errc;
    
        %get next line
        tline=fgetl(fid);
    end
    %G_temp d_obs_temp and d_err_temp have extra rows... reduce
    [I,J,S]=find(G_temp);
    G_sparse=sparse(I,J,S);
    d_obs=d_obs_temp(1:nproc);
    d_err=d_err_temp(1:nproc);
    clear G_temp d_obs_temp d_err_temp;

    %Write out G matrix and data to file
    fprintf('Writing out G matrix\n');
    save(outputmatrix,'G_sparse','d_obs','d_err');
    toc;
    
    %Make a plot of the sampling of the model space
    sampling=1.0;
    mhit=zeros(nblk,1);
    for i=1:nblk
        mhit(i)=numel(find(G_sparse(:,i)));
    end
    [modlat,modlon,hz]=blks_resample(nblk,bsize,nlat,mlat,hsize,mhit,sampling);%sample model on regular grid
    figure;%Pacific centered plot
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    m_pcolor(modlon,modlat,log10(hz)); shading flat;
    plotplates(360);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title({['log_{10} of hit count for period ' dataperiod]});

end %finished G matrix calculation


%% BLOCK 2: TV inversion

%varm=1; %a priori model variance - smaller numbers mean a more highly damped model
varm = input('Set an a priori model variance (smaller numbers mean a more highly damped model, try 0.1 to 10) = '); 
smooth = input('Set a smoothing lengthscale (in degrees, try blocksize/2 to 40) = '); 
%smooth=5; 

%Tarantola & Valette inversion
outputmodelTV=['TV' dataperiod '.' num2str(blocksize) '.' num2str(emult) '.' num2str(varm) '.' num2str(smooth) '.mat'];
force_TV_calculate=1; %set to 1 if recalculation of model is required even if outputmodelTV exists, otherwise 0
plot_TV_model=1; %set to 1 to 


if(exist(fullfile(cd,outputmodelTV),'file')==2&&force_TV_calculate==0)
    fprintf('Reading saved TV model\n');
    load(outputmodelTV);
    %redo thresholding here... allows for changing threshold without
    %recalculating svd
    
    
else
    fprintf('Calculating Tarantola and Valette style inversion\n');
    % Set up data and model covariance matrices for T&V style inversion
    fprintf('Setting up covariance matrices\n');
    tic;
   % Cdinv=diag(1./((emult*d_err).^2)); %vector representing diagonal of inverse data cov matrix
    Cdinv_sparse = sparse(1:length(d_err),1:length(d_err),(emult*d_err).^-2); 
    if(smooth==0)
        Cminv=(1/varm)*eye(nblk); %no covariance between model pars (i.e. no smoothing)
    else
        %Cm=sigma^2*exp(-1/2
        [blat,blon]=blks_latlon(nblk,bsize,nlat,mlat,hsize);
        delta=calc_dist(blat,blon);
        Cm=zeros(nblk);
        Cm(delta<3*smooth)=varm*exp(-0.5*(delta(delta<3*smooth).^2/(smooth^2)));
        %invert Cm using threshold svd
        [U,S,V]=svd(Cm);
        lambda=diag(S);
        lambdap=lambda(lambda>0.000001*lambda(1));
        p=size(lambdap,1);
        Spinv=diag(1./lambdap);
        Up=U(:,1:p);
        Vp=V(:,1:p);
        Cminv=Vp*Spinv*Up';
        %Cminv=inv(Cm);
    end
    toc;

    % TV inversion
    fprintf('Inverting for the model\n');
    tic;
    Ginvg=(G_sparse'*Cdinv_sparse*G_sparse+Cminv)\G_sparse'*Cdinv_sparse;
    mestTV=Ginvg*d_obs;
    %mestTV=(G_sparse'*Cdinv*G_sparse+Cminv)\G_sparse'*Cdinv*d_obs;
    toc;
    
    %Calculate resolution matrix
    RmatTV=Ginvg*G_sparse;
    
  
    
    %Write out model to file
    fprintf('Writing out TV model\n');
    %sparsify matrices before saving
    %Cdinv_sparse=sparse(Cdinv);
    %Cminv_sparse=sparse(Cminv);
    
    %Calculate a posteriori model covariance (simple error propagation
    %fprintf('Calculating a posteriori model covariance\n');
    %tic;
    %Cd=inv(Cdinv_sparse);
    %Cmpost=Ginvg*Cd*Ginvg';
    %toc;
    %Alternatively estimate error with multiple realizations of data vector
    fprintf('Calculating error from 100 realizations of data vector\n');
    tic;
    drand=repmat(d_obs,1,100)+sparse(1:length(d_err),1:length(d_err),d_err)*rand(size(d_obs,1),100);
    mrand=Ginvg*drand;
    merrTV=std(mrand,0,2); %take the standard deviation of the 100 realizations of each model par
    toc;
    
    save(outputmodelTV,'mestTV','Cdinv_sparse','Cminv','RmatTV','merrTV');
    %free up some space
    clear Cm Cminv Cdinv U S V Spinv Up Vp lambda lambdap
end
    
    
%% BLOCK 3: Plot up mestTV
if(plot_TV_model==1)
    sampling=1.0;
    [modlat,modlon,mz]=blks_resample(nblk,bsize,nlat,mlat,hsize,mestTV,sampling);%sample model on regular grid
    %plotting using mapping toolbox is commented out
    %[Z,refvec]=geoloc2grid(modlat,modlon,mz,bsize);%convert to format for plotting with mapping toolbox
    figure;%Pacific centered plot
    subplot(2,1,1);
    %worldmap([-90 90],[0 360]);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    %meshm(Z,refvec);
    m_pcolor(modlon,modlat,mz); shading flat;
    plotplates(360);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title({['T and V Model \sigma^2: ' num2str(varm) ' Smoothing: ' num2str(smooth)];'Pacific centered'});
    subplot(2,1,2);%Atlantic centered plot
    %modlon(modlon>180)=modlon(modlon>180)-360.0;
    %worldmap([-90 90],[-180 180]);
    m_proj('Mollweide','lat',[-90 90],'lon',[-180 180]);
    colormap('jet');
    colormap(flipud(colormap));
    %[Z,refvec]=geoloc2grid(modlat,modlon,mz,bsize);
    %meshm(Z,refvec);
    m_pcolor([modlon-360.,modlon],[modlat,modlat],[mz,mz]); shading flat; %trick to handle wraparound
    plotplates(180);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title('Atlantic centered');
end

%% BLOCK: 4 Plot up uncertainty estimates with TV method

%Checkerboard resolution test parameters
chk_l = input('Set the angular order of spheral hamornic for checkerboard pattern (try 2 to 180/blocksize) = '); 
%chk_l=12; %angular order of spherical harmonic based checkerboard, length scale is ~180/chk_l;


if(plot_TV_model==1)
    %Also plot up resolution matrix test and error
    ckmodel=make_sh_checkerboard(chk_l,nblk,bsize,nlat,mlat,hsize);
    ckmodel_out=RmatTV*ckmodel;
    [modlat,modlon,chkinz]=blks_resample(nblk,bsize,nlat,mlat,hsize,ckmodel,sampling);
    [modlat,modlon,chkoutz]=blks_resample(nblk,bsize,nlat,mlat,hsize,ckmodel_out,sampling);
    figure;
    subplot(3,1,1);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    m_pcolor(modlon,modlat,chkinz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title('T and V Resolution test input checkerboard');
    subplot(3,1,2);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    m_pcolor(modlon,modlat,chkoutz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    caxis([-10 10]);
    title('Resolution test output');
    %modelerr=diag(Cmpost);
    [modlat,modlon,errz]=blks_resample(nblk,bsize,nlat,mlat,hsize,merrTV,sampling);
    subplot(3,1,3);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    %colormap(flipud(colormap));
    m_pcolor(modlon,modlat,errz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    caxis('auto');
    title('Error map');
end


%% BLOCK 5: SVD inversion

%SVD inversion threshold parameter
thresh = input('Set a threshold for SVD inversion relative to largest eigenvalue = '); 
%thresh=0.01; %threshold is set relative to largest singular value (thresh*S(1,1))

%SVD based inversion
outputmodelSVD=['SVD' dataperiod '.' num2str(blocksize) '.' num2str(emult) '.mat'];
force_SVD_calculate=0;
plot_SVD_model=1;

if(exist(fullfile(cd,outputmodelSVD),'file')==2&&force_SVD_calculate==0)
    fprintf('Reading saved SVD model\n');
    load(outputmodelSVD);
    d_obs_scale=Cdinv_sparse*d_obs;

    %Calculate threshold
    fprintf('Calculating from saved SVD matrices and threshold = %f\n',thresh);
    p=size(lambda(lambda>thresh*lambda(1)),1);
    fprintf('Keeping %d singular values using a %f threshold\n',p,thresh);
    Spinv=diag(1./lambda(1:p));
    Up=U(:,1:p);
    Vp=V(:,1:p);
    GinvgSVD=Vp*Spinv*Up';
    mestSVD=GinvgSVD*d_obs_scale;
    
    %Resolution matrix
    RmatSVD=Vp*Vp';
    
    %Estimate error with multiple realizations of data vector
    fprintf('Calculating error from 100 realizations of data vector\n');
    tic;
    drand=Cdinv_sparse*(repmat(d_obs,1,100)+sparse(1:length(d_err),1:length(d_err),d_err)*rand(size(d_obs,1),100));
    mrand=GinvgSVD*drand;
    merrSVD=std(mrand,0,2); %take the standard deviation of the 100 realizations of each model par
    toc;

    

else
    %First use inverse data covariance matrix for weighting of both G
    %matrix and data vector (should be calculated or read from TV inversion
    %above
    fprintf('Applying weighting to G matrix and data vector\n');
    G_sparse=Cdinv_sparse*G_sparse;
    d_obs_scale=Cdinv_sparse*d_obs;
    
    %calculate SVD of G matrix
    fprintf('Calculating SVD of G matrix');
    tic;[U,S,V]=svds(G_sparse,nblk);toc; %throws out data null vectors assuming overdetermined problem
    lambda=diag(S);
    %Calculate threshold
    p=size(lambda(lambda>thresh*lambda(1)),1);
    fprintf('Keeping %d singular values using a %f threshold\n',p,thresh);
    Spinv=diag(1./lambda(1:p));
    Up=U(:,1:p);
    Vp=V(:,1:p);
    GinvgSVD=Vp*Spinv*Up';
    mestSVD=GinvgSVD*d_obs_scale;
    
    %Resolution matrix
    RmatSVD=Vp*Vp';
    
    %Estimate error with multiple realizations of data vector
    fprintf('Calculating error from 100 realizations of data vector\n');
    tic;
    drand=Cdinv_sparse*(repmat(d_obs,1,100)+sparse(1:length(d_err),1:length(d_err),d_err)*rand(size(d_obs,1),100));
    mrand=GinvgSVD*drand;
    merrSVD=std(mrand,0,2); %take the standard deviation of the 100 realizations of each model par
    toc;
   
    
    %Write out model to file
    %To allow changing threshold without recalculating SVD, save U, V, and
    %lambda as well.  Throw out all vectors in data null space (assuming
    %overdetermined problem)
   
    fprintf('Writing out SVD model\n');
    save(outputmodelSVD,'mestSVD','V','U','lambda');
end

%% BLOCK 6: Plot up mestSVD
if(plot_SVD_model==1)
    sampling=1.0;
    [modlat,modlon,mz]=blks_resample(nblk,bsize,nlat,mlat,hsize,mestSVD,sampling);%sample model on regular grid
    %plotting using mapping toolbox is commentedt out
    %[Z,refvec]=geoloc2grid(modlat,modlon,mz,bsize);%convert to format for plotting with mapping toolbox
    figure;%Pacific centered plot
    subplot(2,1,1);
    %worldmap([-90 90],[0 360]);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    %meshm(Z,refvec);
    m_pcolor(modlon,modlat,mz); shading flat;
    plotplates(360);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title({['SVD Model threshold: ' num2str(thresh)];'Pacific centered'});
    subplot(2,1,2);%Atlantic centered plot
    %modlon(modlon>180)=modlon(modlon>180)-360.0;
    %worldmap([-90 90],[-180 180]);
    m_proj('Mollweide','lat',[-90 90],'lon',[-180 180]);
    colormap('jet');
    colormap(flipud(colormap));
    %[Z,refvec]=geoloc2grid(modlat,modlon,mz,bsize);
    %meshm(Z,refvec);
    m_pcolor([modlon-360.,modlon],[modlat,modlat],[mz,mz]); shading flat; %trick to handle wraparound
    plotplates(180);
    %plotcoasts;
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title('Atlantic centered');
end

%% BLOCK 7: Plot up mestSVD resolution matrix test

%Checkerboard resolution test parameters
chk_l = input('Set the angular order of spheral hamornic for checkerboard pattern (try 2 to 180/blocksize) = '); 
%chk_l=12; %angular order of spherical harmonic based checkerboard, length scale is ~180/chk_l;

if(plot_SVD_model==1)
    ckmodel=make_sh_checkerboard(chk_l,nblk,bsize,nlat,mlat,hsize);
    ckmodel_out=RmatSVD*ckmodel;
    [modlat,modlon,chkinz]=blks_resample(nblk,bsize,nlat,mlat,hsize,ckmodel,sampling);
    [modlat,modlon,chkoutz]=blks_resample(nblk,bsize,nlat,mlat,hsize,ckmodel_out,sampling);
    figure;
    subplot(3,1,1);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    m_pcolor(modlon,modlat,chkinz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    title('SVD Resolution test input checkerboard');
    subplot(3,1,2);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    colormap(flipud(colormap));
    m_pcolor(modlon,modlat,chkoutz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    caxis([-10 10]);
    title('Resolution test output');
    [modlat,modlon,errz]=blks_resample(nblk,bsize,nlat,mlat,hsize,merrSVD,sampling);
    subplot(3,1,3);
    m_proj('Mollweide','lat',[-90 90],'lon',[0 360]);
    colormap('jet');
    %colormap(flipud(colormap));
    m_pcolor(modlon,modlat,errz); shading flat;
    plotplates(360);
    m_coast('color','black');
    m_grid('xaxislocation','middle');
    colorbar;
    caxis('auto');
    title('Error map');

end

%% BLOCK 8: Clean up

%Remove m_map tools from path
rmpath(fullfile(cd,'m_map'));
    


%% ACM 154 project
% Goal: use Kalman Filter for a better phase estimate in phase retrieval
%
% Kalman filter: line 221 - line 276
%
% This code is developed based on the paper/code by Wu et. al: 
% Tengfei Wu, Ori Katz, Xiaopeng Shao, and Sylvain Gigan, Opt. Lett. 41, 5003-5006 (2016)
% Source code: https://dx.doi.org/10.6084/m9.figshare.3593370.v1 
%
% Bispecturm is used for the initial (1D) phase estimatiom. Kalman filter
% is applied to the unwrapped (1D) phase estimate.
%
% Toolbox required:
% Image Processing Toolbox, Statistics and Machine Learning Toolbox, Curve Fitting Toolbox
%% parameters
showComparison  = true;
autoPatchSize   = false;        % whether to choose the patch size and sigma automatically 

projection_theta    = 0:10:179; % angles used in radon transformation
overlap_ratio       = 0.85;     % overlap ratio between two adjacent patches

if ~autoPatchSize
    envelopeDimension = 240;    % match up to the original code
    sigma             = 50;     % match up to the original code
    overlap_ratio     = 0.85;   % match up to the original code
end

%% Load experimental data
Large_Speckle_intensity = double(imread('Digit 4 experimental data.tif'));
Large_Speckle_intensity = flip(Large_Speckle_intensity,1);

projectionNum = length(projection_theta);
[xsize,ysize] = size(Large_Speckle_intensity);
xc            = floor(xsize/2+1); % origin for FT
yc            = floor(ysize/2+1); % origin for FT
[Y,X]         = meshgrid(1:ysize,1:xsize);
R             = abs((X-xc)+1i*(Y-yc));
%% calculate autocorrelation
LowPass_Filter  = zeros(xsize,ysize);
Filter_Length   = 15;
LowPass_Filter((xc-Filter_Length):(xc+Filter_Length),...
               (yc-Filter_Length):(yc+Filter_Length)) = 1;
Large_Speckle_intensity_FFT = fftshift(fft2(Large_Speckle_intensity));
Large_Speckle_FFT_lowpass = Large_Speckle_intensity_FFT.*LowPass_Filter;
Large_Speckle_intensity_LowPassVer = ifft2(ifftshift(Large_Speckle_FFT_lowpass));

speckle_filteredFT    = abs(fft2(fftshift(Large_Speckle_intensity./Large_Speckle_intensity_LowPassVer))).^2;
Speckle_intensity_large_AutoCorr = ifftshift(ifft2((speckle_filteredFT)));
Speckle_intensity_large_AutoCorr = Speckle_intensity_large_AutoCorr - min(Speckle_intensity_large_AutoCorr(:));
Speckle_intensity_large_AutoCorr = Speckle_intensity_large_AutoCorr/Speckle_intensity_large_AutoCorr(xc,yc); %% normalization

clear Large_Speckle_FFT_lowpass Large_Speckle_intensity_LowPassVer speckle_filteredFT LowPass_Filter

%% estimate the object size from autocorrelation
lineProfile     = Speckle_intensity_large_AutoCorr(xc,:);
threshold_AC    = 1.5*mean(lineProfile);
[y_start]       = find(lineProfile>threshold_AC,1);
area2cover      = 5*abs(y_start-yc); 

% image (cropped) to be used in calculating the radial profile
autoCorr2use = Speckle_intensity_large_AutoCorr(xc-area2cover:xc+area2cover,...
                                               yc-area2cover:yc+area2cover);
radialProfile = radialAverage(autoCorr2use,area2cover+1,area2cover+1,1:area2cover);
objSize = find(radialProfile < min(radialProfile)*1.5,1); % object size estimate

clear Speckle_intensity_large_AutoCorr
%% estimate amplitude of object's FT
virtualWindow = R< (objSize/2);
% the virtual window is used for estimate the amplitudes in low spatial frequencies 

ftWindow    = fftshift(fft2(virtualWindow));
lineProfile = abs(ftWindow(xc,yc:end));
width       = find(lineProfile./lineProfile(1)<0.3,1);
img2fit     = abs(Large_Speckle_intensity_FFT(xc-width:xc+width,yc-width:yc+width));
[Y,X]       = meshgrid(-width:width);
Rtemp       = abs(X + 1i*Y);
mask        = Rtemp>0.7*width;
X           = X.*mask;
Y           = Y.*mask;

% vectorize X, Y, and img2fit
X       = X(:);
Y       = Y(:);
ydata   = img2fit(:);
ydata   = ydata(X~=0);
X       = X(X~=0);
Y       = Y(Y~=0);

% sample from X, Y, and img2fit
idx     = randsample(1:length(X),2000);
X       = X(idx);
Y       = Y(idx);
ydata   = ydata(idx);

lineProf2 = abs(ftWindow(xc-width:xc+width,yc));
f1 = fit( (-width:width).',lineProf2,'gauss1');
lineProf2 = abs(ftWindow(xc,yc-width:yc+width));
f2 = fit( (-width:width).',lineProf2.','gauss1');

covAmp0= [f1.c1^2,0,f2.c1^2,mean(ydata)*2];
covAmpFit = lsqcurvefit(@Gaussian2D,covAmp0,[X,Y].',ydata.');

% replace the amplitude for low spatial frequency
[Y,X]   = meshgrid(-width:width);
bestFit = Gaussian2D(covAmpFit,[X(:),Y(:)].');
bestFit = reshape(bestFit,size(Y));

mask        = img2fit<mean(ydata)*7;
ampLowFreq  = img2fit.*mask;
normF = max(ampLowFreq(:));

% As DC is the dominate part, amplitude of the zero frequency is the biggest
ampLowFreq = ampLowFreq + (1-mask).*bestFit*max([normF/bestFit(width+1,width+1),1]);

% amplitude estimation
AmplitudeEst = abs(Large_Speckle_intensity_FFT);
AmplitudeEst(xc-width:xc+width,yc-width:yc+width) = ampLowFreq;

clear img2fit img2fit2 lineProfile lineProf2 ftWindow
%% estimate phase of object's FT

%% calculate the patch size

if autoPatchSize
    envelopeDimension = 10*objSize;   % related to the size of the object

    % Note: selection of parameters "envelopeDimension" and "sigma" is very important.
    % Except for the principle metioned above, the parameters should make the 
    % envelope in one frame as large as possible, to ensure the Fourier transform
    % of the Gaussian window function ~delta function (small enough).
    % Too small envope would smoothen the bispectrum phase information (lose
    % phase information of object).

    sigma   = 1.8*objSize;
    % sigma~(2D-3D) D is size of object, larger sigma ensure more information, but introduce more computational complexity  
end
x           = (-(floor(envelopeDimension/2)-1):(floor(envelopeDimension/2)));
[X,Y]       = meshgrid(x);
envelope    = exp(-(X.^2+Y.^2)/(sigma^2));   
interval        = envelopeDimension;
overlap_size    = interval*overlap_ratio;

num_sqrt_1 = floor(1+(xsize-interval)/((1-overlap_ratio)*interval));
num_sqrt_2 = floor(1+(ysize-interval)/((1-overlap_ratio)*interval));
num = num_sqrt_1*num_sqrt_2;

subspeckle = cell(num_sqrt_1,num_sqrt_2);

% patches to use
for i = 1:num_sqrt_1
    for j = 1:num_sqrt_2
        subspeckle{i,j} = (Large_Speckle_intensity((uint16((i-1)*(1-overlap_ratio)*interval+1):uint16(interval*i-overlap_ratio*(i-1)*interval)),...
            (uint16((j-1)*(1-overlap_ratio)*interval+1):uint16(interval*j-overlap_ratio*(j-1)*interval))));
    end
end
subspeckle = subspeckle(:);

clear Large_Speckle_intensity Large_Speckle_intensity_FFT

for q = 1:num

    Speckle_intensity = subspeckle{q};
    Speckle_intensity = Speckle_intensity.*envelope;       
    % filter in the spatial domain with Gaussian window, and its kernel
    % should be between 2D and 3D, where D is the size of object

    R_Speckle = radon(Speckle_intensity,projection_theta);
    [projectionLength_Speckle,~] = size(R_Speckle);
    
    if q == 1
    RealLengthCurrentSpeckle = projectionLength_Speckle-1;         
    R_choose_Speckle = 1:RealLengthCurrentSpeckle;
    end
    
    R_Speckle = R_Speckle(R_choose_Speckle,:);
    if q == 1
        R_Speckle_Store = zeros(RealLengthCurrentSpeckle,length(projection_theta),num);
        R_Speckle_Store(:,:,q) = R_Speckle;
    else
        R_Speckle_Store(:,:,q) = R_Speckle;
    end
end

%% start phase estimation
objPhase_proj = zeros(RealLengthCurrentSpeckle,projectionNum);

nfft    = RealLengthCurrentSpeckle;
Wind    = 1;
nsamp   = RealLengthCurrentSpeckle;
overlap = 0;
mask    = hankel(1:nfft,[nfft,1:nfft-1]);

for projection = 1:projectionNum
    
    fprintf('Start of No.%d projection, %d projections in total.\n\n',projection,projectionNum);

    if nfft<nsamp
        nfft = 2^nextpow2(nsamp);
    end
    Bspec_Speckle = zeros(nfft,nfft);
    Intialization = zeros(nfft,nfft);
    
    amplitude_esti = zeros(1,RealLengthCurrentSpeckle);
    
    for Frame = 1:num
        
%         [Bspec_Speckle_1,~] = Calbispectrum2D(R_Speckle_Store(:,projection,Frame).',nfft,Wind,nsamp,overlap,Intialization);

        [Bspec_Speckle_1,ftY] = Calbispectrum2DV2(R_Speckle_Store(:,projection,Frame),nfft,mask);
        Bspec_Speckle = Bspec_Speckle + Bspec_Speckle_1;                                                                                                                % strong background noise    
    end

    Bspec_Speckle = Bspec_Speckle/num;
    Bspec_Speckle = rot90(Bspec_Speckle);
    Bspec_Speckle = circshift(Bspec_Speckle',[0,1]);

    objPhase_proj(:,projection) = RecursiveProcessV2(Bspec_Speckle,RealLengthCurrentSpeckle,ftY);

end

%%   %%%%%%%%%%%%%% use kalman filter to estimate the phase %%%%%%%%%%%%%%
% prior, p_0~Gaussian(0,0.001); (phase of the zero freq is 0) C0 = 0

% p_(k) -> p_(k+1)+xi:
% p_(k+1) = M*p_k, where  M = 1, xi~N(0,sigma) in our case
% and sigma can be estimated using the object size.

% Why is xi Gaussian?
% This be derived using random phasor sum + a constant phasor. When the
% constant phasor is large, xi is approximately Gaussian.

% y_(k+1) = H*p_(k+1) + yita, where H = 1, yita ~ Norm(0, gamma)

% Why is yita Gaussian?
% Possion noise can be approximated as Gaussian when the signal is large.

% measurement y_k+1 = estimatedPhase_(k+1)

% predict correlation from object size
ratioFT = envelopeDimension/(xsize+ysize/2);
% why we need this factor: intervals in the grid of the DFT for different
% images (different sizes) correspond to different spatial frequencies

windowInFT = Rtemp< (width*ratioFT);
temp = circshift(windowInFT,[1,0]);
indNum = sum(sum(windowInFT - temp>0));

sigma   = (indNum/(sum(windowInFT(:)) - 2*indNum))^2;
gamma   = 0.2;
M       = 1;
H       = 1;
C0      = 0;

phaseMeasurementSize = size(Bspec_Speckle,1);
pc          = floor(phaseMeasurementSize/2+1); % coordinate of origin in the FT of the phaseMeasurement
J    = phaseMeasurementSize - pc + 1;          % largest time index 

objPhase_projEst = zeros(size(objPhase_proj));

for projection = 1:projectionNum
    % C: variance, phaseKF: phase estimates in KF
    C       = zeros(1,J);
    phaseKF = C;
    C(1)    = C0;
    y       = unwrap(objPhase_proj(pc:end,projection));
    
    for idx = 1:J-1
        [phaseKF(idx+1),C(idx+1)] = KF_phase(phaseKF(idx),y(idx+1),sigma,gamma,H,M,C(idx));
    end
    objPhase_projEst(pc:end,projection) = phaseKF;
    if mod(phaseMeasurementSize,2) == 0  % use property of FT of a real signal
        objPhase_projEst(2:pc-1,projection) = -flip(phaseKF(2:end));
    else
        objPhase_projEst(1:pc-1,projection) = -flip(phaseKF(2:end));
    end
end

%% reconstruction
Object_ReSpeckleProjFFT = zeros(RealLengthCurrentSpeckle,projectionNum);
for projection = 1:projectionNum
    Object_ReSpeckleProjFFT(:,projection) = exp(1i*objPhase_projEst(:,projection));
end

% transform the polar coordinate to the Cartesian coordinate
phase = PolarToCartesian(RealLengthCurrentSpeckle,projection_theta,Object_ReSpeckleProjFFT);
if mod(size(phase,2),2) == 0
    phase = circshift(phase,[0,1]);
end

% amplitude
amp = imresize(AmplitudeEst,size(phase),'box'); % medfilt2(AmplitudeEst,[10,10])
if mod(size(amp,1),2) == 0
    amp = amp + circshift(flip(flip(amp,1),2),[1,1]);
else
    amp = amp + flip(flip(amp,1),2);
end
recoveredFT     = amp.*exp(1i*phase);
recoveredImgKF  = abs(ifftshift(ifft2(recoveredFT))); % image reconsturction with KF

if ~showComparison
    figure;imagesc(recoveredImgKF);
end

%% for comparison
if showComparison
    Object_ReSpeckleProjFFT_Org = zeros(RealLengthCurrentSpeckle,projectionNum);
    for projection = 1:projectionNum
        Object_ReSpeckleProjFFT_Org(:,projection) = exp(1i*objPhase_proj(:,projection));
    end

    % transform the polar coordinate to the Cartesian coordinate
    phaseOrg = PolarToCartesian(RealLengthCurrentSpeckle,projection_theta,Object_ReSpeckleProjFFT_Org);
    if mod(size(phaseOrg,2),2) == 0
        phaseOrg = circshift(phaseOrg,[0,1]);
    end
    recoveredFTOrg  = amp.*exp(1i*phaseOrg);
    recoveredImgOrg = abs(ifftshift(ifft2(recoveredFTOrg)));
    
    [xsize,ysize] = size(recoveredImgOrg);
    xc            = floor(xsize/2+1); % origin for FT
    yc            = floor(ysize/2+1); % origin for FT
    myhot = hot(256);
    figure;colormap(myhot(10:246,:));
    subplot(2,2,1);imagesc(recoveredImgKF);axis image;title('Reconstructed Image, w/ KF');
    axis([xc-3*objSize,xc+3*objSize,yc-3*objSize,yc+3*objSize]);
    subplot(2,2,2);imagesc(recoveredImgOrg);axis image;title('Reconstructed Image, wo KF');
    axis([xc-3*objSize,xc+3*objSize,yc-3*objSize,yc+3*objSize]);
    subplot(2,2,3);imagesc(phase);axis image;title('Reconstructed Phase, w/ KF');
    subplot(2,2,4);imagesc(phaseOrg);axis image;title('Reconstructed Phase, wo KF');
end


%% local functions
function profile = radialAverage(IMG, cx, cy, w)
    % computes the radial average of the image IMG around the cx,cy point
    % w is the vector of radii starting from zero
    [a,b] = size(IMG);
    [Y, X] = meshgrid( (1:b)-cy, (1:a)-cx);
    R = abs(X + 1i*Y);
    profile = [];
    for i = w % radius of the circle
        mask = (i-1<=R & R<i+1); % smooth 1 px around the radius
        values = IMG(mask); % without smooth
        profile(end+1) = mean( values(:) );
    end
end



function y = Gaussian2D(CovAmp, x)
    
    myCov = [CovAmp(1),CovAmp(2);CovAmp(2),CovAmp(3)];
    nData = size(x,2);
    y = zeros(1,nData);
%     invCov = inv(myCov);
    for idx = 1:nData
        y(idx) = CovAmp(4)*exp(-1/2*x(:,idx).'/myCov*x(:,idx) );
    end

end

% Kalman Filter for estimating the phase
function [p_kp,C_kp] = KF_phase(p_k,y_kp,sigma,gamma,H,M,C_k)
    temp = inv(sigma + M*C_k*M.');
    C_kp = inv(temp + H/gamma*H.');
    p_kp = C_kp*(temp*M*p_k + H.'/gamma*y_kp);
end
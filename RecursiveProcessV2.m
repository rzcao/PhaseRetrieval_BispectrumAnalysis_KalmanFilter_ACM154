function Phi = RecursiveProcessV2(Bispectrum,N,ftY)
% This is a modified version from the original function.
% This function is used to recover the 1D-object phase from its bispectrum
% with the recursie method.The iterative equation is as follows:
%         ¦×(r) = ¦×(q)+¦×(r-q)-¦Â(r-q,q)
% The main parameters are listed below:

% Bispectrum: the bispectrum of object calculated by HOSA toolbox

% N = size(Bispectrum,1);      % generally the same size as the object,
                               % even number in simulations
% ftY: the Fourier transform of the signal Y, which is used for calculating
% the input "Bispecturm"

Phi = zeros(1,N);
% Phi(1) = 0;

Num         = ceil(N/2);
PhiHalf     = zeros(1,Num);
PhiHalf(1)  = 0;
betaCenter  = floor(N/2+1);
Beta        = angle(Bispectrum(1:betaCenter-1,betaCenter+1:end));
Beta        = flip(Beta,1);
[xsize,ysize] = size(Beta);

useAmp = false;
if nargin == 3
    useAmp      = true;
    mask        = abs(ftY);
    mask(1:5)   = mean(mask(10:25));
end

for r = 2:Num
    idx2use     = (1:floor(r/2))';
    beta2use    = sub2ind([xsize,ysize] , idx2use,(r-idx2use));
    Temp        = PhiHalf(idx2use) + PhiHalf(r - idx2use) - Beta(beta2use).';
    
    if useAmp
        w           = 1./mask(idx2use) + 1./mask(r - idx2use); % weights
        w           = sqrt(1./w);
        PhiHalf(r)  = angle(sum(exp(1i*Temp).*w.')); % weighted "average"
    else
        PhiHalf(r)  = angle(sum(exp(1i*Temp)));
    end
end

Phi(betaCenter:end) = PhiHalf;
if mod(N,2) == 0
    Phi(2:betaCenter-1) = -flip(PhiHalf(2:end));
else
    Phi(1:betaCenter-1) = -flip(PhiHalf(2:end));
end

        
end


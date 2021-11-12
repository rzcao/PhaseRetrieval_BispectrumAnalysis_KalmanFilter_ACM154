function [Bspec,Xf] = Calbispectrum2DV2(y, nfft, mask)
Xf     = fft(y,nfft);    % y-mean(y)
CXf    = conj(Xf);
Bspec  = (Xf * Xf.') .*reshape(CXf(mask), nfft, nfft);
Bspec  = fftshift(Bspec);
end
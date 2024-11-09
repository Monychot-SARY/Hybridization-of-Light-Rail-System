
function [S,k,fssdim,freelle,S2,k2,fece]=spectrece(xx,fece)

% [S,k]=spectre(ce,fe);
%
% Matlab function to estimate the spectrum of a function xx
% 
% Inputs :
%  - xx  : sampled function [N x 1]
%  - fece : sampling frequency (in Hertz)
%
% Outputs :
%  - S : spectrum value    [Nfft x 1]
%  - k : Sampling elements of the spectrum [Nfft x 1]


xx=xx(:);
N=length(xx);
g=hamming(N);% Hamming window

xx=xx.*g;
Nfft=8*length(xx);  
TFD=fft(xx,Nfft); 
S=abs(TFD);      

%% Calibration of the abscissa
S2=S(1:Nfft/2+1);
k=0:Nfft-1; 
k2=k(1:Nfft/2+1);
fssdim=k2/Nfft; 
freelle=k2*fece/Nfft; 

S2=2*S2/abs(sum(g));

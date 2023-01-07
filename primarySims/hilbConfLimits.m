function [confLimits] = hilbConfLimits(data,filtData, Fs,filtBandwidth,alpha)
% this function takes a filtered data 'data' and then genereates
% confidence limits for this using Lepage Kramer Eden results on the 
% asymptotic properties of the sampling distribution for the phase
% under a hilbert transform
% INPUTS
% data: filtered data 
% Fs: sampling frequency
% filtBandwidth: the difference between stopband edges of narrowband filter
% alpha: the probability for confidence limits
% OUTPUTS
% confLimits: 2xN vector with upper and lower bounds on phase

if isempty(alpha)
    alpha = 0.05;
end

phaseMLE = angle(hilbert(filtData));
ampSquare = abs(hilbert(filtData)).^2;

N = length(filtData);
p = N*filtBandwidth*(1/Fs);

dataVar = mean((data - filtData).^2);

phaseSTD = sqrt((p*dataVar)./(2*N*ampSquare));
% under gaussian assumptions then:
mult = norminv(1-alpha/2,0,1);
confLimits(:,1) = phaseMLE - mult*phaseSTD;
confLimits(:,2) = phaseMLE + mult*phaseSTD;
% ensuring max width is 2*pi
inds = find(mult*phaseSTD>pi);
confLimits(inds,1) = -pi;
confLimits(inds,2) = pi;


function [noisim,nois] =  addnoise_outall(im,avg,stddev)% 161105 
% only add noise to real if image is only real
%151203
% add gaussian white noise with mean "avg" and deviation "stddev" (or root
% of variance) to each component of complex data
% outputs noisy image "noisim" and noise alone "nois"

[numrow,numcol,numimage] = size(im);
imim = imag(im); reim = real(im);
for n = 1:numimage
    nsre(:,:,n) = stddev*randn(numrow,numcol);
    ns_reim(:,:,n) = reim(:,:,n) + nsre(:,:,n) +avg;
    if max(imim(:)) > 0 %if there's no imaginary image
        nsim(:,:,n) = stddev*randn(numrow,numcol);
        ns_imim(:,:,n) = imim(:,:,n) + nsim(:,:,n) +avg;
    else ns_imim = 0;nsim = 0;
    end
end
noisim = ns_reim + 1i*ns_imim;
nois = nsre + 1i*nsim;


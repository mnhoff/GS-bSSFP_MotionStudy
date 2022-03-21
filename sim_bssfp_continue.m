function [amat,bmat,mmat,theta,cdata,E1mat,rotM,noisonly,puredata] = sim_bssfp_continue(cyc,noisval,th_pts,flip,nul,highres)
%140727 add noise analysis 
if nargin == 4,nul = 1;end % can nul the M rotation/TE dephasing
npi = 1;%1*pi for -pi to pi, 2*pi for -2pi to 2pi
%% base val vectors: th, a, b, and mag
th = linspace(-npi*pi,npi*pi*(1-2/th_pts),floor(th_pts));
a = [0.91:0.01:0.99 0.9999];
% a = [0.91:0.001:0.999 0.9999];
% a = [0.91:0.0005:0.999 0.9999];
% a = [0.91:0.0001:0.999 0.9999];%best for row calculation plots, but slow
if highres, b = (0.795:-0.0005:0.10)'; else, b = (0.795:-0.005:0.10)'; end % use higher density for actual, lower for trial data%  
% a = [0.1:0.1:0.9 0.9999];
% b = (0.8:-0.02:0.10)';
mag = 4000;
%% matrices: theta, amat, bmat, and mmat
npc = length(cyc);
nr = length(b);
nc = length(a)*th_pts;
theta = repmat(th,nr,length(a));
amat = zeros(nr,nc);
for t = 1:length(a)
    amat(:,th_pts*(t-1)+1:t*th_pts) = repmat(a(t),nr,th_pts);
end
bmat = repmat(b,1,nc);
% don't know T1, T2, so determine M from a, b, and the flip angle
term = amat*(1+cos(flip));
mmat = mag*bmat.*sin(flip)./term;
E1mat = (term-bmat.*(1+amat.^2.*cos(flip)))./(term-bmat.*(amat.^2+cos(flip)));
%% build the datasets
cdata = zeros(nr,nc,npc);
pc = cyc*2*pi;%phase cycles
rotM = mmat.*exp(nul*1i.*theta/2);
for j = 1:npc
    cdata(:,:,j) = rotM.*(1-amat.*exp(-1i.*(theta+pc(j))))./(1-bmat.*cos(theta+pc(j)));
end
puredata = cdata;
%% add gaussian white noise
if ~isempty(noisval)
    [cdata,noisonly] = addnoise_outall(cdata,0,noisval);%data, mean, noise factor
end      

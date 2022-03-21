function [tempa,tempb,tempm,theta,cdata,tempE,tempr,noisonly,motdata,puredata,row1,rowlst,col1,collst,masque] = sim_bssfp_SpecNoisMot(cyc,noisval,th_pts,flip,nul,agran,TR,motdir,freq,ampl,PEdirn)
%171105 simulate bSSFP with motion at a particular noise level, ampl, freq 
% (to be run for many noise vals)
%161031 rev3 switch to San's simulation method
% 160824 rev2 adds image replicas at varied shifts to a 4-cycle bssfp simulation 
% to simulate motion
% amends sim_bssfp_continuous_noisout_agran.m

if nargin == 4,nul = 1;end % can nul the M rotation/TE dephasing
npi = 2;%1*pi for -pi to pi, 2*pi for -2pi to 2pi
%% base val vectors: th, a, b, and mag
% Create a factor to add more th_pts per 2pi if very few a values are chosen
agrnTHfact = 1;
% if agran == 1, agrnTHfact = 3;
% elseif agran == 2, agrnTHfact = 2;
% elseif agran == 3, agrnTHfact = 1.75;
% elseif agran == 4, agrnTHfact = 1.5;
% elseif agran == 5, agrnTHfact = 1.25;
if agran == 1, agrnTHfact = 3;
elseif agran == 2, agrnTHfact = 2;
elseif agran == 3, agrnTHfact = 1.75;
elseif agran == 4, agrnTHfact = 1.5;
elseif agran == 5, agrnTHfact = 1.25;
elseif agran == 9, agrnTHfact = 2;
end
th_pts = agrnTHfact*th_pts;
th = linspace(-npi*pi,npi*pi*(1-2/th_pts),floor(th_pts));
a1 = 0.91;%1st 'a' value, standard
% a1 = 0.98; %focus on water, csf- for the GAS paper accumulated noise plot
al = 0.999;%last 'a' value
resol = (al-a1)/agran; % agran = 2 gives minimal cycles of 2pi
if agran == 9
    a = [a1 0.9999];%special case if I only want 2 a vals
else
    a = [a1:resol:al 0.9999];
end
b = (0.795:-0.005:0.10)';%standard
% b = (0.795:-0.001:0.65)';%accentuate the assymetry of the complex noise at
%higher b values - for the GAS paper accumulated noise plot (water-CSF
%region)
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
%% add a background border and some stucture to accentuate artifacts
%shapes is a dataset that will work to mask the compdata such that data
%only exists within the shapes
shapes = imbinarize(imresize(double(imread('basic_shapes.png')),0.45),0.105);
shapes(:,69:78) = [];
[mnr,mnc] = size(shapes);
row1 = floor((256-nr)/2+1); col1 = floor((256-nc)/2+1);
rowlst = row1+nr-1;  collst = col1+nc-1;
maskrow1 = floor((256-mnr)/2+1);
maskcol1 = floor((256-mnc)/2+1);
maskrowlst = maskrow1+mnr-1;  maskcollst = maskcol1+mnc-1;
zeromatrix = zeros(256,256);temp1 = repmat(zeromatrix,[1 1 4]);masque = zeromatrix;

% set aside compdata with correct dimensions
temp1(row1:rowlst,col1:collst,:) = cdata; clear cdata; cdata = temp1;
masque(maskrow1:maskrowlst,maskcol1:maskcollst) = shapes; %shape-mask

%perform masking
for nn = 1:npc
    cdata(:,:,nn) =  temp1(:,:,nn).*masque;
end                                                       

% regenerate parameter maps
tempm = zeromatrix;tempa = zeromatrix;tempb = zeromatrix;tempr = zeromatrix;tempE = zeromatrix;
tempm(row1:rowlst,col1:collst) = mmat;
tempa(row1:rowlst,col1:collst) = amat;
tempb(row1:rowlst,col1:collst) = bmat;
tempr(row1:rowlst,col1:collst) = rotM;
tempE(row1:rowlst,col1:collst) = E1mat;

%mask new parameter maps
tempm = tempm.*masque;tempa = tempa.*masque;tempb = tempb.*masque;tempr = tempr.*masque;tempE = tempE.*masque;

%% run a motion artifact generator
TR = TR/1000; %switch units to [s]
if strcmpi(motdir,'h')
    [motdata,motKSP] = transMotGen(PEdirn,cdata(:,:,1), TR, 'sinXampl', ampl, 'sinXfreq', freq );
%     disp('Motion is horizontal');
else %default is PE motion
    [motdata,motKSP] = transMotGen(PEdirn,cdata(:,:,1), TR, 'sinYampl', ampl, 'sinYfreq', freq );
%     disp('Motion is vertical');
end

% [PSF, ksp] = simCartesianMRI( 'TR', TR, 'matrix', [nc nr], 'sinYampl', ampl, 'sinYfreq', freq );

%% add gaussian white noise
puredata = cdata;
if ~isempty(noisval)
    [puredata,noisonly] = addnoise_outall(puredata,0,noisval);%data, mean, noise factor
    cdata(:,:,1) = motdata;
    cdata = cdata + noisonly;
end


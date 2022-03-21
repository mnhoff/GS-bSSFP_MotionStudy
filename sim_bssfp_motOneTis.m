function [tempa,tempb,tempm,theta,cdata,tempE,tempr,noisonly,motdata,puredata,row1,rowlst,col1,collst,masque] = sim_bssfp_motOneTis(cyc,noisval,T1,T2,th_pts,flip,nul,TR,motdir,freq,ampl,PEdirn)
%171105 simulate bSSFP with motion for a particular tissue, ampl, freq 
% (to be run for many tissues)
% 161031 rev3 switch to San's simulation method
% 160824 rev2 adds image replicas at varied shifts to a 4-cycle bssfp simulation 
% to simulate motion
% amends sim_bssfp_continuous_noisout_agran.m

if nargin == 4,nul = 1;end % can nul the M rotation/TE dephasing
npi = 2;%1*pi for -pi to pi, 2*pi for -2pi to 2pi
%% base val vectors: th, a, b, and mag
th_pts = 2*th_pts;% this is a granularity decision usually arrived at by input agran
th = linspace(-npi*pi,npi*pi*(1-2/th_pts),floor(th_pts));
%% parameters
E2 = exp(-TR/T2);% = E2
E1 = exp(-TR/T1);
mag0 = 4000;
Q = (1-E1)/(1-E1*cos(flip)-E2^2*(E1-cos(flip)));
b = E2*(1+cos(flip))*Q;
mag = mag0*sin(flip)*Q;
npc = length(cyc);
term = E2*(1+cos(flip));
mm = mag*b*sin(flip)/term;
E1mat = (term-b*(1+E2^2*cos(flip)))/(term-b*(E2^2+cos(flip)));

%% now buld the datasets
nr = 256; nc = 256; threp = nc/length(th);
bmat = repmat(b,nr,nc);
amat = repmat(E2,nr,nc);
mmat = repmat(mm,nr,nc);
theta = repmat(th,nr,threp);
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

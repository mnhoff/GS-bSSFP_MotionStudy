function [tempa,tempb,tempm,theta,cdata,tempE,tempr,noisonly,nonmotdata,motdata,row1,rowlst,col1,collst,freqs,amplitudes,masque] = sim_bssfp_shiftshape(cyc,noisval,th_pts,flip,nul,agran,TR,motdir,PEdirn)
%170731 rev4 vary frequency instead of range (shft4)
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
% cdata = rot90(cdata);% this rotation should make a PE motion artifact issue along y coherently

%% run a motion artifact generator
% freqs  = 40;    amplitudes = 3;%dfrequency and amplitude used in fig2
freqs  = [5 10 25 50];    amplitudes = [1 3 15 25];%debug res
% freqs = [1 5:5:60]; amplitudes = [0 1 2 3:2:27 30];% 30x14med res
% freqs = 1:2:60; amplitudes = 0:30;% 30x31high res

nfreq = length(freqs);
namp = length(amplitudes);
[nr,nc,~] = size(cdata);
motdata = zeros(nr,nc,nfreq,namp);
% for aa = 1:
TR = TR/1000; %switch units to [s]
% ampl  =  30;  % amplitude in Px
cntf = 0;
%  if strcmpi(motdir,'F'), disp('motion is horizontal, as is FE')
%  else disp('Motion is vertical, as is PE'); end
for fr = freqs
    cntf = cntf+1;
    cnta = 0;
    for amp = amplitudes
        cnta = cnta+1;
        if strcmpi(motdir,'h')
            [motdata(:,:,cntf,cnta),motKSP] = transMotGen(PEdirn,cdata(:,:,1), TR, 'sinXampl', amp, 'sinXfreq', fr );
        else %default is PE motion
            [motdata(:,:,cntf,cnta),motKSP] = transMotGen(PEdirn,cdata(:,:,1), TR, 'sinYampl', amp, 'sinYfreq', fr );
        end
    end
end
%print out a few figures at disparate vals
% selF = [10 50];    selA = [3 25];
% % selF = [5 55];    selA = [3 28];%for 30x31 ampsxfreqs
% % selF = [20 50];    selA = [3 20];
% for ff = 1:length(selF)
%     f_ind = find(freqs == selF(ff));
%     for aa = 1:length(selA)
%         a_ind = find(amplitudes == selA(aa));
%         figure; imagesc(abs(motdata(:,:,f_ind,a_ind)));colormap(gray); axis image
%         title([num2str(selF(ff)),'-Cycled ',num2str(selA(aa)),'-Pixel Motion Noiseless image']);
%         motstr = [num2str(selF(ff)),'Cyc_',num2str(selA(aa)),'Pix_' motdir 'EdirnMotNoiselessIm'];
%     end
% end

% [PSF, ksp] = simCartesianMRI( 'TR', TR, 'matrix', [nc nr], 'sinYampl', ampl, 'sinYfreq', freq );

%% add gaussian white noise
if ~isempty(noisval)
    [nonmotdata,noisonly] = addnoise_outall(cdata,0,noisval);%data, mean, noise factor
    % since we are creating multiple motion scenarios, an ideal performance
    % comparison requires that the noise is the same in each scenario.
    % Thus the cdata replication should happen after noise generation
    cdata = repmat(nonmotdata,[1,1,1,nfreq,namp]);%essentially initialize the matrix of multiple motion scenarios
    for nf = 1:nfreq
        for na = 1:namp
            motdata(:,:,nf,na) = motdata(:,:,nf,na) + noisonly(:,:,1);%if you only generate one shift
            cdata(:,:,1,nf,na) = motdata(:,:,nf,na); %replace the origianl dataset with the motion one
        end
    end
end

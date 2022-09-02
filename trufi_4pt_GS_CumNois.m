%% trufi_4pt_GS_CumNois.m 20200412 plot cumulative noise and average deviation
%% vs normalized radiality 
% originally devised as the GS Noise radiality simulation for ISMRM2015

%% define some parameters, housekeeping
clear all
tic;
format long g
cyc = [0 1/4 2/4 3/4];
npc = length(cyc);
flip = 70; flipr = flip*pi/180; %in radians
TR = 4.2; TE = TR/2;

%% generate data AND random noise
noisval = input('Enter noise level (press return for default of 50) = ');
if isempty(noisval), noisval = 50;end
highres = input('Do you want highres (1) data or quick (0) = ');
if highres, th_pts = 1600; else, th_pts = 160; end% allow lower density for troubleshooting 
[amat,bmat,mmat,theta,compdata,E1mat,~,noisonly,puredata] = sim_bssfp_continue(cyc,noisval,th_pts,flipr,1,highres);

%% Calculate M-spoke Reference Frame
% calculate angles that bisecting lines make relative to the horizontal
[bisang13,bisang24] = calc_spoke_angle(puredata);
% rotate noise data clockwise by bisang13 and bisang24 angles  
[noisonly_rot] = noisrotate(noisonly,bisang13,bisang24);
% noisonly_rot has radial noise along the horizontal/real axis, while 
% tangential noise is along the vertical/imaginary axis

%% Calculate Radiality of Noise using reference frame
% Radiality is a radial - tangential normalized  noise metric with
% -1 (100% tangential) to +1 (100% % radial) range
[~,rtmet,radtotnorm] = comp_radiality_nois(noisonly_rot);
% this computes the average and the sum-norm (radtotnorm) radiality of the
% noise for the four phase cycles of a given pixel.
[nr,nc,~,~] = size(compdata);


%% compute GS, CS, LGS
pr = [1 3;2 4];
CS = abs(mean(compdata,3));%compsum 
[GS,GS_sing] = geosoln(compdata);
pairings = [1 3;2 4];    reg_siz = 5;  
[LGS,~] = linearizeGS(compdata,GS,reg_siz,pairings);

%% noise = solutions - gold standards
CSgold = mmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);
CSerr = abs(CS) - CSgold; CSmagnois = abs(CSerr);
GSerr = abs(GS) - mmat; GSmagnois = abs(GSerr);
LGSerr = abs(LGS) - mmat; LGSmagnois = abs(LGSerr);

%% bin the noise radiality data for TRE calculations
numbin = 50;
%define edges and centres of bins
binEdges = linspace(-1,1,numbin+1);
binCtrs = (binEdges(1:end-1) + binEdges(2:end))./2;
%compute counts and indices of radiality vals in each bin
[bincnt,binIdx] = histc(radtotnorm(:), binEdges);%for the normed

%% initialize radiality cells and arrays
% solutions 
GSradcell{numbin} = []; CSradcell{numbin} = []; LGSradcell{numbin} = [];
% gold standards
Mgoldradcell{numbin} = [];CSgoldradcell{numbin} = [];
%TRE
GSerrPerBin = zeros(1,numbin); CSerrPerBin = zeros(1,numbin);LGSerrPerBin = zeros(1,numbin); 
%STDDEV cells
GS_ERRcell{numbin} = [];CS_ERRcell{numbin} = [];LGS_ERRcell{numbin} = []; 
%STDDEV arrays
GS_STD = zeros(1,numbin); CS_STD = zeros(1,numbin);LGS_STD = zeros(1,numbin);
%Bias
GSbias = zeros(1,numbin); CSbias = zeros(1,numbin);LGSbias = zeros(1,numbin); 
%Confdence intervals
GS_CI95 = zeros(1,numbin); LGS_CI95 = zeros(1,numbin); CS_CI95 = zeros(1,numbin); 
%% compute mean sqare error and variation for each noise radiality bin
for bin = 1:numbin
    % find bin indices for the radiality vals
    indvec = find(binIdx == bin);
    % check that the counts are the same as calculated above
    if length(indvec) ~= bincnt(bin)
        error('there is a problem, sir!!')
    end
    %% compute MSE for each normed radiality bin
    if ~isempty(indvec)
        % allocate data for each bin into cells
        for radcnt = 1:bincnt(bin) 
            radind = indvec(radcnt);
            GSradcell{bin}(radcnt) = abs(GS(radind));
            LGSradcell{bin}(radcnt) = abs(LGS(radind));
            CSradcell{bin}(radcnt) = abs(CS(radind));
            Mgoldradcell{bin}(radcnt) = abs(mmat(radind));
            CSgoldradcell{bin}(radcnt) = abs(CSgold(radind));
    %         also set aside cell nois vals for STD calculations
            GS_ERRcell{bin}(radcnt) = GSradcell{bin}(radcnt)-Mgoldradcell{bin}(radcnt);
            LGS_ERRcell{bin}(radcnt) = LGSradcell{bin}(radcnt)-Mgoldradcell{bin}(radcnt);
            CS_ERRcell{bin}(radcnt) = CSradcell{bin}(radcnt)-CSgoldradcell{bin}(radcnt);
        end
        [~,nugo]=size(GS_ERRcell{bin});
        [~,nulo]=size(LGS_ERRcell{bin});
        [~,nuco]=size(CS_ERRcell{bin});
        if nugo ~= nulo || nugo ~= nuco
             error('There a different numbers of solutions? ');
        end
        if nuco > 2 
        
            %calculate the RMSE for each bin
            GSerrPerBin(bin)=  sqrt(sum((GS_ERRcell{bin}).^2)/nuco);
            LGSerrPerBin(bin)=  sqrt(sum((LGS_ERRcell{bin}).^2)/nuco);
            CSerrPerBin(bin)= sqrt(sum((CS_ERRcell{bin}).^2)/nuco);
            % calculate error bias
            GSbias(bin)=  mean(GS_ERRcell{bin});
            LGSbias(bin)=  mean(LGS_ERRcell{bin});
            CSbias(bin)= mean(CS_ERRcell{bin});
            %calculate the STD of each solutions error
            GS_STD(bin) = std(GS_ERRcell{bin});
            LGS_STD(bin) = std(LGS_ERRcell{bin});
            CS_STD(bin) = std(CS_ERRcell{bin});
            % calculate confidence intervals
            GS_CI95(bin) = GS_STD(bin)/sqrt(nuco);
            LGS_CI95(bin) = LGS_STD(bin)/sqrt(nuco);
            CS_CI95(bin) = CS_STD(bin)/sqrt(nuco);
        end
    end        
end

%% finishing plots
prnt = input(' Enter 1 if you want to print the figures = ');
maxmg = 2700;% maxmg = max(abs(GS(:)));  %set a max for plots to avoid singularities

%% plot parameters
close all
lw = 1.5; ptsize = 1; fsize = 17;
% since edge bins are sparsely populated, consider dropping from plots
first = ceil(numbin*0.03);last = floor(numbin*0.99);
subrange = input('Plot the full radiality range(0), or just the meaningful range (1) = ');
if subrange,rang = first:last; else, rang = 1:numbin;end

%% plot bias w 95% CI vs normalized noise radiality
BIASrad = figure(64);
if ~subrange,rang = (~GSbias==0);end
errorbar(binCtrs(rang),GSbias(rang),GS_CI95(rang),'LineWidth',lw,'Color',[0 0 1]);
hold on
if ~subrange,rang = (~LGSbias==0);end
errorbar(binCtrs(rang),LGSbias(rang),LGS_CI95(rang),'LineWidth',lw,'Color',[0 0.6 0]);
if ~subrange,rang = (~CSbias==0);end
errorbar(binCtrs(rang),CSbias(rang),CS_CI95(rang),'LineWidth',lw,'Color',[0.7 0 0]);
ylim([-0.5 2.5])
hleg = legend('GS','LGS','CS');
set(hleg,'FontSize',fsize-2,'FontWeight','bold');
% xlabel('Normalized Pixel Noise Radiality','FontWeight','Bold','FontSize',fsize);
% ylabel({'Bias'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
% annotation(BIASrad,'textbox',[0.15 0.13 0.047 0.109],'String',{'More',...
%     'Tangential','Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
% annotation(BIASrad,'textbox',[0.76 0.13 0.038 0.1057],'String',{'More',...
%     'Radial','Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
hold off

%% plot STD vs normalized noise radiality 
STDrad = figure(65);
if ~subrange,rang = (~GS_STD==0);end
plot(binCtrs(rang),GS_STD(rang),'LineWidth',lw,'Color',[0 0 1]);
hold on
if ~subrange,rang = (~LGS_STD==0);end
plot(binCtrs(rang),LGS_STD(rang),'LineWidth',lw,'Color',[0 0.6 0]); 
if ~subrange,rang = (~CS_STD==0);end
plot(binCtrs(rang),CS_STD(rang),'LineWidth',lw,'Color',[0.7 0 0]); 
hleg = legend('GS','LGS','CS');
set(hleg,'FontSize',fsize-2,'FontWeight','bold');
% xlabel('Normalized Pixel Noise Radiality','FontWeight','Bold','FontSize',fsize);
hold off

%% plot STD & Bias w 95% CI vs normalized noise radiality 
STDnBIAS = figure(66);
if ~subrange,rang = (~GS_STD==0);end
figGSTD = plot(binCtrs(rang),GS_STD(rang),'LineWidth',lw,'Color',[0 0 1]);
hold on
figGBIAS = errorbar(binCtrs(rang),GSbias(rang),GS_CI95(rang),'LineWidth',lw,'Color',[0 0 1]);
if ~subrange,rang = (~LGS_STD==0);end
figLSTD = plot(binCtrs(rang),LGS_STD(rang),'LineWidth',lw,'Color',[0 0.6 0]); 
figLBIAS = errorbar(binCtrs(rang),LGSbias(rang),LGS_CI95(rang),'LineWidth',lw,'Color',[0 0.6 0]);
if ~subrange,rang = (~CS_STD==0);end
figCSTD = plot(binCtrs(rang),CS_STD(rang),'LineWidth',lw,'Color',[0.7 0 0]); 
figCBIAS = errorbar(binCtrs(rang),CSbias(rang),CS_CI95(rang),'LineWidth',lw,'Color',[0.7 0 0]);
ylim([-3 47]);
% xlim([-0.8 0.8]);
legsub = [figGSTD figLSTD figCSTD];
hleg = legend(legsub,'GS','LGS','CS');
set(hleg,'FontSize',fsize-2,'FontWeight','bold');
% xlabel('Phase-cycled Normalized Pixel Radiality','FontWeight','Bold','FontSize',fsize);
hold off

%% pix-by-pix scatterplot of abs(GS/CS) error vs normalized noise radiality 
typicerror = 230;
% CS cumulative error plot
% CSnfig = figure(67);scatter(radtotnorm(:),CSerr(:),ptsize,[0.7 0 0],'.');axis equal
% xlim([-1 1]);ylim([-typicerror typicerror]); axis square

CSnfig = figure;scatter3(radtotnorm(:),CSerr(:),ptsize,[0.7 0 0],'.');axis equal
xlim([-1 1]);ylim([-typicerror typicerror]); axis square


% title('CS Error vs Normalized Pix-Noise Radiality'); 
% xlabel('Norm Pixel Noise Radiality \rho','FontWeight','Bold','FontSize',fsize);
% ylabel({'Complex','Sum','error'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
% annotation(CSnfig,'textbox',[0.21 0.11 0.047 0.109],'String',{'Tangential',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
% annotation(CSnfig,'textbox',[0.76 0.11 0.038 0.1057],'String',{'Radial',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
% GS cumulative error plot
GSnfig = figure(68);scatter(radtotnorm(:),GSerr(:),ptsize,'b','.');axis equal
xlim([-1 1]);ylim([-typicerror typicerror]); axis square
% title('GS Error vs Norm Pixel Noise Radiality'); 
% xlabel('Norm Pixel Noise Radiality \rho','FontWeight','Bold','FontSize',fsize);
% ylabel({'Geometric','Solution','error'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
% annotation(GSnfig,'textbox',[0.21 0.11 0.047 0.109],'String',{'Tangential',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
% annotation(GSnfig,'textbox',[0.76 0.11 0.038 0.1057],'String',{'Radial',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
%LGS cumulative error plot
LGSnfig = figure(69);scatter(radtotnorm(:),LGSerr(:),ptsize,[0 0.6 0],'.');axis equal
xlim([-1 1]);ylim([-typicerror typicerror]); axis square
% title('LGS Error vs Norm Pixel Noise Radiality'); 
% xlabel('Norm Pixel Noise Radiality \rho','FontWeight','Bold','FontSize',fsize);
% ylabel({'Linearized','Geometric','Solution','error'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
% annotation(LGSnfig,'textbox',[0.21 0.11 0.047 0.109],'String',{'Tangential',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);
% annotation(LGSnfig,'textbox',[0.76 0.11 0.038 0.1057],'String',{'Radial',...
%     'Noise'},'FitBoxToText','off','LineStyle','none','FontSize',fsize-2);

%% Source data plots
%phase cycles
for se = 1:4
    figure(1+se);colormap(gray);imagesc(abs(compdata(:,:,se)), [0 maxmg]);
    deg = sprintf([int2str(cyc(se)*360),'%c'], char(176));
    title(['\Delta\theta = ',deg],'FontSize',15);axis image;set(gca,'xtick',[],'ytick',[]);
    set(gca,'XTick',1:80*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
    [~,h1]=suplabel('Magnitude');
end

%% Gold standards
%M (for GS and LGS)
figure(14);colormap(gray);imagesc(mmat,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
title(['Base M, Noise=',int2str(noisval)],'FontSize',15);
set(gca,'XTick',1:80*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
% Center-of-mass (for CS)
figure(26);colormap(gray);imagesc(CSgold,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
title(['CS N = \inf, Noise=',int2str(noisval)],'FontSize',15);
set(gca,'XTick',1:80*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% CS
figure(16);colormap(gray);imagesc(CS,[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
title('CS','FontSize',15);
set(gca,'XTick',1:80*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% GS
figure(17);colormap(gray);imagesc(abs(GS),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);
title('GS','FontSize',15);
set(gca,'XTick',1:80*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% LGS
figure(18);colormap(gray);imagesc(abs(LGS),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);
title('LGS','FontSize',15);
set(gca,'XTick',1:80*th_pts:nc)
set(gca,'XTickLabel',amat(1,1:80*th_pts:nc));
set(gca,'YTick',nr/7:6*nr/7:nr)
set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% print images?

if prnt
    print -f2 -djpeg100 0PC.jpg
    print -f3 -djpeg100 90PC.jpg
    print -f4 -djpeg100 180PC.jpg
    print -f5 -djpeg100 270PC.jpg
    print -f16 -djpeg100 CS.jpg
    print -f17 -djpeg100 gs.jpg
    print -f18 -djpeg100 lgs.jpg
    print -f14 -djpeg100 Mgold.jpg
    print -f26 -djpeg100 CSgold.jpg
    print -f64 -djpeg100 BIASvsNoisRadiality.jpg
    print -f65 -djpeg100 STDvsNoisRadiality.jpg
    print -f66 -djpeg100 BiasSTDvsNoisRadiality.jpg
    print -f67 -djpeg100 CSerrorvsNormNoisRadiality.jpg
    print -f68 -djpeg100 GSerrorvsNoisRadiality.jpg
    print -f69 -djpeg100 LGSerrorvsNoisRadiality.jpg
 end

toc
%% SimMotionbSSFP_TREvsAmplFreq.m 20220124 
% this script simulates phase cycled bSSFP data with basic shapes, 1 with 
% translational motion.  This is used for figures (2, 3, 4) in the LGS 
% motion manuscript 
clear all
format long g
cyc = [0 1/4 2/4 3/4]';npc = length(cyc);

%% generate dataset parameters
rr = 1;cc = 4;
flip = 70; flipr = flip*pi/180; %in radians
TR = 4.2; TE = TR/2;

%% generate data AND random noise
noisval = input('Noise level? ENTER gives a standard deviation of 25 = ');
if isempty(noisval),noisval = 25;end
th_pts = 32;% th_pts on circle for 2pi rotation
agran = 9; % granularity of the "a" term
    
%% motion parameters "vers"
% default (0) varies motion amplitude AND frequency for the Fig 3 surface plot
% (1) tissue type and (2) noise for 4b/a (respectively) line plots
motdir = input(' Would you like horizontal (h) or vertical (v) motion? Default is vertical ','s');
if strcmpi(motdir,'h'), disp('motion is horizontal'); else disp('Motion is vertical'); end
PEdirn = input(' Which direction would you like to run phase encoding, (h), or (v)? ','s'); 
if strcmpi(PEdirn,'h'), disp('Phase encoding is horizontal'); else disp('Phase encoding is vertical'); end
vers = input('Would you like to test motion vs. motion amplitude and frequency(0), tissue type T1/T2(1), or noise(2)? ');
if vers == 1 %Fig 4b
    T1 = [3300 4000 280 1200 720 950 660 550 900 900]; T2 = [3300 2000 70 200 90 100 60 44 60 45];
    % water, CSF, Fat, Blood(oxy), White Matter, Gray matter, Kidney, liver, Muscle(heart), and muscle(skeletal), respectively
    T1T2rat = T1./T2; numparm = length(T1T2rat); noisval = 5;
    disp('Run a motion simulation with variable tissue type (T1/T2 ratio)');
    amplitude = input(' How many pixels of motion range would you like?  Click enter for 3 = ');if isempty(amplitude), amplitude = 3; end
    freq = input(' How many motion cycles would you like?  Default is 40 = ');if isempty(freq), freq = 40; end   
    for tt = 1:numparm
        [amat,bmat,mmat,theta,compdata(:,:,:,tt),E1mat,rotmat,noisonly,motdata(:,:,tt),puredata(:,:,:,tt),row1,rowlst,col1,collst,mask] = sim_bssfp_motOneTis(cyc,noisval,T1(tt),T2(tt),th_pts,flipr,1,TR,motdir,freq,amplitude,PEdirn);
    end
    [nr,nc,npc,numparm] = size(compdata);
elseif vers == 2 %Fig 4a
    disp('Run a motion simulation at variable noise levels');
    amplitude = input(' How many pixels of motion range would you like?  Default is 3 = ');if isempty(amplitude), amplitude = 3; end
    freq = input(' How many motion cycles would you like?  Default is 40 = ');if isempty(freq), freq = 40; end 
    noisval = [1 5:5:100]; numparm = length(noisval);
    for tt = 1:numparm
        [amat,bmat,mmat,theta,compdata(:,:,:,tt),E1mat,rotmat,noisonly,motdata(:,:,tt),puredata(:,:,:,tt),row1,rowlst,col1,collst,mask] = sim_bssfp_SpecNoisMot(cyc,noisval(tt),th_pts,flipr,1,agran,TR,motdir,freq,amplitude,PEdirn);
    end  
    [nr,nc,npc,numparm] = size(compdata);
else
    vers = 0; %Fig 3
    disp('Run a continuous variation of motion frequency and amplitude');             
    [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,nonmotdata,motdata,row1,rowlst,col1,collst,freqs,amps,mask] = sim_bssfp_shiftshape(cyc,noisval,th_pts,flipr,1,agran,TR,motdir,PEdirn); 
    [nr,nc,npc,parmf,parma] = size(compdata);
end 

%% Housekeeping: String and parameter assignment
if strcmpi(motdir,'h'), motdir = 'H'; else motdir = 'V'; end  
if exist('PEdirn')
    if strcmpi(PEdirn,'h'), PEdirn = 'H'; else PEdirn = 'V'; end 
end
T2 = -TR./log(amat);
E1 = (amat.*(1+cos(flipr))- bmat.*(1+amat.^2.*cos(flipr)))./(amat.*(1+cos(flipr))-bmat.*(amat.^2+cos(flipr)));
T1 = -TR./log(E1);

%% Complex Sum
CS = squeeze(abs(mean(compdata,3)));%compsum 
CSgold = mmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);
CSgold(isnan(CSgold)) = 0; 

%% Linearized Geometric Solution
[GS,~] = geosoln(compdata);%this should work for up to 6D data
reg_siz = 5; ptpairs = [1 3;2 4]; % region size and point pairs for linearization
[LGS,~] = linearizeGS(compdata,GS,reg_siz,ptpairs);

%% set aside a non-motion dataset for the T1/T2 and noise variability versions
if vers == 1 || vers == 2
    CS1 = squeeze(abs(mean(puredata,3)));
    [GS1,~] = geosoln(puredata);
    [LGS1,~] = linearizeGS(puredata,GS1,reg_siz,ptpairs); 
end

%% Total Relative Error (TRE)
% compute error for motion-corrupted solutions relative to "gold standards" 
% formed from solutions to non-motion data
if vers == 1 || vers == 2 %line plots
    for tt = 1:numparm
        GSt = GS(:,:,tt);CSt = CS(:,:,tt);LGSt = LGS(:,:,tt);% temp motion-corrupted solns for TRE
        GS1t = GS1(:,:,tt);CS1t = CS1(:,:,tt);LGS1t = LGS1(:,:,tt);% temp motion-uncorrupted solns for TRE
        TRE_GS(tt) = totrelerrorTRE(GSt(:),GS1t(:));
        TRE_CS(tt) = totrelerrorTRE(CSt(:),CS1t(:));
        TRE_LGS(tt) = totrelerrorTRE(LGSt(:),LGS1t(:)); 
    end
else %(vers = 0), surface plot
    GS1 = GS(:,:,1,1);% temp motion-uncorrupted 
    CS1 = CS(:,:,1,1);% temp motion-uncorrupted 
    LGS1 = LGS(:,:,1,1); %% temp motion-uncorrupted LGS
    for nf = 1:parmf
        for na = 1:parma 
            % temp motion-corrupted solns for TRE
            GSt = GS(:,:,nf,na);CSt = CS(:,:,nf,na);LGSt = LGS(:,:,nf,na);
            %TRE calculation for motion error
            TRE_GS(nf,na) = totrelerrorTRE(GSt(:),GS1(:));
            TRE_CS(nf,na) = totrelerrorTRE(CSt(:),CS1(:));
            TRE_LGS(nf,na) = totrelerrorTRE(LGSt(:),LGS1(:)); 
        end
    end
end

%% TRE surface and line plots
lw = 3; ptsize = 250; fsize = 13; lastpt = length(TRE_CS); %plot parameters
if vers == 0  %variation of motion amplitude AND frequency
    figure(82);surf(amps(2:end),freqs,TRE_CS(:,2:end),'Facecolor','r')
    alpha(0.75)
    hold on
    surf(amps(2:end),freqs,TRE_LGS(:,2:end),'Facecolor','g');alpha(0.75)
    zlabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);%;
    ylabel({'Motion Frequency';'[cycles / dataset]'},'FontWeight','Bold','FontSize',fsize-2);
    xlabel({'Motion Amplitude';'[Max. Pixel Shift]'},'FontWeight','Bold','FontSize',fsize-2);
    title(['1/4 Image ',motdir,'-dirn Motion TRE vs. Cycles & Pixel Shift, Flip',int2str(flip),'\circ'],'FontSize',fsize-2);
    hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['Surf_TREvs' motdir 'dirnFlip' int2str(flip) 'PixelShiftnMotionCycles'];
elseif vers == 1
    figure(82); %plot TRE vs tissue type T1/T2
    plot(T1T2rat,TRE_CS,'LineWidth',lw,'Color','r');
    hold on
    plot(T1T2rat,TRE_LGS,'LineWidth',lw,'Color','g');
    ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
    xlabel('T1/T2','FontWeight','Bold','FontSize',fsize);
    axis square
    hold off
    hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['TREvsT1T2' PEdirn 'dirnPE_' motdir 'dirn'  int2str(amplitude) 'Pixel' int2str(freq) 'Cyc_Motion'];   
elseif vers == 2
    figure(82); %plot TRE vs noise
    [~,~,nzmmat] = find(mmat);%non-zero elements of pure GS 
    nrm = mean(nzmmat(:));% an average pure GS over all signal areas
    plot(100*noisval./nrm,TRE_CS,'LineWidth',lw,'Color','r');
    hold on
    plot(100*noisval./nrm,TRE_LGS,'LineWidth',lw,'Color','g');
    y = ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
    set(y, 'position', get(y,'position')-[0,0,0]); 
    x = xlabel('Noise % (Standard Deviation/Mean Signal)','FontWeight','Bold','FontSize',fsize);
    set(x, 'position', get(x,'position')-[0,0,0]); 
    xlim([0 max(noisval*100/nrm)]);
    axis square
    hold off
    hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['TREvsNoise' PEdirn 'dirnPE_' motdir 'dirn'  int2str(amplitude) 'Pixel' int2str(freq) 'Cyc_Motion']; 
end

ifplot = input('Enter 1 to plot base images = ');
if ifplot
%% finishing plots
maxmg = 0.65*max(abs(compdata(:)));%set a maxval for optimal window
%base image data
if vers == 0
    for na = 1:parma
        for nf = 1:parmf 
            for se = 1:4
                figure(100+10*nf+na);subplot(rr,cc,se);colormap(gray);imagesc(abs(compdata(:,:,se,nf,na)), [0 maxmg]);
                title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);
            end
            figure;colormap(gray);imagesc(CS(:,:,nf,na));axis image; set(gca,'xtick',[],'ytick',[]);
            title([int2str(amps(na)),' PixShift, ',int2str(freqs(nf)),' Cyc CS']);
            figure;colormap(gray);imagesc(abs(LGS(:,:,nf,na)));axis image;set(gca,'xtick',[],'ytick',[]);
            title([int2str(amps(na)),' PixShift, ',int2str(freqs(nf)),' Cyc LGS']);
        end
    end

%% just one
else
    %base image data
    for np = 1:numparm 
        for se = 1:4
            figure(100+np);subplot(rr,cc,se);colormap(gray);imagesc(abs(compdata(:,:,se,np)));
            title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);  
        end
    end
    
    figure;colormap(gray);imagesc(CS(:,:,6));axis image;
    set(gca,'xtick',[],'ytick',[]);
    figure;colormap(gray);imagesc(abs(LGS(:,:,6)));axis image;
    set(gca,'xtick',[],'ytick',[]);


%% CS
    figure;colormap(gray);imagesc(CS(:,:,np));axis image;
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% GS
    figure;colormap(gray);imagesc(abs(GS(:,:,np)));
    axis image;set(gca,'xtick',[],'ytick',[]);
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

%% LGS = Linearized GS
    figure;colormap(gray);imagesc(abs(LGS(:,:,np)));
    axis image;set(gca,'xtick',[],'ytick',[]);
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);

end
end
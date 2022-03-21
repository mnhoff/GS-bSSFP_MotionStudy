%% trufi_4pt_GS_TREvsMotionSim.m 170104
% this version of trufi_4pt_GSvsMotion cuts out extra stuff, and focuses on 
% streamlining for plotting multiple recon graphs

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
    numtis = 0;
%     numtis = input('Continuous (0), or Multi-Tissue (1) simulations = ');
    pr = input('Do you want to print images? (1 or 0)');
    % pr = 1; 
    noispt = 1;
    if numtis %tri-tissue
        rowrep = 30; nc = 192;
        T1 = [420 240 2400];%liver, fat, CSF
        T2 = [50 70 1400];%liver, fat, CSF
        tiss = length(T1);
%         [amat,bmat,mmat,theta,compdata] = sim_bssfp_3tis(cyc,noisval,rowrep,T1,T2,noispt,flipr,TR);
        [amat,bmat,mmat,theta,compdata,rotmat] = sim_bssfp_MultiTisStacked(cyc,noisval,T1,T2,flipr,TR,rowrep,nc);
    else %continuous data
        th_pts = 32;% %formerly 16 when I had 4 "a" vals, but I want a 2pi rotation
        agran = 9; 
%         agran = input('Enter the ''a'' granularity (9->900, hit enter for default of 9)  = ');
        if isempty(agran), agran = 9; end % agran = 2 gives 4 bands with each th_pts*npi large, nice and sparse. 9 actually gives 2 avals, better for npi = 2
        motdir = input(' Would you like horizontal (h) or vertical (v) motion? Default is vertical ','s');
        vers = input('Would you like to test motion amp(0),vs tissue type(1), vs noise(2), vs flip(3), freq(4), amp and freq (5)? ');
        %ghost intesity?  Datasets corrupted?;
%         vers = 5;% temporary!!!
        motsusc = 'y';
%         motsusc = input('Calculate only the motion error(y) or (n)? ','s');
        splitsignois = input('Would you like to make separate calculations for motion artifact in signal and noise regions? (y or n)? ','s');
        if strcmpi(motsusc,'y'), MS = 'Mot'; else; MS = 'MotSusc'; end  
        if vers == 1
            noisval = 5;
            disp('Let us run a motion simulations for variable tissue type(T1/T2 ratio)');
            amplitude = input(' How many pixels of motion range would you like?  Default is 3 = ');if isempty(amplitude), amplitude = 3; end
            freq = input(' How many motion cycles would you like?  Default is 40 = ');if isempty(freq), freq = 40; end 
            PEdirn = input(' Which direction would you like to run phase encoding, (h), or (v)? ','s');
            T1 = [3300 4000 280 1200 720 950 660 550 900 900]; T2 = [3300 2000 70 200 90 100 60 44 60 45];
            T1T2rat = T1./T2; numparm = length(T1T2rat);
            % these combos are water, CSF, Fat, Blood(oxy), White Matter, Gray matter, Kidney, liver, Muscle(heart), and muscle(skeletal), respectively
            for tt = 1:numparm
                [amat,bmat,mmat,theta,compdata(:,:,:,tt),E1mat,rotmat,noisonly,motdata(:,:,tt),puredata(:,:,:,tt),row1,rowlst,col1,collst,mask] = sim_bssfp_1tisMot(cyc,noisval,T1(tt),T2(tt),th_pts,flipr,1,TR,motdir,freq,amplitude,PEdirn);
            end
        elseif vers == 2
            disp('Let us run a motion simulation at variable noise levels');
            amplitude = input(' How many pixels of motion range would you like?  Default is 3 = ');if isempty(amplitude), amplitude = 3; end
            freq = input(' How many motion cycles would you like?  Default is 40 = ');if isempty(freq), freq = 40; end 
            PEdirn = input(' Which direction would you like to run phase encoding, (h), or (v)? ','s');
            noisval = [1 5:5:100]; numparm = length(noisval);
            for tt = 1:numparm
                [amat,bmat,mmat,theta,compdata(:,:,:,tt),E1mat,rotmat,noisonly,motdata(:,:,tt),puredata(:,:,:,tt),row1,rowlst,col1,collst,mask] = sim_bssfp_continuous_noisMot(cyc,noisval(tt),th_pts,flipr,1,agran,TR,motdir,freq,amplitude,PEdirn);
            end            
        elseif vers == 3
            disp('Let us run variable-flip angle motion simulation');% can't do this till you fix the nightmare
            [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,~] = sim_bssfp_continuous_flip(cyc,noisval,th_pts,1,agran);
        elseif vers == 4
            disp('Let us run variable-motion cycle (ghost separation) simulation');            
            amplitude = input(' How many pixels of motion range would you like?  Default is 5 = ');
            if isempty(amplitude), amplitude = 5; end
            [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,motdata,~,row1,rowlst,col1,collst,freqs,mask] = sim_bssfp_continuous_freq(cyc,noisval,th_pts,flipr,1,agran,TR,motdir,amplitude,pr); 
        elseif vers == 5
            disp('Let''s run a continuous variation of motion frequency and amplitude');             
            if strcmpi(motdir,'h'), motdir = 'F'; else; motdir = 'P'; end 
            [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,nonmotdata,motdata,row1,rowlst,col1,collst,freqs,amplitudes,mask] = sim_bssfp_wmot(cyc,noisval,th_pts,flipr,1,agran,TR,motdir,pr); 
        else
            vers = 0;
            disp('Let us run a variable motion amplitude simulation, potentially with T1/T2 variability in analysis ');
            freq = input(' How many motion cycles would you like?  Default is 30 = '); if isempty(freq), freq = 30; end
            PEdirn = input(' Which direction would you like to run phase encoding, (h), or (v)? ','s');%I don't think PE can be run horizontally unless you change simsanmotion from row sampling to column sampling
%             [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonlys,motdata,~] = sim_bssfp_continuous_shft2(cyc,noisval,th_pts,flipr,1,agran,TR);
%             [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonlys,motdata,~,row1,rowlst,col1,collst] = sim_bssfp_continuous_shft3(cyc,noisval,th_pts,flipr,1,agran,TR);
            [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,motdata,~,row1,rowlst,col1,collst,amplitudes,mask] = sim_bssfp_continuous_shft4(cyc,noisval,th_pts,flipr,1,agran,TR,motdir,freq,PEdirn);
%         disp('Let us corrupt a variable amount of datasets with motion');
%             [amat,bmat,mmat,theta,compdata,E1mat,rotmat,noisonly,~] = sim_bssfp_continuous_corrsets(cyc,noisval,th_pts,flipr,1,agran);
        end           
        if strcmpi(motdir,'h'), motdir = 'H'; else motdir = 'V'; end  
        if exist('PEdirn')
            if strcmpi(PEdirn,'h'), PEdirn = 'H'; else PEdirn = 'V'; end 
        end
        %% what is the corresponding T1 and T2?     
        T2 = -TR./log(amat);
        E1 = (amat.*(1+cos(flipr))- bmat.*(1+amat.^2.*cos(flipr)))./(amat.*(1+cos(flipr))-bmat.*(amat.^2+cos(flipr)));
        T1 = -TR./log(E1);
        T1T2 = T1./T2;
    end
    if vers == 0
        [~,~,npcc,parma] = size(compdata);% this is sloppy, better to 
% freeze each parameter in place in the simulations and then squeeze later
    elseif vers == 5
        [~,~,npcc,parmf,parma] = size(compdata);%nfreq and namp make up a 2D set of 3D phase-cycled, motion corrupted datasets
    else
        [~,~,npcc,~,~] = size(compdata);
    end
if npc ~= 4, error('This script only runs quad-phase cycles'); end

    %limit FOV?
%     rows = 1:nr;   cols = 1:nc; %full
%     compdata = compdata(rows,cols,:,:);

if npc ~= npcc, error('There appears to be a mismatch with the number of phase cycles'); end
  [nr,nc,~,~,~] = size(compdata);  

%% compute the GS and other quantities
CScomp = squeeze(mean(compdata,3));
CS = squeeze(abs(mean(compdata,3)));%compsum 

tic
% % [GS,GS_sing] = xs(compdata);%GS
[GS,~] = xs_multicoil_multiphase(compdata);%this should work for up to 6D data

%%??? continue edits here
%% LGS 
% reg_siz = input(' Please enter the desired length of the square ROI = ');
reg_siz = 5;
isop = pi/1;
% [LGS,~] = secpass(compdata,GS,reg_siz,[1 3;2 4]);  %LGS 
[LGS,~] = secpassMultiCoilPhase(compdata,GS,reg_siz,[1 3;2 4],isop);

% if you don't have a 0-motion dataset in the 1st motion-amplitude position
if vers == 1 || vers == 2
    CS1 = squeeze(abs(mean(puredata,3)));
    [GS1,~] = xs(puredata);
    [LGS1,~] = secpass(puredata,GS1,reg_siz,[1 3;2 4]); 
end

toc

%% plotter?
% paus = 0.2;
% % % % % [numpt] = plotter_adv(compdata,abs(compdata(:,:,1)),GS,paus);
% % % % % minr = -500;maxr = 500;
% % % % %plot GS
% % % % plotslice = 18; %this happens to be the slice with the minimum motion
plotslice = 1;%for simulations
% [numpt] = plotter_adv2(compdata(row1:rowlst,col1:collst,:,plotslice),abs(compdata(row1:rowlst,col1:collst,3,plotslice)),LGS(row1:rowlst,col1:collst,plotslice),paus);
%or plot 2nd base image
% [numpt] = plotter_adv2(compdata(row1:rowlst,col1:collst,:),abs(compdata(row1:rowlst,col1:collst,1)),abs(compdata(row1:rowlst,col1:collst,2)),paus,minr,maxr);

%% Gold Standards
CSgold = mmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);
CSgoldrot = rotmat.*(1+(bmat-amat).*(1./sqrt(1-bmat.^2)-1)./bmat);
CSgold(isnan(CSgold)) = 0; CSgoldrot(isnan(CSgoldrot)) = 0;

%% Total Relative Error
if numtis
    for tis = 1:length(T1)
        %section tissue recons and their gold standards
        GS_tis = GS((tis-1)*rowrep+1:rowrep*tis,1:noispt*nc);
        CS_tis = CS((tis-1)*rowrep+1:rowrep*tis,1:noispt*nc);
        mmat_tis = mmat((tis-1)*rowrep+1:rowrep*tis,1:noispt*nc);   
        CSgold_tis = CSgold((tis-1)*rowrep+1:rowrep*tis,1:noispt*nc);
        LGS_tis = LGS((tis-1)*rowrep+1:rowrep*tis,1:noispt*nc);
        if strcmpi(motsusc,'n')%calculate TRE for susceptibility and motion
            %compute corresponding TRE and store in array for each tissue
            TRE_GS(tis) = totrelerrorTRE(GS_tis,mmat_tis);
            TRE_CS(tis) = totrelerrorTRE(CS_tis,CSgold_tis);
            TRE_LGS(tis) = totrelerrorTRE(LGS_tis,mmat_tis);
        else %calculate TRE for motion only
            % do later: mirror continuous dataset below
        end
        clear GS_tis CS_tis mmat_tis CSgold_tis 
    end
else   %for the whole dataset
    %  form "t" truncated datasets which disinclude edge rows and columns
    % frow = 10; lrow = nr-9; fcol = 10; lcol = nc-9;%I normally use this edge
    % filter, but I'm thinking that all motion should be considered
    %or include them all
    frow = 1; lrow = nr; fcol = 1; lcol = nc;
    mmat = mmat(frow:lrow,fcol:lcol);
    CSgold = CSgold(frow:lrow,fcol:lcol);
    GS = GS(frow:lrow,fcol:lcol,:,:);
    LGS = LGS(frow:lrow,fcol:lcol,:,:);
    CS = CS(frow:lrow,fcol:lcol,:,:);
%     TRE_GS = zeros(parmf,parma); TRE_CS = zeros(parmf,parma);
%     TRE_LGS = zeros(parmf,parma);
    %the following are supposed to be the non-motion solutions, the first
    %dataset.  But then you have to reconstruct it!!
%     if vers ~= 5        % see above, I'leave this as standard for all
%     versions
% here are the solutions w/o motion to be used as gold standards for motion
% calculations
    if vers == 1 || vers == 2% simple calculation of TRE for each of the T1/T2 datasets
        for tt = 1:numparm
            GSt = GS(:,:,tt);CSt = CS(:,:,tt);LGSt = LGS(:,:,tt);
            GS1t = GS1(:,:,tt);CS1t = CS1(:,:,tt);LGS1t = LGS1(:,:,tt);
            TRE_GS(tt) = totrelerrorTRE(GSt(:),GS1t(:));
            TRE_CS(tt) = totrelerrorTRE(CSt(:),CS1t(:));
            TRE_LGS(tt) = totrelerrorTRE(LGSt(:),LGS1t(:)); 
        end
    else
        GS1 = GS(:,:,1,1);
        CS1 = CS(:,:,1,1);
        LGS1 = LGS(:,:,1,1);
        % separate signal and noise regions
        %1st for the non-motion data
        inGS1 = GS1.*mask; outGS1 = GS1.*imcomplement(mask);
        inCS1 = CS1.*mask; outCS1 = CS1.*imcomplement(mask);
        inLGS1 = LGS1.*mask; outLGS1 = LGS1.*imcomplement(mask);
%         inmmat = mmat.*mask; outmmat = mmat.*imcomplement(mask);
%         inCSgold = CSgold.*mask; outCSgold = CSgold.*imcomplement(mask);        
        
        % then for the recon datasets of noisy settings
        inGS = GS.*mask; outGS = GS.*imcomplement(mask);
        inCS = CS.*mask; outCS = CS.*imcomplement(mask);
        inLGS = LGS.*mask; outLGS = LGS.*imcomplement(mask);        
%     end
        if strcmpi(motsusc,'n')%calculate TRE for susceptibility and motion
            for nf = 1:parmf
                for na = 1:parma        
                    GSt = GS(:,:,nf,na);CSt = CS(:,:,nf,na);LGS(:,:,nf,na);
                    TRE_GS(nf,na) = totrelerrorTRE(GSt(:),mmat(:));
                    TRE_CS(nf,na) = totrelerrorTRE(CSt(:),CSgold(:));
                    TRE_LGS(nf,na) = totrelerrorTRE(LGSt(:),mmat(:));
                end
            end
        else %calculate TRE for motion only
            if splitsignois == 'y' % make calculations for bothe signal and noise regions separately
                for nf = 1:parmf
                    for na = 1:parma 
                        %regions
                        inGSt = inGS(:,:,nf,na);outGSt = outGS(:,:,nf,na);
                        inCSt = inCS(:,:,nf,na);outCSt = outCS(:,:,nf,na);
                        inLGSt = inLGS(:,:,nf,na);outLGSt = outLGS(:,:,nf,na);
                        %inner, signal TRE calculations motion only
                        TRE_inGS(nf,na) = totrelerrorTRE(inGSt(:),inGS1(:));
                        TRE_inCS(nf,na) = totrelerrorTRE(inCSt(:),inCS1(:));
                        TRE_inLGS(nf,na) = totrelerrorTRE(inLGSt(:),inLGS1(:));
    %                     motion and Susceptibility
    %                     TRE_inGSms(nf,na) = totrelerrorTRE(inGSt(:),inmmat(:));
    %                     TRE_inCSms(nf,na) = totrelerrorTRE(inCSt(:),inCSgold(:));
    %                     TRE_inLGSms(nf,na) = totrelerrorTRE(inLGSt(:),inmmat(:));
                        %outer, noise regions motion only
                        TRE_outGS(nf,na) = totrelerrorTRE(outGSt(:),outGS1(:));
                        TRE_outCS(nf,na) = totrelerrorTRE(outCSt(:),outCS1(:));
                        TRE_outLGS(nf,na) = totrelerrorTRE(outLGSt(:),outLGS1(:));
    %                     motion and Susceptibility
    %                     TRE_outGSms(nf,na) = totrelerrorTRE(outGSt(:),outmmat(:));
    %                     TRE_outCSms(nf,na) = totrelerrorTRE(outCSt(:),outCSgold(:));
    %                     TRE_outLGSms(nf,na) = totrelerrorTRE(outLGSt(:),outmmat(:));
                    end
                end
            end
            %don't split the signal and noise
            if vers == 0% I ultimately abandon this method in favor of generating new datasets for each T1T2 ratio
                T1T2tre = input('Would you like to calculate TRE as a function T1/T2 ratio (y) or (n)? ','s');
                if T1T2tre == 'y'% replace the number of frequencies with the T1/T2 ratio of different rows
                    for na = 1:parma
                        for r = 1:nr
                            GSt = GS(r,:,na);CSt = CS(r,:,na);LGSt = LGS(r,:,na);
                            GS1t = GS1(r,:); CS1t = CS1(r,:); LGS1t = LGS1(r,:);
                            TRE_GS(r,na) = totrelerrorTRE(GSt(:),GS1t(:));
                            TRE_CS(r,na) = totrelerrorTRE(CSt(:),CS1t(:));
                            TRE_LGS(r,na) = totrelerrorTRE(LGSt(:),LGS1t(:)); 
                        end
                    end
                else 
                    for na = 1:parma 
                        GSt = GS(:,:,na);CSt = CS(:,:,na);LGSt = LGS(:,:,na);
                        TRE_GS(na) = totrelerrorTRE(GSt(:),GS1(:));
                        TRE_CS(na) = totrelerrorTRE(CSt(:),CS1(:));
                        TRE_LGS(na) = totrelerrorTRE(LGSt(:),LGS1(:)); 
                    end
                end
            else
                for nf = 1:parmf
                    for na = 1:parma 
                        GSt = GS(:,:,nf,na);CSt = CS(:,:,nf,na);LGSt = LGS(:,:,nf,na);
                        TRE_GS(nf,na) = totrelerrorTRE(GSt(:),GS1(:));
                        TRE_CS(nf,na) = totrelerrorTRE(CSt(:),CS1(:));
                        TRE_LGS(nf,na) = totrelerrorTRE(LGSt(:),LGS1(:)); 
                    end
                end
            end
        end
    end
end

% %% for a test, let's see what just the susceptibity looks like
%     
% TRE_GS_susc = totrelerrorTRE(GS1(:),mmatt(:));
% TRE_CS_susc = totrelerrorTRE(CS1(:),CSgoldt(:));
% TRE_LGS_susc = totrelerrorTRE(LGS1(:),mmatt(:));

%% TRE plots
% note that my TRE calculations have taken on the format of: 
% amplitudes(2:end),freqs,TRE_"""(:,2:end)
% The TRE at 0 amplitude is 0...ultimately with transparent surfaces, I
% will include it.  There is no need for a freq = 0 so I start at a number,
% so this may be included
lw = 3;
ptsize = 250;
fsize = 13;
lastpt = length(TRE_CS);
if vers == 5  %meshplot for the ampfreq dual variation
    if splitsignois == 'y'
        figure(80);%signal-region motion correction mesh-plot
        surf(amplitudes(2:end),freqs,TRE_inCS(:,2:end),'Facecolor','r')
        alpha(0.75)
        hold on
    %     surf(amplitudes,freqs,TRE_inGS,'Facecolor','b')
        surf(amplitudes(2:end),freqs,TRE_inLGS(:,2:end),'Facecolor','g');alpha(0.75);
        z = zlabel({'Motion';'TRE';'Signal Region'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);%;
        set(z, 'position', get(z,'position')-[0,0,0]); 
        y = ylabel({'Motion Frequency';'[cycles / dataset]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',az
        set(y, 'position', get(y,'position')-[0,0,0]); 
        x = xlabel({'Motion Amplitude';'[Max. Pixel Shift]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',el
        set(x, 'position', get(x,'position')-[0,0,0]); 
        title(['1/4 Image ',motdir,'-dirn Signal-Region-Motion TRE vs. Cycles & Pixel Shift, Flip',int2str(flip),'\circ'],'FontSize',fsize-2);
        hleg = legend('Complex Sum','Linearized Geometric Solution');
        set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
        sigmot_str = ['SurfSig' MS 'TREvs' motdir 'dirnFlip' int2str(flip) 'PixelShiftvsMotionCycles'];
        print(gcf,'-djpeg100',sigmot_str);
        
        figure(81);%noise-region motion correction mesh-plot
        surf(amplitudes(2:end),freqs,TRE_outCS(:,2:end),'Facecolor','r')
        alpha(0.75)
        hold on
    %     surf(amplitudes,freqs,TRE_outGS,'Facecolor','b')
        surf(amplitudes(2:end),freqs,TRE_outLGS(:,2:end),'Facecolor','g');alpha(0.75)
        zlabel({'Motion';'TRE';'Noise Region'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);%;
        ylabel({'Motion Frequency';'[cycles / dataset]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',az
        xlabel({'Motion Amplitude';'[Max. Pixel Shift]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',el
        title(['1/4 Image ',motdir,'-dirn Noise-Region-Motion TRE vs. Cycles & Pixel Shift, Flip',int2str(flip),'\circ'],'FontSize',fsize-2);
        hleg = legend('Complex Sum','Linearized Geometric Solution');
        set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');        
        noismot_str = ['SurfNois' MS 'TREvs' motdir 'dirnFlip' int2str(flip) 'PixelShiftvsMotionCycles'];
        print(gcf,'-djpeg100',noismot_str);
    end
    figure(82);%entire image motion correction mesh-plot
    surf(amplitudes(2:end),freqs,TRE_CS(:,2:end),'Facecolor','r')
    alpha(0.75)
    hold on
%     surf(amplitudes,freqs,TRE_GS,'Facecolor','b')
    surf(amplitudes(2:end),freqs,TRE_LGS(:,2:end),'Facecolor','g');alpha(0.75)

    zlabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);%;
    ylabel({'Motion Frequency';'[cycles / dataset]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',az
    xlabel({'Motion Amplitude';'[Max. Pixel Shift]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',el
    title(['1/4 Image ',motdir,'-dirn Motion TRE vs. Cycles & Pixel Shift, Flip',int2str(flip),'\circ'],'FontSize',fsize-2);
    hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['Surf' MS 'TREvs' motdir 'dirnFlip' int2str(flip) 'PixelShiftnMotionCycles'];
elseif vers == 0 
    figure(82);
    if T1T2tre == 'y'%image motion correction mesh-plot varying T1T2 and amplitude
        surf(amplitudes(2:end),1:nr,TRE_CS(:,2:end),'Facecolor','r')% you may need to drop off the non-signal rows
        alpha(0.75)
        hold on
        %     surf(amplitudes,1:nr,TRE_GS,'Facecolor','b')
        surf(amplitudes(2:end),1:nr,TRE_LGS(:,2:end),'Facecolor','g');alpha(0.75)

        zlabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);%;
        ylabel({'T1/T2 Ratio';'[units]'},'FontWeight','Bold','FontSize',fsize-2);% need to scale this
        xlabel({'Motion Amplitude';'[Max. Pixel Shift]'},'FontWeight','Bold','FontSize',fsize-2);%,'Rotation',el
        title(['1/4 Image ',motdir,'-dirn Motion TRE vs. T1/T2 Ratio & Pixel Shift, Flip',int2str(flip),'\circ'],'FontSize',fsize-2);
        hleg = legend('Complex Sum','Linearized Geometric Solution');
        set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
        mot_str = ['Surf' MS 'TREvs' motdir 'dirn' PEdirn 'dirnPE_Flip' int2str(flip) 'T1T2nPixelShift'];
    else %just plot TRE vs amplitude
        plot(amplitudes(1:lastpt),TRE_CS(1:lastpt),'LineWidth',lw,'Color','r');
        % set(0,'DefaultAxesColorOrder',[0 0 0]);
        hold on
        plot(amplitudes(1:lastpt),TRE_GS(1:lastpt),'LineWidth',lw,'Color','b');
        plot(amplitudes(1:lastpt),TRE_LGS(1:lastpt),'LineWidth',lw,'Color','g');
    %     ....'Color','k');
    %    ...'Color','c');
        ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
        xlabel('Pixel Shift','FontWeight','Bold','FontSize',fsize);
        title(['TRE vs ',int2str(freq),'-Cycle Pixel Shift in ',motdir,'-Direction'],'FontSize',fsize); axis square
        hold off
        hleg = legend('Complex Sum','Geometric Solution','Linearized Geometric Solution');
        set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
        mot_str = ['TREvs' motdir 'dirn' int2str(freq) 'CyclePixelShift'];
    end
elseif vers == 4 
    figure(82); %plot TRE vs frequency
    plot(freqs(1:lastpt),TRE_CS(1:lastpt),'LineWidth',lw,'Color','r');
    % set(0,'DefaultAxesColorOrder',[0 0 0]);
    hold on
    plot(freqs(1:lastpt),TRE_GS(1:lastpt),'LineWidth',lw,'Color','b');
    plot(freqs(1:lastpt),TRE_LGS(1:lastpt),'LineWidth',lw,'Color','g');
    ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
    xlabel('Motion Cycles','FontWeight','Bold','FontSize',fsize);
    title(['TRE vs ',int2str(amplitude),'-Pixel-Shift Motion Cycles in ',motdir,'-Direction'],'FontSize',fsize); axis square
    hold off
    hleg = legend('Complex Sum','Geometric Solution','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['TREvs' motdir 'dirn' int2str(amplitude) 'PixelShiftMotionCycles'];
elseif vers == 1
    figure(82); %plot TRE vs frequency
    plot(T1T2rat,TRE_CS,'LineWidth',lw,'Color','r');
    % set(0,'DefaultAxesColorOrder',[0 0 0]);
    hold on
%     plot(T1T2rat,TRE_GS,'LineWidth',lw,'Color','b');
    plot(T1T2rat,TRE_LGS,'LineWidth',lw,'Color','g');
    ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
    xlabel('T1/T2','FontWeight','Bold','FontSize',fsize);
%     title(['TRE vs T1/T2, ',int2str(amplitude),'-Pixel, ' int2str(freq) '-Cycle ' motdir '-dirn Motion w ' PEdirn '-dirn PE'],'FontSize',fsize); 
    axis square
    hold off
%     hleg = legend('Complex Sum','Geometric Solution','Linearized Geometric Solution');
    hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['TREvsT1T2' PEdirn 'dirnPE_' motdir 'dirn'  int2str(amplitude) 'Pixel' int2str(freq) 'Cyc_Motion'];   
elseif vers == 2
    figure(82); %plot TRE vs noise
    [~,~,nzmmat] = find(mmat);%non-zero elements of pure GS 
    nrm = mean(nzmmat(:));% an average pure GS over all signal areas
    plot(100*noisval./nrm,TRE_CS,'LineWidth',lw,'Color','r');
    % set(0,'DefaultAxesColorOrder',[0 0 0]);
    hold on
%     plot(100*noisval./nrm,TRE_GS,'LineWidth',lw,'Color','b');
    plot(100*noisval./nrm,TRE_LGS,'LineWidth',lw,'Color','g');
    y = ylabel({'Motion';'TRE'},'FontWeight','Bold','FontSize',fsize,'Rotation',0);
    set(y, 'position', get(y,'position')-[0,0,0]); 
    x = xlabel('Noise % (Standard Deviation/Mean Signal)','FontWeight','Bold','FontSize',fsize);
    set(x, 'position', get(x,'position')-[0,0,0]); 
    xlim([0 max(noisval*100/nrm)]);
%     title(['TRE vs Noise, ',int2str(amplitude),'-Pixel, ' int2str(freq) '-Cycle ' motdir '-dirn Motion w ' PEdirn '-dirn PE'],'FontSize',fsize); 
    axis square
    hold off
%     hleg = legend('Complex Sum','Geometric Solution','Linearized Geometric Solution');
     hleg = legend('Complex Sum','Linearized Geometric Solution');
    set(hleg,'FontSize',fsize-4,'FontWeight','bold','Location','northwest');
    mot_str = ['TREvsNoise' PEdirn 'dirnPE_' motdir 'dirn'  int2str(amplitude) 'Pixel' int2str(freq) 'Cyc_Motion']; 
end
print(gcf,'-djpeg100',mot_str)
save(mot_str,'amplitudes', 'flip', 'freqs', 'motdir','MS', 'splitsignois', 'TRE_CS', 'TRE_GS', 'TRE_LGS', 'vers');
plotty = input('Enter 1 if you want to plot all images = ');
if plotty
%% finishing plots
% maxmg = max(abs(GS(:)));  %set a max for plots to avoid singularities
%originals
maxmg = 0.65*max(abs(compdata(:)));%liver phantom and sim data at least

% %% plot M, for simdata
% 
%     figure;colormap(gray);imagesc(mmat,[0 maxmg]);axis image;
%     set(gca,'xtick',[],'ytick',[]);
%      title(['Base M, Noise=',int2str(noisval)],'FontSize',15);
%     if numtis
%         format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
%             {[num2str(T1(1)),'/',num2str(T2(1))],...
%              [num2str(T1(2)),'/',num2str(T2(2))],...
%              [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
%             [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
% %         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
% %         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
%     else
% %         set(gca,'XTick',1:th_pts:nc)
% %         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
% %         set(gca,'YTick',nr/7:nr/7:nr)
% %         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
% %         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
% %         set(get(gca,'YLabel'),'Rotation',0.0)
%     set(gca,'XTick',1:8*th_pts:nc)
%     set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
%     set(gca,'YTick',nr/7:6*nr/7:nr)
%     set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
%     end
% 
%     if pr
%             M_str = 'M';
%             print(gcf,'-djpeg100',M_str)
%                 print(gcf,'-deps',M_str)
%                 saveas(gcf,M_str)
%             close
%     end
% 
% % %% plot CS gold, for simdata
%     figure;colormap(gray);imagesc(CSgold,[0 maxmg]);axis image;
%     set(gca,'xtick',[],'ytick',[]);
%      title(['CS N = \inf, Noise=',int2str(noisval)],'FontSize',15);
%     if numtis
%         format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
%             {[num2str(T1(1)),'/',num2str(T2(1))],...
%              [num2str(T1(2)),'/',num2str(T2(2))],...
%              [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
%             [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
% %         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
% %         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
%     else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
%     end
%     if pr
%         CSg_str = 'CS Gold Standard';
%         print(gcf,'-djpeg100',CSg_str)
% %             print(gcf,'-deps',CSg_str)
% %             saveas(gcf,CSg_str)
%         close
%     end
    
%iterate through each motion scenario
if vers == 5
%plot vers 5 stuff, including surfplot above
else
for sl = 1:parmf 
for se = 1:4,figure(100+sl);subplot(rr,cc,se);colormap(gray);imagesc(abs(compdata(:,:,se,sl)), [0 maxmg]);
title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end 
if vers == 0
    [~,h1]=suplabel(['Magnitude, ',int2str(amplitudes(sl)),'-pixel Motion, Noise=',int2str(noisval)],'t');
elseif vers == 4
    [~,h1]=suplabel(['Magnitude, ',int2str(freqs(sl)),'-cycle Motion, Noise=',int2str(noisval)],'t');
end

    for se = 1:4,figure(150+sl);subplot(rr,cc,se);colormap(gray);imagesc(angle(compdata(:,:,se,sl)),[-pi pi]);
    title([int2str(cyc(se)*360),'cyc']);axis image;set(gca,'xtick',[],'ytick',[]);end
if vers == 0
    [~,h2]=suplabel(['Magnitude, ',int2str(amplitudes(sl)),'-pixel Motion, Noise=',int2str(noisval)],'t');
elseif vers == 4
    [~,h1]=suplabel(['Magnitude, ',int2str(freqs(sl)),'-cycle Motion, Noise=',int2str(noisval)],'t');
end
%% just one
if vers == 2 % if you're running the noise plot, calculate the base images for the 2018 abstract.  Most of this loops contents are redundant, and can be safely deleted
    pr = 1;
    for se = 1:4
       figure;colormap(gray);imagesc(abs(compdata(:,:,se,6)));
     axis image;set(gca,'xtick',[],'ytick',[]);
         if pr
                PC_str = ['PC' int2str(360*cyc(se))];
                print(gcf,'-djpeg100',PC_str)
                    print(gcf,'-deps',PC_str)
                    saveas(gcf,PC_str)
                close
         end
    end
     figure;colormap(gray);imagesc(CS(:,:,6));axis image;
    set(gca,'xtick',[],'ytick',[]);
    if pr
            CS_str = ['CS'];
            print(gcf,'-djpeg100',CS_str)
                print(gcf,'-deps',CS_str)
                saveas(gcf,CS_str)
            close
    end
    figure;colormap(gray);imagesc(abs(LGS(:,:,6)));
    axis image;set(gca,'xtick',[],'ytick',[]);
    if pr       
         LGS_str = 'LGS';
         print(gcf,'-djpeg100',LGS_str)
    %          print(gcf,'-deps',LGS_str)
    %          saveas(gcf,LGS_str)
         close
    end
end
% title([int2str(amplitudes(sl)),'-pixel Motion, \Delta\theta=',int2str(cyc(se)*360),'\circ Original, {\alpha}=',...
%     int2str(flip),'\circ, TR=',num2str(TR),'ms']);
% 
%     if numtis
%         format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
%             {[num2str(T1(1)),'/',num2str(T2(1))],...
%              [num2str(T1(2)),'/',num2str(T2(2))],...
%              [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
%             [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
% %         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
% %         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
%     else
% %         set(gca,'XTick',1:th_pts:nc)
% %         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
% %         set(gca,'YTick',nr/7:nr/7:nr)
% %         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
% %         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
% %         set(get(gca,'YLabel'),'Rotation',0.0)
% %         for shft = 1:10
% %             text(2+16.1*(shft-1),137,'-\pi \rightarrow \pi','Color',[1 1 1]);
% %         end
%         set(gca,'XTick',1:8*th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
%         set(gca,'YTick',nr/7:6*nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
% %         for shft = 1:5%-2pi to 2pi
% %             text(4+32.2*(shft-1),137,'-2\pi\rightarrow2\pi','Color',[1 1 1],'FontSize',20,'FontWeight','bold');
% %         end
% 
%     end
% 
% if pr
%     COMP_str = [num2str(cyc(se)*2*180) 'PhaseCycle'];
%     print(gcf,'-djpeg100',COMP_str)
%             print(gcf,'-deps',COMP_str)
%             saveas(gcf,COMP_str)
%     close
% end
% end

%% CS
figure;colormap(gray);imagesc(CS(:,:,sl),[0 maxmg]);axis image;
set(gca,'xtick',[],'ytick',[]);
if vers == 0
    title(['CS, ',int2str(amplitudes(sl)),'-pix Motion, Noise = ',int2str(noisval),' TRE = ',num2str(TRE_CS(sl))],'FontSize',12);
elseif vers == 4
    title(['CS, ',int2str(freqs(sl)),'-cycle Motion, Noise = ',int2str(noisval),' TRE = ',num2str(TRE_CS(sl))],'FontSize',12);
end   

    if numtis
        format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
            {[num2str(T1(1)),'/',num2str(T2(1))],...
             [num2str(T1(2)),'/',num2str(T2(2))],...
             [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
            [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
    else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
    end

if pr
        CS_str = ['CS'];
        print(gcf,'-djpeg100',CS_str)
            print(gcf,'-deps',CS_str)
            saveas(gcf,CS_str)
        close
end
%% GS
figure;colormap(gray);imagesc(abs(GS(:,:,sl)),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);

if vers == 0
    title(['GS, ',int2str(amplitudes(sl)),'-pix Motion, Noise=',int2str(noisval),' TRE = ',num2str(TRE_GS(sl))],'FontSize',12);
elseif vers == 4
    title(['GS, ',int2str(freqs(sl)),'-cycle Motion, Noise = ',int2str(noisval),' TRE = ',num2str(TRE_GS(sl))],'FontSize',12);
end 
    if numtis
        format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
            {[num2str(T1(1)),'/',num2str(T2(1))],...
             [num2str(T1(2)),'/',num2str(T2(2))],...
             [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
            [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
    else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
    end
if pr       
    GS_str = ['GS'];
    print(gcf,'-djpeg100',GS_str)
        print(gcf,'-deps',GS_str)
        saveas(gcf,GS_str)
         close
end
%     format_ticks(gca,{'-2\pi','-\pi','0','\pi','2\pi'},[],[-2*pi -pi 0 pi 2*pi]);

    %
%Phase
%     figure;colormap(gray);imagesc(angle(GS(:,:,sl)));axis image;set(gca,'xtick',[],'ytick',[]);
%     title(['GS Phase, Noise=',int2str(noisval)],'FontSize',15);
%     %         print -djpeg100 gs_phase.jpg
%     %         crop ('gs_phase.jpg')
%% LGS = Linearized GS
   figure;colormap(gray);imagesc(abs(LGS(:,:,sl)),[0 maxmg]);
axis image;set(gca,'xtick',[],'ytick',[]);

if vers == 0
    title(['LGS, ',int2str(amplitudes(sl)),'-pix Motion, Noise=',int2str(noisval),' TRE = ',num2str(TRE_LGS(sl))],'FontSize',12);
elseif vers == 4
    title(['LGS, ',int2str(freqs(sl)),'-cycle Motion, Noise = ',int2str(noisval),' TRE = ',num2str(TRE_LGS(sl))],'FontSize',12);
end 
    if numtis
        format_ticks_image(gca,{'-2\pi','-\pi','0','\pi','2\pi'},...
            {[num2str(T1(1)),'/',num2str(T2(1))],...
             [num2str(T1(2)),'/',num2str(T2(2))],...
             [num2str(T1(3)),'/',num2str(T2(3))]},[1 1+nc/4 1+nc/2 1+nc*3/4 nc],...
            [rowrep/2 rowrep*3/2 rowrep*5/2],[],90); 
%         ylabel('T1/T2 (ms/ms)','Position',[-2 16],'FontSize',12);
%         xlabel('\theta off-resonance','Position',[33 33],'FontSize',12);
    else
%         set(gca,'XTick',1:th_pts:nc)
%         set(gca,'XTickLabel',amat(1,1:th_pts:nc));
%         set(gca,'YTick',nr/7:nr/7:nr)
%         set(gca,'YTickLabel',bmat(nr/7:nr/7:nr));
%         xlabel('E_2','FontSize',12);ylabel('b','FontSize',12);
%         set(get(gca,'YLabel'),'Rotation',0.0)
    set(gca,'XTick',1:8*th_pts:nc)
    set(gca,'XTickLabel',amat(1,1:8*th_pts:nc));
    set(gca,'YTick',nr/7:6*nr/7:nr)
    set(gca,'YTickLabel',bmat(nr/7:6*nr/7:nr),'FontSize',18);
    end
if pr       
     LGS_str = 'LGS';
     print(gcf,'-djpeg100',LGS_str)
%          print(gcf,'-deps',LGS_str)
%          saveas(gcf,LGS_str)
     close
end
end
end
%% print data and TRE
if pr  

        save('TRE.mat','TRE_GS','TRE_CS','TRE_LGS');

end
% save('snr.mat','snrCS','snrGS');
end

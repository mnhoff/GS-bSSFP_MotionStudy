% MotionGradientEnergyQuantification_MotionGSpaper 20210208 
% studies slices for Motion GS paper in 2021
% ISMRM2020 version just ran slice 167 for ISMRM abstract

% subtract motion corrupted data from non (complex and magnitude)
% compute gradient energy to quantify

clear all
format long g
sl = 1;
off = 1; %off TE data?
cyc = [0 1/4 2/4 3/4]';

rr = 1; cc = 4;
% slind = input(' Please enter the slice range of interest (just enter for all slices) = ');

%% Image load - ultimately I'd like to load all images in a .mat file automatically, but for now I formed the arrays manually for ISMRM
numfile = 4;
% dir = dirfile;
% files = dirfile(3:end);
% for fl = 1:length(files)
[flowname,pname] = uigetfile('*.mat','Select desired flow *.mat file');
load(flowname);
[noflowname,nopname] = uigetfile('*.mat','Select desired non-flow *.mat file');
load(noflowname);
% how many PCs suffer from motion?
pcflow = floworder{1}(1);


%% image masking
%     nslice = 2;
%     compdata = rebin_compslice(nslice,compdata);
[nr,nc,npc,nsl,nfl] = size(compmc_flow);

% slind = [33 34 57 58 75];% best coronals from 210131 data.  33 and 75 are so-so on the 4pc motions
% slind = 24:30;% best axials from 210131 data.  24 ok on 4pc, 26 so-so on
% slind = [1:10 80:90];% peripheral coronals from 210131 data. 87 all the way
% 1pc, 26-28 are best
% slind = 87;% best no-mask coronal
% slind = 45;%
% slind = 15:75;%
% slind = 11:15;
% slind = [31 34:36 39 40 47 49 61 70:72 76 87];
% slind = [34 35 70 71];
slind = 71;

if ~exist('slind')
    slind = [];
end
if ~isempty(slind)
    slrange = slind; numleft = length(slind);
    disp(['There are ',num2str(nsl),' total slices!! But I''m only keeping ',num2str(numleft),' slices between slice '...
        ,num2str(slrange(1)),' to ',num2str(slrange(end))]); 
else
    slrange = 1:nsl; 
end
compmc_flow = compmc_flow(1:nr,1:nc,:,slrange,:);
compmc_noflow = compmc_noflow(1:nr,1:nc,:,slrange);
CSMC_flow = CSMC_flow(1:nr,1:nc,slrange,:);
CSMC_noflow = CSMC_noflow(1:nr,1:nc,slrange);
LGS_flow = LGS_flow(1:nr,1:nc,slrange,:);
LGS_noflow = LGS_noflow(1:nr,1:nc,slrange);

%in-plane mask
load masks.mat
maskmat = masks(:,:,slrange);
%either mask the original data
% CSMC_noflow = CSMC_noflow.*maskmat;
% LGS_noflow = LGS_noflow.*maskmat;
% for fl = 1:nfl
%     CSMC_flow(:,:,:,fl) = CSMC_flow(:,:,:,fl).*maskmat;
%     LGS_flow(:,:,:,fl) = LGS_flow(:,:,:,fl).*maskmat;
% end
%or set aside masked results so that unmasked build-up images may still be
%shown

CSMC_noflow_msk = CSMC_noflow.*maskmat;
LGS_noflow_msk = LGS_noflow.*maskmat;
CSMC_flow_msk = CSMC_flow; LGS_flow_msk = LGS_flow;
for fl = 1:nfl
    CSMC_flow_msk(:,:,:,fl) = CSMC_flow(:,:,:,fl).*maskmat;
    LGS_flow_msk(:,:,:,fl) = LGS_flow(:,:,:,fl).*maskmat;
end
[nr,nc,nsl,nfl] = size(LGS_flow);


%% image differences
for fl = 1:nfl
%     CSmagdiff(:,:,:,fl) = abs(abs(CSMC_flow(:,:,:,fl))-abs(CSMC_noflow));
    CScompdiff(:,:,:,fl) = abs(CSMC_flow(:,:,:,fl) - CSMC_noflow);
    CScompdiff_msk(:,:,:,fl) = abs(CSMC_flow_msk(:,:,:,fl) - CSMC_noflow_msk);
%     LGSmagdiff(:,:,:,fl) = abs(abs(LGS_flow(:,:,:,fl))-abs(LGS_noflow));
    LGScompdiff(:,:,:,fl) = abs(LGS_flow(:,:,:,fl) - LGS_noflow);
    LGScompdiff_msk(:,:,:,fl) = abs(LGS_flow_msk(:,:,:,fl) - LGS_noflow_msk);
    for sl = 1:nsl
        %%gradient energies
%         CSmagGE(sl,fl) = sum(sum((CSmagdiff(:,:,sl,fl)).^2));
        CScompGE(sl,fl) = sum(sum((CScompdiff(:,:,sl,fl)).^2));
        CScompGE_msk(sl,fl) = sum(sum((CScompdiff_msk(:,:,sl,fl)).^2));
%         LGSmagGE(sl,fl) = sum(sum((LGSmagdiff(:,:,sl,fl)).^2));
        LGScompGE(sl,fl) = sum(sum((LGScompdiff(:,:,sl,fl)).^2));
        LGScompGE_msk(sl,fl) = sum(sum((LGScompdiff_msk(:,:,sl,fl)).^2));
%         %plot normalization
% %         summatrix = [CSmagdiff(:,:,sl,fl) CScompdiff(:,:,sl,fl) LGSmagdiff(:,:,sl,fl) LGScompdiff(:,:,sl,fl)];maxval = max(max(summatrix));
%         summatrix = [CScompdiff(:,:,sl,fl) LGScompdiff(:,:,sl,fl)];maxval = max(max(summatrix));
%         figure;
% %         subplot(2,2,1);colormap(gray);imagesc(CSmagdiff(:,:,sl,fl),[0 maxval]);axis image; title(['CSmagdif ',floworder{fl},'Sl',int2str(slrange(sl)),', GE = ',num2str(sprintf('%.2d',CSmagGE(sl,fl)))]);
%         subplot(1,2,1);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image; title(['CScompdif ',floworder{fl},'Sl',int2str(slrange(sl)),', GE = ',num2str(sprintf('%.2d',CScompGE(sl,fl)))]);
% %         subplot(2,2,3);colormap(gray);imagesc(LGSmagdiff(:,:,sl,fl),[0 maxval]);axis image; title(['GSmagdif ',floworder{fl},'Sl',int2str(slrange(sl)),', GE = ',num2str(sprintf('%.2d',LGSmagGE(sl,fl)))]);
%         subplot(1,2,2);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image; title(['GScompdif ',floworder{fl},'Sl',int2str(slrange(sl)),', GE = ',num2str(sprintf('%.2d',LGScompGE(sl,fl)))]);
    end
end

%% Some sample image plots - for now, normalize PCs, LGS, CS, and difference images separately
pr = input('Do you want to print images? (1 or 0)');
% pr = 1;
for sl = 1:nsl
    %no motion phase cycles

% figure;
% for pcc = 1:npc
%     subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_noflow(:,:,pcc,sl)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
% end
% [~,h1]=suplabel(['No-Flow Phase Cycles, Sl',int2str(slrange(sl))],'t');
%or separate images
sumPCs = squeeze(cat(4,abs(squeeze(compmc_noflow(:,:,:,sl))),squeeze(abs(compmc_flow(:,:,:,sl,:))))); maxval = 0.5*max(sumPCs(:));%normalize the PC images
for pcc = 1:npc
    figure;colormap(gray);imagesc(abs(compmc_noflow(:,:,pcc,sl)),[0 maxval]);axis image;axis off;title(['NoFlowPC' int2str(cyc(pcc)*360),'cyc']);
    if pr
        bSSFP_str = ['bSSFPpc' int2str(cyc(pcc)*360) '_noflow_sl' sprintf('%02i',slrange(sl)) ];
        print(gcf,'-djpeg100',bSSFP_str)
    %         close
    end
end
% base motion phase cycles
% figure;
% for pcc = 1:npc
%     subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,1)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
% end
% [~,h1]=suplabel([floworder{1} '-Flow Phase Cycles, Sl',int2str(slrange(sl))],'t');
%or separate images
% for pcc = 1:npc
%     figure;colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,1)),[0 maxval]);axis image;axis off;title([floworder{1} '-Flow PC' int2str(cyc(pcc)*360),'cyc']);
%     if pr
%         bSSFP1PC_str = ['bSSFPpc' int2str(cyc(pcc)*360) '- Slice' sprintf('%02i',slrange(sl)) '_' floworder{1} 'flow'];
%         print(gcf,'-djpeg100',bSSFP1PC_str)
%     %         close
%     end
% end
% figure;
% for pcc = 1:npc
%     subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,4)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
% end
% [~,h1]=suplabel([floworder{4} '-Flow Phase Cycles, Sl',int2str(slrange(sl))],'t');
%or separate images
for pcc = 1:npc
    figure;colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,2)),[0 maxval]);axis image;axis off;title([floworder{2} '-Flow PC' int2str(cyc(pcc)*360),'cyc']);
    if pr
    bSSFPallPC_str = ['bSSFPpc'  int2str(cyc(pcc)*360) '_' floworder{2} 'flow_Sl' sprintf('%02i', slrange(sl))];
        print(gcf,'-djpeg100',bSSFPallPC_str)
%         close
    end
end

%CS flow
% % figure;
% % subplot(1,4,1);colormap(gray);imagesc(abs(CSMC_noflow(:,:,sl)));axis image;title('CS - NoFlow');
% % subplot(1,4,2);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,1)));axis image;title(['CS - ',floworder{1}]);
% % subplot(1,4,3);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,2)));axis image;title(['CS - ',floworder{2}]);
% % subplot(1,4,4);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,3)));axis image;title(['CS - ',floworder{3}]);
% % [~,h1]=suplabel(['CS in ' pcflow '-PC Flow, Sl',int2str(slrange(sl ))],'t');
% % if pr
% %     baseCS1PC_str = ['Slice' sprintf('%02i',slrange(sl)) '_baseCS_' pcflow 'PCflow'];
% %     print(gcf,'-djpeg100',baseCS1PC_str)
% % %         close
% % end

%or separate images
sumCS = squeeze(cat(3,abs(squeeze(CSMC_noflow(:,:,sl))),squeeze(abs(CSMC_flow(:,:,sl,:))))); maxval = max(sumCS(:));%normalize the CS images
figure;colormap(gray);imagesc(abs(CSMC_noflow(:,:,sl)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slrange(sl)) ' - NoFlow']);
if pr
    CSnoflow_str = ['CS_NoFlow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',CSnoflow_str)
%         close
end
figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,1)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slrange(sl)) ' - ',floworder{1}]);
if pr
    CSflow1_str = ['CS_' floworder{1} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',CSflow1_str)
%         close
end
figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,2)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slrange(sl)) ' - ',floworder{2}]);
if pr
    CSflow2_str = ['CS_' floworder{2} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',CSflow2_str)
%         close
end
figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,3)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slrange(sl)) ' - ',floworder{3}]);
if pr
    CSflow3_str = ['CS_' floworder{3} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',CSflow3_str)
%         close
end
% % % % figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,4)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slrange(sl)) ' - ',floworder{4}]);
% % % % if pr
% % % %     CSflow4_str = ['CS_' floworder{4} '_Sl' sprintf('%02i',slrange(sl))];
% % % %     print(gcf,'-djpeg100',CSflow4_str)
% % % % %         close
% % % % end
% LGS flow
% % % % % % figure;
% % % % % % subplot(1,4,1);colormap(gray);imagesc(abs(LGS_noflow(:,:,sl)));axis image;title('LGS - NoFlow');
% % % % % % subplot(1,4,2);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,1)));axis image;title(['LGS - ',floworder{1}]);
% % % % % % subplot(1,4,3);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,2)));axis image;title(['LGS - ',floworder{2}]);
% % % % % % subplot(1,4,4);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,3)));axis image;title(['LGS - ',floworder{3}]);
% % % % % % [~,h1]=suplabel(['LGS in ' pcflow '-PC Flow, Sl',int2str(slrange(sl))],'t');
% % % % % % if pr
% % % % % %     baseLGS1PC_str = ['Slice' sprintf('%02i',slrange(sl)) '_baseLGS_' pcflow 'PCflow'];
% % % % % %     print(gcf,'-djpeg100',baseLGS1PC_str)
% % % % % % %         close
% % % % % % end

sumLGS = squeeze(cat(3,abs(squeeze(LGS_noflow(:,:,sl))),squeeze(abs(LGS_flow(:,:,sl,:))))); maxval = max(sumLGS(:));%normalize the LGS
figure;colormap(gray);imagesc(abs(LGS_noflow(:,:,sl)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slrange(sl)) ' - NoFlow']);
if pr
    LGSnoflow_str = [ 'LGS_NoFlow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',LGSnoflow_str)
%         close
end
figure;colormap(gray);imagesc(abs(LGS_flow(:,:,sl,1)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slrange(sl)) ' - ',floworder{1}]);
if pr
    LGSflow1_str = ['LGS_' floworder{1} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',LGSflow1_str)
%         close
end
figure;colormap(gray);imagesc(abs(LGS_flow(:,:,sl,2)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slrange(sl)) ' - ',floworder{2}]);
if pr
    LGSflow2_str = ['LGS_' floworder{2} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',LGSflow2_str)
%         close
end
figure;colormap(gray);imagesc(abs(LGS_flow(:,:,sl,3)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slrange(sl)) ' - ',floworder{3}]);
if pr
    LGSflow3_str = ['LGS_' floworder{3} 'flow_Sl' sprintf('%02i',slrange(sl))];
    print(gcf,'-djpeg100',LGSflow3_str)
%         close
end
%superplots - put all flows for a given slice in the same plot.  THis is a
%bit too busy, need to separate the allpc and 1pcsss
% % for sl = 1:nsl
% %     summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));
% %     figure;
% %     for fl = 1:nfl
% %         subplot(2,nfl,fl);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['CS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))]);
% %         subplot(2,nfl,fl+nfl);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['GS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))]);
% %     end
% %     [~,h1]=suplabel(['Recon Flow Artifact Difference Plots with Gradient Energy (GE), Sl',int2str(slrange(sl))],'t');
% %     if pr
% %         flow_str = ['Slice' sprintf('%02i',slrange(sl)) '_FlowArt'];
% %         print(gcf,'-djpeg100',flow_str)
% % %         close
% %     end
% % end
    

%superplots - put all flows for a given slice in the same plot
% % % %     summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));
% % % %     figure;
% % % %     for fl = 1:nfl
% % % %         subplot(2,3,fl);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['CS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))], 'FontSize', 7);
% % % %         subplot(2,3,fl+nfl);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['GS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))], 'FontSize', 7);
% % % %         %note: suplabel doesn't work with axis off!  so silly
% % % % %         [~,~]=suplabel(['Flow Artifact Difference Plots with Artifact Energy, Sl',int2str(slrange(sl))],'t');
% % % % %         if pr
% % % % %             flow1PC_str = ['Slice' sprintf('%02i',slrange(sl)) '_1PCFlowArt'];
% % % % %             print(gcf,'-djpeg100',flow1PC_str)
% % % % %     %         close
% % % % %         end   
% % % %     end
% % % %     [~,h1]=suplabel([pcflow '-PC Flow Artifact Difference Plots with Gradient Energy (GE), Sl',int2str(slrange(sl))],'t');
% % % %     if pr
% % % %         flow_str = ['Slice' sprintf('%02i',slrange(sl)) '_' pcflow 'pcFlowArt'];
% % % %         print(gcf,'-djpeg100',flow_str)
% % % %     end
% print the mask for the chosen slice
figure;colormap(gray);imagesc(maskmat(:,:,sl));axis image;axis off; title(['Flowmask_Sl' sprintf('%02i',slrange(sl))]);
Flowmask_str = ['Flowmask-Sl' sprintf('%02i',slrange(sl))];
            print(gcf,'-djpeg100',Flowmask_str)
%or separate difference images
    summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));%normalize the difference images
    for fl = 1:nfl
        figure;colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['CSdiff ',floworder{fl},', En = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))]);
        if pr
%             CSdiffFlow_str = ['CSdiff_Flow' floworder{fl} '_Sl' sprintf('%02i',slrange(sl)) '_ArtEn_' num2str(sprintf('%.1d',CScompGE(sl,fl)))];
            CSdiffFlow_str = ['CSdiff_' floworder{fl} 'flow_Sl' sprintf('%02i',slrange(sl))];
            print(gcf,'-djpeg100',CSdiffFlow_str)
        %         close
        end
        figure;colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['LGSdiff ',floworder{fl},', En = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))]);
        if pr
%             LGSdiffFlow_str = ['LGS_Flow' floworder{fl} '_Sl' sprintf('%02i',slrange(sl)) '_ArtEn_' num2str(sprintf('%.1d',LGScompGE(sl,fl)))];
            LGSdiffFlow_str = ['LGSdiff_' floworder{fl} 'flow_Sl' sprintf('%02i',slrange(sl))];%can't put decimals in file name
            print(gcf,'-djpeg100',LGSdiffFlow_str)
        %         close
        end
        figure;colormap(gray);imagesc(CScompdiff_msk(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['CSdiffMASK ',floworder{fl},', En = ',num2str(sprintf('%.1d',CScompGE_msk(sl,fl)))]);
        if pr
%             CSdiffFlow_str = ['CSdiff_Flow' floworder{fl} '_Sl' sprintf('%02i',slrange(sl)) '_ArtEn_' num2str(sprintf('%.1d',CScompGE(sl,fl)))];
            CSdiffFlowMSK_str = ['CSdiffMSK_' floworder{fl} 'flow_Sl' sprintf('%02i',slrange(sl))];
            print(gcf,'-djpeg100',CSdiffFlowMSK_str)
        %         close
        end
        figure;colormap(gray);imagesc(LGScompdiff_msk(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['LGSdiffMASK ',floworder{fl},', En = ',num2str(sprintf('%.1d',LGScompGE_msk(sl,fl)))]);
        if pr
%             LGSdiffFlow_str = ['LGS_Flow' floworder{fl} '_Sl' sprintf('%02i',slrange(sl)) '_ArtEn_' num2str(sprintf('%.1d',LGScompGE(sl,fl)))];
            LGSdiffFlowMSK_str = ['LGSdiffMSK_' floworder{fl} 'flow_Sl' sprintf('%02i',slrange(sl))];%can't put decimals in file name
            print(gcf,'-djpeg100',LGSdiffFlowMSK_str)
        %         close
        end
    end
end
     



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
% slind = [33 34 57 58 75];% best coronals from 210131 data.  33 and 75 are so-so on the 4pc motions
slind = 24:30;% best axials from 210131 data.  24 ok on 4pc, 26 so-so on
% 1pc, 26-28 are best

numfile = 4;
% dir = dirfile;
% files = dirfile(3:end);
% for fl = 1:length(files)
[flowname,pname] = uigetfile('*.mat','Select desired flow *.mat file');
load(flowname);
[noflowname,pname] = uigetfile('*.mat','Select desired non-flow *.mat file');
load(noflowname);

    %% reduce the dataset if it's too big
%     nslice = 2;
%     compdata = rebin_compslice(nslice,compdata);
[nr,nc,npc,nsl,nfl] = size(compmc_flow);

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

[nr,nc,~,nsl,~] = size(compmc_flow);

pr = input('Do you want to print images? (1 or 0)');
% pr = 1;
%% image differences
for fl = 1:nfl
%     CSmagdiff(:,:,:,fl) = abs(abs(CSMC_flow(:,:,:,fl))-abs(CSMC_noflow));
    CScompdiff(:,:,:,fl) = abs(CSMC_flow(:,:,:,fl) - CSMC_noflow);
%     LGSmagdiff(:,:,:,fl) = abs(abs(LGS_flow(:,:,:,fl))-abs(LGS_noflow));
    LGScompdiff(:,:,:,fl) = abs(LGS_flow(:,:,:,fl) - LGS_noflow);
    for sl = 1:nsl
        %%gradient energies
%         CSmagGE(sl,fl) = sum(sum((CSmagdiff(:,:,sl,fl)).^2));
        CScompGE(sl,fl) = sum(sum((CScompdiff(:,:,sl,fl)).^2));
%         LGSmagGE(sl,fl) = sum(sum((LGSmagdiff(:,:,sl,fl)).^2));
        LGScompGE(sl,fl) = sum(sum((LGScompdiff(:,:,sl,fl)).^2));
%         %plot normalization
% %         summatrix = [CSmagdiff(:,:,sl,fl) CScompdiff(:,:,sl,fl) LGSmagdiff(:,:,sl,fl) LGScompdiff(:,:,sl,fl)];maxval = max(max(summatrix));
%         summatrix = [CScompdiff(:,:,sl,fl) LGScompdiff(:,:,sl,fl)];maxval = max(max(summatrix));
%         figure;
% %         subplot(2,2,1);colormap(gray);imagesc(CSmagdiff(:,:,sl,fl),[0 maxval]);axis image; title(['CSmagdif ',floworder{fl},'Sl',int2str(slind(sl)),', GE = ',num2str(sprintf('%.2d',CSmagGE(sl,fl)))]);
%         subplot(1,2,1);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image; title(['CScompdif ',floworder{fl},'Sl',int2str(slind(sl)),', GE = ',num2str(sprintf('%.2d',CScompGE(sl,fl)))]);
% %         subplot(2,2,3);colormap(gray);imagesc(LGSmagdiff(:,:,sl,fl),[0 maxval]);axis image; title(['GSmagdif ',floworder{fl},'Sl',int2str(slind(sl)),', GE = ',num2str(sprintf('%.2d',LGSmagGE(sl,fl)))]);
%         subplot(1,2,2);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image; title(['GScompdif ',floworder{fl},'Sl',int2str(slind(sl)),', GE = ',num2str(sprintf('%.2d',LGScompGE(sl,fl)))]);
    end
end

%% Some sample image plots - for now, normalize PCs, LGS, CS, and difference images separately
%no motion phase cycles
% for sl = 7 % interested in 167. ATM, ONLY running one slice
for sl = 1:nsl
figure;
for pcc = 1:npc
    subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_noflow(:,:,pcc,sl)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
end
[~,h1]=suplabel(['No-Flow Phase Cycles, Sl',int2str(slind(sl))],'t');
%or separate images
% % % % % % % sumPCs = squeeze(cat(4,abs(squeeze(compmc_noflow(:,:,:,sl))),squeeze(abs(compmc_flow(:,:,:,sl,:))))); maxval = 0.5*max(sumPCs(:));%normalize the PC images
% % % % % % % for pcc = 1:npc
% % % % % % %     figure;colormap(gray);imagesc(abs(compmc_noflow(:,:,pcc,sl)),[0 maxval]);axis image;axis off;title(['NoFlowPC' int2str(cyc(pcc)*360),'cyc']);
% % % % % % %     if pr
% % % % % % %         bSSFP_str = ['bSSFPpc' int2str(cyc(pcc)*360) ' - Slice' sprintf('%02i',slind(sl)) '_noflow'];
% % % % % % %         print(gcf,'-djpeg100',bSSFP_str)
% % % % % % %     %         close
% % % % % % %     end
% % % % % % % end
% base motion phase cycles
figure;
for pcc = 1:npc
    subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,1)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
end
[~,h1]=suplabel([floworder{1} '-Flow Phase Cycles, Sl',int2str(slind(sl))],'t');
%or separate images
% % % % % % % for pcc = 1:npc
% % % % % % %     figure;colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,1)),[0 maxval]);axis image;axis off;title([floworder{1} '-Flow PC' int2str(cyc(pcc)*360),'cyc']);
% % % % % % %     if pr
% % % % % % %         bSSFP1PC_str = ['bSSFPpc' int2str(cyc(pcc)*360) '- Slice' sprintf('%02i',slind(sl)) '_' floworder{1} 'flow'];
% % % % % % %         print(gcf,'-djpeg100',bSSFP1PC_str)
% % % % % % %     %         close
% % % % % % %     end
% % % % % % % end
% figure;
% for pcc = 1:npc
%     subplot(1,numfile,pcc);colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,4)));axis image;title([int2str(cyc(pcc)*360),'cyc']);
% end
% [~,h1]=suplabel([floworder{4} '-Flow Phase Cycles, Sl',int2str(slind(sl))],'t');
%or separate images
% % % % for pcc = 1:npc
% % % %     figure;colormap(gray);imagesc(abs(compmc_flow(:,:,pcc,sl,4)),[0 maxval]);axis image;axis off;title([floworder{4} '-Flow PC' int2str(cyc(pcc)*360),'cyc']);
% % % %     if pr
% % % %         bSSFPallPC_str = ['bSSFPpc' int2str(cyc(pcc)*360) ' - Slice' sprintf('%02i',slind(sl)) '_' floworder{4} 'flow'];
% % % %         print(gcf,'-djpeg100',bSSFPallPC_str)
% % % % %         close
% % % %     end
% % % % end

%CS flow
figure;
subplot(1,4,1);colormap(gray);imagesc(abs(CSMC_noflow(:,:,sl)));axis image;title('CS - NoFlow');
subplot(1,4,2);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,1)));axis image;title(['CS - ',floworder{1}]);
subplot(1,4,3);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,2)));axis image;title(['CS - ',floworder{2}]);
subplot(1,4,4);colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,3)));axis image;title(['CS - ',floworder{3}]);
[~,h1]=suplabel(['CS in 4-PC Flow, Sl',int2str(slind(sl))],'t');
if pr
    baseCS1PC_str = ['Slice' sprintf('%02i',slind(sl)) '_baseCS1PCflow'];
    print(gcf,'-djpeg100',baseCS1PC_str)
%         close
end

%or separate images
% % % % % % % % sumCS = squeeze(cat(3,abs(squeeze(CSMC_noflow(:,:,sl))),squeeze(abs(CSMC_flow(:,:,sl,:))))); maxval = max(sumCS(:));%normalize the CS images
% % % % % % % % figure;colormap(gray);imagesc(abs(CSMC_noflow(:,:,sl)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slind(sl)) ' - NoFlow']);
% % % % % % % % if pr
% % % % % % % %     CSnoflow_str = ['CS_NoFlow_Slice' sprintf('%02i',slind(sl))];
% % % % % % % %     print(gcf,'-djpeg100',CSnoflow_str)
% % % % % % % % %         close
% % % % % % % % end
% % % % % % % % figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,1)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slind(sl)) ' - ',floworder{1}]);
% % % % % % % % if pr
% % % % % % % %     CSflow1_str = ['CS_' floworder{1} '_Slice' sprintf('%02i',slind(sl)) ];
% % % % % % % %     print(gcf,'-djpeg100',CSflow1_str)
% % % % % % % % %         close
% % % % % % % % end
% % % % figure;colormap(gray);imagesc(abs(CSMC_flow(:,:,sl,4)),[0 maxval]);axis image;axis off; title(['CS Sl' int2str(slind(sl)) ' - ',floworder{4}]);
% % % % if pr
% % % %     CSflow4_str = ['CS_' floworder{4} '_Slice' sprintf('%02i',slind(sl))];
% % % %     print(gcf,'-djpeg100',CSflow4_str)
% % % % %         close
% % % % end
% LGS flow
figure;
subplot(1,4,1);colormap(gray);imagesc(abs(LGS_noflow(:,:,sl)));axis image;title('LGS - NoFlow');
subplot(1,4,2);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,1)));axis image;title(['LGS - ',floworder{1}]);
subplot(1,4,3);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,2)));axis image;title(['LGS - ',floworder{2}]);
subplot(1,4,4);colormap(gray);imagesc(abs(LGS_flow(:,:,sl,3)));axis image;title(['LGS - ',floworder{3}]);
[~,h1]=suplabel(['LGS in 4-PC Flow, Sl',int2str(slind(sl))],'t');
if pr
    baseLGS1PC_str = ['Slice' sprintf('%02i',slind(sl)) '_baseLGS1PCflow'];
    print(gcf,'-djpeg100',baseLGS1PC_str)
%         close
end

% % % % sumLGS = squeeze(cat(3,abs(squeeze(LGS_noflow(:,:,sl))),squeeze(abs(LGS_flow(:,:,sl,:))))); maxval = max(sumLGS(:));%normalize the LGS
% % % % figure;colormap(gray);imagesc(abs(LGS_noflow(:,:,sl)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slind(sl)) ' - NoFlow']);
% % % % if pr
% % % %     LGSnoflow_str = [ 'LGS_NoFlow_Slice' sprintf('%02i',slind(sl))];
% % % %     print(gcf,'-djpeg100',LGSnoflow_str)
% % % % %         close
% % % % end
% % % % figure;colormap(gray);imagesc(abs(LGS_flow(:,:,sl,1)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slind(sl)) ' - ',floworder{1}]);
% % % % if pr
% % % %     LGSflow1_str = ['LGS_' floworder{1} '_Slice' sprintf('%02i',slind(sl)) ];
% % % %     print(gcf,'-djpeg100',LGSflow1_str)
% % % % %         close
% % % % end
% % % figure;colormap(gray);imagesc(abs(LGS_flow(:,:,sl,4)),[0 maxval]);axis image;axis off; title(['LGS Sl' int2str(slind(sl)) ' - ',floworder{4}]);
% % % if pr
% % %     LGSflow4_str = ['LGS_' floworder{4} '_Slice' sprintf('%02i',slind(sl)) ];
% % %     print(gcf,'-djpeg100',LGSflow4_str)
%         close
% % % end
%superplots - put all flows for a given slice in the same plot.  THis is a
%bit too busy, need to separate the allpc and 1pcsss
% % for sl = 1:nsl
% %     summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));
% %     figure;
% %     for fl = 1:nfl
% %         subplot(2,nfl,fl);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['CS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))]);
% %         subplot(2,nfl,fl+nfl);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['GS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))]);
% %     end
% %     [~,h1]=suplabel(['Recon Flow Artifact Difference Plots with Gradient Energy (GE), Sl',int2str(slind(sl))],'t');
% %     if pr
% %         flow_str = ['Slice' sprintf('%02i',slind(sl)) '_FlowArt'];
% %         print(gcf,'-djpeg100',flow_str)
% % %         close
% %     end
% % end
    

%superplots - put all flows for a given slice in the same plot
    summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));
    figure;
    for fl = 1:nfl
        subplot(2,3,fl);colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['CS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))], 'FontSize', 7);
        subplot(2,3,fl+nfl);colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;title(['GS ',floworder{fl},', GE = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))], 'FontSize', 7);
        %note: suplabel doesn't work with axis off!  so silly
%         [~,~]=suplabel(['Flow Artifact Difference Plots with Artifact Energy, Sl',int2str(slind(sl))],'t');
%         if pr
%             flow1PC_str = ['Slice' sprintf('%02i',slind(sl)) '_1PCFlowArt'];
%             print(gcf,'-djpeg100',flow1PC_str)
%     %         close
%         end   
    end
    [~,h1]=suplabel(['Recon Flow Artifact Difference Plots with Gradient Energy (GE), Sl',int2str(slind(sl))],'t');
    if pr
        flow_str = ['Slice' sprintf('%02i',slind(sl)) '_1pcFlowArt'];
        print(gcf,'-djpeg100',flow_str)
    end

%or separate difference images
% % % %     summatrix = [CScompdiff(:,:,sl,:) LGScompdiff(:,:,sl,:)];maxval = max(max(max(summatrix)));%normalize the difference images
% % % %     for fl = 1:nfl
% % % %         figure;colormap(gray);imagesc(CScompdiff(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['CS ',floworder{fl},', En = ',num2str(sprintf('%.1d',CScompGE(sl,fl)))]);
% % % %         if pr
% % % % %             CSdiffFlow_str = ['CSdiff_Flow' floworder{fl} '_Slice' sprintf('%02i',slind(sl)) '_ArtEn_' num2str(sprintf('%.1d',CScompGE(sl,fl)))];
% % % %             CSdiffFlow_str = ['CSdiff_Flow' floworder{fl} '_Slice' sprintf('%02i',slind(sl))];
% % % %             print(gcf,'-djpeg100',CSdiffFlow_str)
% % % %         %         close
% % % %         end
% % % %         figure;colormap(gray);imagesc(LGScompdiff(:,:,sl,fl),[0 maxval]);axis image;axis off; title(['LGS ',floworder{fl},', En = ',num2str(sprintf('%.1d',LGScompGE(sl,fl)))]);
% % % %         if pr
% % % % %             LGSdiffFlow_str = ['LGS_Flow' floworder{fl} '_Slice' sprintf('%02i',slind(sl)) '_ArtEn_' num2str(sprintf('%.1d',LGScompGE(sl,fl)))];
% % % %             LGSdiffFlow_str = ['LGSdiff_Flow' floworder{fl} '_Slice' sprintf('%02i',slind(sl))];%can't put decimals in file name
% % % %             print(gcf,'-djpeg100',LGSdiffFlow_str)
% % % %         %         close
% % % %         end
% % % %     end
end   



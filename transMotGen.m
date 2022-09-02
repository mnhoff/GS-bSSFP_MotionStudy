function [motdata,motKSP] = transMotGen(PEdirn,cdata,TR,varargin )
% mnh this uses the Volker MRI motion simulation as a shell, guts it, and 
% adds San's method first scripted in "simimshift.m":
% (1) Take the magnitude of your object as I0.
% (2) Move I0 according to your "motion schedule" to get (x-y-t) space data as a "movie".
% (3) 2DFT every movie frame to get to (Kx-Ky-t) space (using shift theorem, you can save some not so valuable processing time here:))
% (4) Sample the 3D (Kx-Ky-t) space obliquely into a 2D data sheet as you do in an MRI scan.
% (5) Perform the inverse 2DFT on the sampled data sheet to get the final image that is motion affected.

% pr = input('Do you want to print motion simulation images? (1 or 0)');
pr = 0;
%% (1) Base Image Magnitude
I0 = cdata;
% figure;imagesc(abs(I0));colormap(gray); title('Original'); axis image
   
%% (2a) Create motion schedule
[nr,nc,~] = size(cdata);
t  = (1:nr).' .* TR; 
motionType = varargin;% generate motion
pos = getMotion(t, motionType{:} ); %Volker's motion function

%% (2b) Create a movie {x-y-t} 
I0frame = zeros(nr,nc,nr);
for tt = 1:nr
    sx = pos(tt,1); sy = pos(tt,2);%ignore z-motion for now
%     I0frame(:,:,tt) = imtranslate(I0,[sx, sy]);% for some reason the x-horizontal imtranslate is not spatial-shift-invariant (non-linear)
    I0frame(:,:,tt) = circshift(I0,[round(sy), round(sx)]);
end

%% (3)2DFT each movie frame {x-y-t} -> {Kx-Ky-t}
ksp = fft2(I0frame);%may have to normalize and shift
shftKSP = fftshift(ksp);
%     I0framereturn = ifft2(ksp);%test that I get the same info when I ifft

% Visualize the shifted, FFTed movie frames, and a k-space surface plot 
% % % figure;
% % % for tt = 1: nr
% % %     pause(0.1)
% % %     im = subplot(1,3,1);imagesc(abs(I0frame(:,:,tt)));colormap(gray); axis image% 
% % %     title(['Image Mag, Frame # ',num2str(tt),', shift = ',num2str(pos(tt,2))]);
% % % % im1 = subplot(1,2,2);imagesc(log(abs(shftKSP(:,:,tt))), [0 10]); axis image
% % % im1 = subplot(1,3,2);imagesc(log(abs(shftKSP(:,:,tt)))); axis image
% % % colormap(gray)
% % % title(['K-space Mag, Frame # ',num2str(tt),', shift = ',num2str(pos(tt,2))]);
% % % im2 = subplot(1,3,3);imagesc(angle(shftKSP(:,:,tt))); axis image
% % % colormap(gray)
% % % title(['K-space Phase, Frame # ',num2str(tt)]);
% % %  
% % % end
% close

%% (4) Sample 3D {Kx-Ky-t} obliquely as real data would be read
% since each movie frame corresponds to a TR separation, taking 1 row per
% frame sequentially is appropriate
motKSP = zeros(nr,nc);
% for tt = nr:-1:1
if PEdirn == 'h'
    for tt = 1:nc
        motKSP(:,tt) = shftKSP(:,tt,tt);%asign the tt-th phase encode col of the tt-th frame
    end
else
    for tt = 1:nr
        motKSP(tt,:) = shftKSP(tt,:,tt);%asign the tt-th phase encode row of the tt-th frame
    end
end
% figure;
% imagesc(log(abs(motKSP)),[1 10]);colormap(gray); axis image
% title('Motion-corrupted k-space');
% if pr
%     print(gcf,'-djpeg100','MotCoruptKsp');
% end

%% (5) inverse 2DFT obliquely sampled data for motion-corrupted image
finalKSP = ifftshift(motKSP);% return the DC-peak to the k-space periphery prior to IFFT
motdata = ifft2(finalKSP);


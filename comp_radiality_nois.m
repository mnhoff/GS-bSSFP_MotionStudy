function [radtanmetric,rtmet, radtotnorm] = comp_radiality_nois(noisonly_rot) 
% compute noise radiality metric, but sum noise of all 4 cycles BEFOREhand
% Radiality is along the real axis,tangentiality along the imaginary
[nr,nc,npc] = size(noisonly_rot);
rtmet = zeros(nr,nc,npc);
radtanmetric = zeros(nr,nc); radtotnorm = zeros(nr,nc);
for rr = 1:nr
    for cc = 1:nc
        sumrad = 0; sumtan = 0;
        for pc = 1:npc
            radnois = abs(real(noisonly_rot(rr,cc,pc)));
            tannois = abs(imag(noisonly_rot(rr,cc,pc)));
            rtmet(rr,cc,pc) = (radnois-tannois)/(radnois+tannois);
            sumrad = sumrad + radnois; sumtan = sumtan + tannois;
        end
        radtotnorm(rr,cc) = (sumrad-sumtan)/(sumrad+sumtan);
        radtanmetric(rr,cc) = mean(rtmet(rr,cc,:));
    end
end
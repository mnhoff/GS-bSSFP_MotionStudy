function [noisonly_rot] = noisrotate(noisonly,bisang13,bisang24)%140901 mnh
% function rotates the pure noise applied to simulated data such that the
% 'radial' noise (along a cross-solution bisector) is horizontal, and the
% 'tangential' noise is vertical.  The bisector angles are positive when
% they need a clockwise rotation, and negative if they need
% counterclockwise. Thus the rotation is negative (clockwise)

[nr,nc,npc] = size(noisonly);
noisonly_rot = zeros(nr,nc,npc);
if npc ~= 4
    error('the number of phase cycles is not 4! Rewrite your code ;)')
end
for rr = 1:nr
    for cc = 1:nc
        noisonly_rot(rr,cc,1) = noisonly(rr,cc,1)*exp(-1i.*bisang13(rr,cc));
        noisonly_rot(rr,cc,2) = noisonly(rr,cc,2)*exp(-1i.*bisang24(rr,cc));
        noisonly_rot(rr,cc,3) = noisonly(rr,cc,3)*exp(-1i.*bisang13(rr,cc));
        noisonly_rot(rr,cc,4) = noisonly(rr,cc,4)*exp(-1i.*bisang24(rr,cc));
    end
end
        
        
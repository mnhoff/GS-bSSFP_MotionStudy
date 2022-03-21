function [bisang13,bisang24] = calc_spoke_angle(cdata)%% 14/08/01 input noiseless bSSFP data
% calculate the angle of each bisecting line relative to the horizontal
% real axis
[nr,nc,~] = size(cdata);
bisang13 = zeros(nr,nc);
bisang24 = zeros(nr,nc);
for r = 1:nr
    for c = 1:nc
        dat = cdata(r,c,:);
        bisang13(r,c) = atan((imag(dat(3)) - imag(dat(1)))/(real(dat(3)) - real(dat(1))));
        bisang24(r,c) = atan((imag(dat(4)) - imag(dat(2)))/(real(dat(4)) - real(dat(2))));
    end
end             
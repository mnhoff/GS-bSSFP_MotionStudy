function TRE = totrelerrorTRE(data,gold)% 111105 calculate the total 
% relative error TRE of an image from a gold standard as in eq [14] from 
% SPEED-ACE
[nrg,ncg] = size(gold);         [nr,nc] = size(data);
if nrg ~= nr || ncg ~= nc
    error('The image and it''s gold standard are different sizes!');
end
data = abs(data(:));            gold = abs(gold(:));
dif = abs(data) - abs(gold);
SOSres = sum(dif.^2);
% if ~all(gold)
%     TRE = sqrt(SOSres)/sum(SOSres);
% else
%     TRE = sqrt(SOSres)/sum(gold);
% end
TRE = sqrt(SOSres)/sum(gold);
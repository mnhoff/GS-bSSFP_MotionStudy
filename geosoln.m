function [master,regmap] = geosoln(compdata)%151203 mnh
% compute 4pt geometric cross-solution (GS)
[nr,nc,npc,nsl,nph,nch] = size(compdata);
if npc == 4
    master = zeros(nr,nc,nsl,nph,nch);
    regmap = zeros(nr,nc,nsl,nph,nch);
    for ch = 1:nch
        for np = 1:nph
            for s = 1:nsl
                for r = 1:nr
                    for c = 1:nc
                        cdat = compdata(r,c,:,s,np,ch); rdat = real(cdat); idat = imag(cdat);
                        numer = (rdat(1)*idat(3)-rdat(3)*idat(1))*(cdat(2)-cdat(4)) ...
                            - (rdat(2)*idat(4)-rdat(4)*idat(2))*(cdat(1)-cdat(3));
                        denom = rdat(1)*idat(2)+rdat(2)*idat(3)+rdat(3)*idat(4)+rdat(4)*idat(1)- ...
                            (rdat(1)*idat(4)+rdat(4)*idat(3)+rdat(3)*idat(2)+rdat(2)*idat(1));
                        master(r,c,s,np,ch) = numer/denom;

                        % regularize with CS
                        if denom == 0 || (abs(master(r,c,s,np,ch)) > max(abs(cdat)))
                            master(r,c,s,np,ch) = mean(cdat);
                            regmap(r,c,s,np,ch) = 1;% note when there is a regularized pixel
                        end
                    end
                end
            end
        end
    end
else
    error('Your dataset in incorrectly formatted!');
end
        
function [DTCmask] = DetectDTC(rPla, iminfo, gmask, pix, mupr)

%%% DTCmask: 
%%%        | 1: mask image on each z-plane 
sizel = iminfo(4); 
t1 = cell(sizel,1);
DTCmask = cell(iminfo(4),2);
fint = cell(iminfo(4),1);
se = strel('disk', pix*10);
gmask1 = imdilate(gmask, se);
% totInt = zeros(iminfo(4),1);
mofm = zeros(iminfo(4),1);
sofs = zeros(iminfo(4),1);
mask1 = cell(sizel,1);
% mupr = 1;

% make general threshold
multipr = 1;
parfor i = 1:sizel
    t1{i} = rPla{i,1};
    
    % make a binary image where background signal is.
    thrs = adaptthresh(t1{i});
    mask1{i} = imbinarize(t1{i}, thrs);

    % calculate background level only inside the gonad using gonadal mask.
    % calculate background for individual z planes and save in 'mofm'
    bAll = t1{i}(mask1{i} == 0 & gmask1 > 0);
    sAll = t1{i}(mask1{i} > 0);
    sAlls = sortrows(sAll, 'descend');
    sAlln = sAlls(round(size(sAlls,1)*1/4):round(size(sAlls,1)*3/4));
    mofm(i) = mean(bAll) * multipr;   
    sofs(i) = mean(sAlln);
end

% remove NaN from 'mofm' to calculate 'mofm'
tempM = mean(mofm(~isnan(mofm)));
mofm(isnan(mofm)) = tempM;
mofm = mofm * mupr;
range = [round(length(mofm)/4) ceil(length(mofm)*3/4)];
mom = mean(mofm(range(1):range(2)));
%%% normalize background for top and bottom slides
mofm(mofm < mom * 0.75) = mom;


for i = 1:sizel
%     t2 = medfilt2(t1{i}, [pix*2 pix*2]);
    t2 = t1{i} .* uint16(mask1{i});
    t2(t2 < sofs(i)* 0.75 ) = 0;
    t3 = bwareaopen(logical(t2), round(pix*2));
    t4 = imclose(t3, true(pix*2));
    
    DTCmask{i,1} = t4;
    DTCmask{i,2} = (t1{i} - mofm(i)) .* uint16(t4);
%     totInt(i) = sum(sum(uint8(t5) .* fint{i}));
    fprintf ('\nProcessing %d-th z-slice of the DTC image', i);
end


%%% visaul DTC
% for i = 1:iminfo(4)
%     figure, imshow(DTCmask{i});
%     pause
%     close all
%     fprintf('\n%d-th zplane.', i);
% end


%%% reconstruct DTC in 3D.
% blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
%       | 6th: matching nuc #(row # in 'nuc')
% tfDTC = cat(3,DTCmask{:});
% temp = bwconncomp(tfDTC, 26);
% connDTC = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
% blDTC = struct2cell(connDTC)';


















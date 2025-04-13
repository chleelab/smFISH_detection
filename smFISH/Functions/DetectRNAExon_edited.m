function [DETrna, derna, mom] = DetectRNAExon_edited(rPla, pixr, iminfo, f, numOfImg, g, thresForExon, gmask)
% This function detects cytoplasmic mRNA spots (transcription sites). The
% output is RNA coordinates, counts and intensity

sizel = iminfo(4);
fint = cell(sizel,1);  
t1 = cell(sizel,1); % cropped RNA images only in the nuclei.
mofm = zeros(sizel,1);
Pmofm = zeros(sizel,1);
pMask = cell(sizel,1);
se = strel('disk', pixr);
se2 = strel('disk', pixr*20);
se3 = strel('disk', pixr*4);
mupr = 1;
mupr2 = thresForExon; % multiplier to define the background level (becomes higher when bg is high)     
if thresForExon <= 0
    mupr2 = 0.001;
end
gmask2 = imerode(gmask, se2);

% ================ %%% sharing for RNA detection ========================= 
% make general threshold
multipr = 0.5;
parfor i = 1:sizel
    t1{i} = rPla{i,1};
    
    % make a binary image where background signal is.
    thrs = multithresh(t1{i});
    mask1 = imquantize(t1{i}, thrs);
    minV = min(min(mask1));
    mask1 = mask1 - minV;
    mask1(mask1 > 0) = 1;
    
    % calculate background level only inside the gonad using gonadal mask.
    % calculate background for individual z planes and save in 'mofm'
    tAll = t1{i}(mask1 > 0 & gmask2 > 0);
    mofm(i) = mean(tAll) * multipr;     
end

% remove NaN from 'mofm' to calculate 'mofm'
tempM = mean(mofm(~isnan(mofm)));
mofm(isnan(mofm)) = tempM;
mofm = mofm * mupr;
range = [round(length(mofm)/4) ceil(length(mofm)*3/4)];
mom = mean(mofm(range(1):range(2)));

%%% normalize background for top and bottom slides
mofm(mofm < mom * 0.85) = mom;
mofm = mofm * 1.5;    %%% <--- IMPORTANT for the function  "FastPeakFind", a.k.a. local peak detection!! See line 57!!!.



%%% threshold images and detect RNA spots
parfor i = 1:sizel 
% =============== peak detection 
    [A, B] = FastPeakFind(t1{i}, mofm(i));
    A = [A(1:2:end) A(2:2:end)];
        
    %%% Overall signal should be 25% brighter than the overall background
    %%% This also looks at individual detected peaks and compare to its surrounding.
    %%% That ratio should be larger than 1.1
    temp = mean(imdilate(t1{i}(B == 1),se));
    f1 = imdilate(B, se3);
    bgTMP = mean(t1{i}(f1 == 0 & gmask2 == 1));
    if isnan(temp)
        Pmofm(i) = 0;
    else
        Pmofm(i) = temp;
    end
    
    if Pmofm(i)/bgTMP < 1.4 % cutoff for overall dots-bg ratio
        B(A(:,2), A(:,1)) = 0;
    else
        ratioSB = zeros(size(A,1),2);
        ratioSoverall = zeros(size(A,1),2);  %%% signal-to-overall background ratio.
        A(:,3) = 0;
        for j = 1:size(A,1)
            rng1 = pixr*2;
            rng2 = pixr*10;

            Xrange1 = A(j,1)-rng1:A(j,1)+rng1 ;
            Yrange1 = A(j,2)-rng1:A(j,2)+rng1;
            Xrange1 = Xrange1(Xrange1 > 0 & Xrange1 <= iminfo(2));
            Yrange1 = Yrange1(Yrange1 > 0 & Yrange1 <= iminfo(3));
            sigMean1 = t1{i}( Yrange1 , Xrange1 );
            sigMean2 = sigMean1(sigMean1 > 0);
            sortN1 = sort(sigMean2(:), 'descend');
            sizeS = round(length(sortN1)/2);
            sigMean2 = mean(sortN1(1:sizeS));

            Xrange2 = A(j,1)-rng2:A(j,1)+rng2;
            Yrange2 = A(j,2)-rng2:A(j,2)+rng2;
            Xrange2 = Xrange2(Xrange2 > 0 & Xrange2 <= iminfo(2));
            Yrange2 = Yrange2(Yrange2 > 0 & Yrange2 <= iminfo(3));
            bgMean1 = t1{i}( Yrange2 , Xrange2 );
            bgMean1(rng2-rng1*2+1:rng1+rng2+1, rng2-rng1*2+1:rng1+rng2+1) = 0;
            bgMean2 = bgMean1(bgMean1 > 0);
            sortN2 = sort(bgMean2(:));
            sizeB1 = round(length(sortN2)/4);
            sizeB2 = round(length(sortN2)*3/4);
            bgMean3 = mean(sortN2(sizeB1:sizeB2));
            ratioSB(j) = sigMean2/bgMean3;
            ratioSoverall(j,1) = sigMean2/bgTMP;

            if ratioSB(j) < 1.15 + (mupr2-1)/2 || ratioSoverall(j,1) < 1.35 + (mupr2-1)/4  % the cutoff signal-local bg ratio of each mRNA (3x3 mRNA vs. 7x7 bg region.
                B(A(j,2), A(j,1)) = 0;
                A(j,3) = 999;
                ratioSB(j,2) = 999;
                ratioSoverall(j,2) = 999;
            end
        end
        A(A(:,3) == 999,:) = [];
        ratioSB(ratioSB(:,2) == 999,:) = [];
        ratioSoverall(ratioSoverall(:,2) == 999,:) = [];
        A(:,3) = [];
    end
%     B = imdilate(B, se);
   
%%% Quick visual -----------------------
    figure, imshow(t1{i}*30)
    hold on
    plot(A(:,1), A(:,2), 'r+');
    numbn = 1:size(A,1);
    for j = 1:size(A,1)
        text(A(j,1), A(j,2), num2str(ratioSB(j)), 'color', 'c');
%         text(A(j,1), A(j,2), num2str(ratioSoverall(j)), 'color', 'c');
%         text(A(j,1), A(j,2), num2str(numbn(j)), 'color', 'c');           %%Signal ID
    end
%%%-------------------------------------

    pMask{i} = B;
    fprintf('\nCh%d: %d/%d images,... %d/%d z-planes.', g, f, numOfImg, i, sizel);
end
fprintf('\n');

%%% Normalize z-slices by mRNA intensities to compensate photobleach or poor read inside sample.
Pmmofm = mean(Pmofm(Pmofm > 0 ));
Pmofm(Pmofm < 1 ) = Pmmofm;

for i = 1:sizel
    timesN = Pmmofm  / Pmofm(i);
    fint{i} = t1{i} * timesN ;
end



%%%%%%%%%%%%%%%%%%% 3D reconstitution %%%%%%%%%%%%%%%%%%%%%%%%%
%%% For mRNA (cytoplasmic spots)
for i=1:size(pMask,1)
    if isempty(pMask{i})
        pMask{i} = zeros(iminfo(3),iminfo(2));
    end
end

tfrna = cat(3,pMask{:});
temp = bwconncomp(tfrna, 26);
conrna = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
derna = struct2cell(conrna)';

% calculates intensity of detected blobs from original images.
iint = zeros(1,length(derna(:,1)));
for i=1:length(derna(:,1))
    iint(i) = 0;
    for j=1:derna{i,3}(6)
        cutimg = fint{round(derna{i,3}(3)),1};
        cutimg = cutimg(round(derna{i,3}(2)):round(derna{i,3}(2))...
            +derna{i,3}(5)-1, round(derna{i,3}(1)): round(derna{i,3}(1))+derna{i,3}(4)-1);

        iint(i) = iint(i) + sum(sum(derna{i,4}(:,:,j) .* double(cutimg)));
    end
end

derna = [derna num2cell(iint')];


if  isempty(derna) == 0
    if isempty(derna{1,1}) == 0
        % 'DETrnaP': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %           | total intensity | x centroid | y centroid | z centroid.
        DETrnaP = zeros(length(derna(:,1)),1);
        DETrnaP(:,1) = 1:length(derna(:,1));
        temp = cell2mat(derna(:,3));
        DETrnaP(:,2:4) = temp(:,4:6);
        temp = cell2mat(derna(:,2));
        DETrnaP(:,7:9) = temp(:,1:3);

        DETrnaP(:,6) = cell2mat(derna(:,5));
        DETrnaP(:,5) = cell2mat(derna(:,1));

        %%% remove false-positively detected blobs
        if ~isempty(DETrnaP)
            % get mean of obj. size (pixels) and intensity (normalized: ind. int. - mean)
            mint = mean(DETrnaP(:,6));
            mpix = mean(DETrnaP(:,5));
        end

        % Remove objects detected less than XX z-planes
        % Remove dim/small objects .
        if ~isempty(DETrnaP)
            XX = 2;
            temp = DETrnaP(:,4) < XX;
            DETrnaP(temp,:) = [];
        end
            
        if ~isempty(DETrnaP)
            temp = DETrnaP(:,4) < XX & DETrnaP(:,6) < mint*0.20 ;
            DETrnaP(temp,:) = [];
        end

        if isempty(DETrnaP) == 0
            temp = DETrnaP(:,4) < XX & DETrnaP(:,5) < mpix*0.20 ;
            DETrnaP(temp,:) = [];
        end
    else
        DETrnaP = [];
    end
else
    DETrnaP = [];
end

% DETrna : |1-3: x,y,z-coordinates | 4: zero | 5: vol ratio to mean vol per spot |
%          |6: total intensity| 7: z-plane of the center
%          Column 4 stays zero if mRNA is not matched to a nucleus (mRNA in rachis).
if ~isempty(DETrnaP)

    DETrna = DETrnaP(:,7:9);
    DETrna(:,3) = DETrna(:,3) * iminfo(7) / iminfo(6);
    DETrna(:,4) = 0;
    DETrna(:,5) = floor(DETrnaP(:,5) / mean(DETrnaP(:,5)) *2/3);
    DETrna(DETrna(:,5) == 0,5) = 1;
    DETrna(:,6) = DETrnaP(:,6);
    DETrna(:,7) = DETrnaP(:,9);

    % remove mRNA with intensity 0
    DETrna(DETrna(:,6) <= 0,:) = [];
else
    DETrna = [];
end



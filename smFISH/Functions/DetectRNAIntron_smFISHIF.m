function [blobrna, blrna, thrrna, mom] = DetectRNAIntron_smFISHIF(rPla, pixr, iminfo, f, numOfImg, g, thresForIntron, gmask)
% This function detects intronic RNA spots (transcription sites). The
% output is RNA coordinates, counts and intensity


showResult =      0        ;  %%% 1 to show, 0 to omit


sizel = iminfo(4); 
t1 = cell(sizel,1); % cropped RNA images only in the nuclei.
mofm = zeros(sizel,1);
thrrna = cell( sizel, 1);
Pmofm = zeros(sizel,1);
% sen = strel('disk', pixr*8); 
se = strel('disk', pixr*1); 
se2 = strel('disk', pixr*20);
se3 = strel('disk', pixr*2); 
mupr = 1;
mupr2 = thresForIntron + 0; % multiplier to define the background level (becomes higher when bg is high)
if thresForIntron <= 0
    mupr2 = 0.001;
end
gmask2 = imerode(gmask, se2);

% make general threshold
multipr = 1;
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
mofm = mofm * 0.3;    %%% <--- IMPORTANT for the function  "FastPeakFind", a.k.a. local peak detection!! See line 57!!!.



%%%%=================================================================================
%%%%=================================================================================



%%% threshold images and detect RNA spots
parfor i = 1:sizel    %%% i is for the corresponding i-th z-plane.
%%%----------------- peak detection method (1/3) ------------------
    [A, B] = FastPeakFind(t1{i}, mofm(i));
    A = [A(1:2:end) A(2:2:end)];
    
    
%%%------------------ Quick visual -----------------------
%     figure, imshow(t1{i}*20)
%     figure, imshow(t1{i}*20)
%     hold on
%     plot(A(:,1), A(:,2), 'r+');
%%%--------------------------------------------------------


    if ~isempty(A)
        %%% Overall signal should be 25% brighter than the overall background
        %%% This also looks at individual detected peaks and compare to its surrounding.
        %%% That ratio should be larger than 1.1
        temp = mean(t1{i}(B == 1));
        f1 = imdilate(B, se3);
        bgTMP = mean(t1{i}(f1 == 0 & gmask2 == 1));
        if isnan(temp)
            Pmofm(i) = 0;
        else
            Pmofm(i) = temp;
        end
        A(:,3) = 0;
        
        if Pmofm(i)/bgTMP < 0.5  % cutoff for overall dots-bg ratio
            B(A(:,2), A(:,1)) = 0;
            A(:,3) = 999;
        else
            ratioSB = zeros(size(A,1),2);   %%% local signal-to-noise ratio.
            ratioSoverall = zeros(size(A,1),2);  %%% signal-to-overall background ratio.
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
                ratioSB(j,1) = sigMean2/bgMean3;
                ratioSoverall(j,1) = sigMean2/bgTMP;
%                 fprintf('\n%.2f for %d', ratioSB(j), j)
                
                if ratioSB(j,1) < 1.05 + (mupr2-1)/2 || ratioSoverall(j,1) < 0.6 + (mupr2-1)/2/2  % the cutoff signal-local bg ratio of each ATS (3x3 ATS vs. 7x7 bg region).
                    B(A(j,2), A(j,1)) = 0;
                    A(j,3) = 999;
                    ratioSB(j,2) = 999;
                    ratioSoverall(j,2) = 999;
                end
            end
            A(A(:,3) == 999,:) = [];
            ratioSB(ratioSB(:,2) == 999,:) = [];
            ratioSoverall(ratioSoverall(:,2) == 999,:) = [];
        end
        
        
%%%------------------ Quick visual -----------------------
%         figure, imshow(t1{i}*40)
%         figure, imshow(t1{i}*40)
%         hold on
%         plot(A(:,1), A(:,2), 'r+');
%         numbn = 1:size(A,1);
%         for j = 1:size(A,1)
% %             text(A(j,1), A(j,2), num2str(ratioSB(j)), 'color', 'c');       %%Local Ratio
% %             text(A(j,1), A(j,2), num2str(ratioSoverall(j)), 'color', 'c'); %%Overall ratio
%             text(A(j,1), A(j,2), num2str(numbn(j)), 'color', 'c');           %%Signal ID
%         end
%%%--------------------------------------------------------
      
        if ~isempty(A)
            % remove peaks too close
            A(:,3) = 0;
            for j = 1:length(A(:,1))
                for k = j+1:length(A(:,1))
                    dist = sqrt((A(k,1)-A(j,1))^2 + (A(k,2)-A(j,2))^2);
                    if dist < pixr*2
                        A(j,3) = 999;
                        B(A(j,2),A(j,1)) = 0;
                    end
                end
            end
            A(A(:,3) == 999,:) = [];
            A(:,3) = [];

%%%---------- Image segmentation (bwconncomp) method (2/3) ---------
%             t2 = t1{i} - mofm(i);
%             t2(t2 < 0 ) = 0;
% 
%             loopn = 2;
%             for j = 1:loopn
%                 t2 = medfilt2(t2, [pixr pixr]);
%                 tt = reshape(t2,[],1);
%                 tt = sort(tt);
%                 ty = tt(tt ~= 0);
%                 tmofm = mean(mean(ty(20:end)));
%                 t2(t2 < tmofm) = 0;
%             end
% 
%             thrs = multithresh(t2);
%             mask1 = imquantize(t2, thrs);
%             minV = min(min(mask1));
%             mask1 = mask1 - minV;
%             mask1(mask1 > 0) = 1;
%             t5 = bwareaopen(mask1, ceil(pixr*2));
%             t6 = imdilate(t5, se);
% 
% %%%------ Cross validate detected spots (bwconncomp + peakfinder) ----
% %%% This removes bwconncomp-detected dots that are not found by peakfinder.
%             A(:,3) = 0;
%             for j=1:size(A,1)
%                 if t6(A(j,2),A(j,1)) == 0
%                     A(j,3) = 999;
%                     B(A(j,2),A(j,1)) = 0;
%                 end
%             end
%             A(A(:,3) == 999,:) = [];
%             
            
          
%=========== visual detected RNAs =======
%             figure,imshow(t1{i}*100);
%             hold on
%             plot(A(:,1), A(:,2), 'r+');
%             pause
%             close all
% ---------------------------------------

            if ~isempty(A)
                A(:,3) = 0;

            %%%-------- signal-to-local background ratio method (3/3) ----------
            % get intensity of spots and background to get the s:n ratio.
                sizC = pixr*2-1; % radius of the center region to get intensity
                cal = zeros(size(A,1), 4);
                f1 = t1{i};
                f1(gmask==0) = 0;

                if ~isempty(A) && sum(sum(A)) ~= 0
                    for winSize = 2:4
                        for j=1:size(A,1)
                            % 'A' col: | x-coor | y-coor | center_intensity | backgd2 int | bkgd3 int | bkgd4 int |.
                            sizB = pixr*2*winSize; 
                            xb = 1:sizB*2+1;

                            [x, y] = meshgrid(xb,xb);
                            Rb = sqrt((x-sizB-1).^2 + (y-sizB-1).^2); % center circle

                            Rc = Rb;
                            Rc(Rc>sizC) = 0;
                            Rc(sizB+1,sizB+1) = 1;
                            Rb(Rb>sizB) = 0;
                            Rb(Rb<sizC+2) = 0;

                            LB = A(j,2) - sizB; % left boundary of the ROI for scanning
                            if LB < 1
                                Rb(:,1:1-LB)=[];
                                Rc(:,1:1-LB)=[];
                                LB = 1;
                            end

                            RB = A(j,2) + sizB;
                            if RB > iminfo(3)
                                Rb(:,end-RB+iminfo(3)+1:end)=[];
                                Rc(:,end-RB+iminfo(3)+1:end)=[];
                                RB = iminfo(3);
                            end

                            UB = A(j,1) - sizB;
                            if UB < 1
                                Rb(1:1-UB,:)=[];
                                Rc(1:1-UB,:)=[];
                                UB = 1;
                            end

                            DB = A(j,1) + sizB;
                            if DB > iminfo(2)
                                Rb(end-DB+iminfo(2)+1:end,:)=[];
                                Rc(end-DB+iminfo(2)+1:end,:)=[];
                                DB = iminfo(2);
                            end

                            cutcen = f1(LB:RB,UB:DB);
                            centr = cutcen(Rc > 0);
                            if sum(centr) ~= 0
                                centr(centr == 0) = [];
                            end

                            cutbck = f1(LB:RB,UB:DB);
                            surrd = cutbck(Rb > 0);
                            if sum(surrd) ~= 0
                                surrd(surrd == 0) = [];
                            end
                            
                            sizeC = round(length(centr)/2);
                            sortC1 = sort(centr, 'descend');
                            ctr = sortC1(1:sizeC);
                            if sum(surrd) ~= 0
                                ctr = ctr(ctr > 0);
                            end
                            ctrMean = mean(ctr);
                            
                            sizeS1 = ceil(length(surrd)/4);
                            sizeS2 = round(length(surrd)*3/4);
                            sortS1 = sort(surrd);
                            srd = sortS1(sizeS1:sizeS2);
                            if sum(srd) ~= 0
                                srd = srd(srd > 0);
                            end
                            srdMean = mean(srd);

                            % subtract detected spots that are dim (not significant spots).
                            if ctrMean ~= 0 && srdMean ~= 0
                                [~,cal(j,1)] = ttest2(double(ctr), double(srd));
                                cal(j,2) = ctrMean/srdMean; %ratio to surrounding bg.
                                cal(j,3) = ctrMean/mofm(i); % ratio to mean background (fixed)
                                cal(j,4) =  j;
                                caller = 1;
                            else
                                cal(j,1:4) = 0;
                            end

                            if winSize == 2
                                if cal(j,1) > 1*10^-5  || cal(j,2) < 1.05  || cal(j,3) < 1.05
                                    B(A(j,2),A(j,1)) = 0;
                                    A(j,3) = 999;
                                    caller = 0;
                                end
                            elseif winSize == 3
                                if cal(j,1) > 1*10^-10  || cal(j,2) < 1.1
                                    B(A(j,2),A(j,1)) = 0;
                                    A(j,3) = 999;
                                    caller = 0;
                                end
                            elseif winSize == 4
                                if cal(j,1) > 1*10^-15  || cal(j,2) < 1.15
                                    B(A(j,2),A(j,1)) = 0;
                                    A(j,3) = 999;
                                    caller = 0;
                                end
                            end

% =========================== show results ========================
                            if showResult == 1
                                fprintf('\nRound%d  %dth-> p: [[%d]] %d, ratio: %d, ratio(pos back): %d', winSize-1, j, caller, cal(j,1), cal(j,2), cal(j,3));
                                if j == size(A,1)
                                    fprintf('\n   < %dXradius applied >',winSize);
                                    fprintf('\n');
                                end
                            end
%------------------------------------------------------------------
                        end
                    end
                end
                
                
                A(A(:,3) == 999,:) = [];
                thrrna{i} = imdilate(B, se);
            end
        end
    end
    fprintf('\nCh%d: %d/%d images,... %d/%d z-planes.', g, f, numOfImg, i, sizel);
end
fprintf('\n');




%%%------------------ Quick visual -----------------------
if showResult == 1
    figure, imshow(t1{i}*15)
    figure, imshow(t1{i}*15)
    hold on
    plot(A(:,1), A(:,2), 'r+');
    numbn = 1:size(A,1);
    for j = 1:size(A,1)
        text(A(j,1), A(j,2), num2str(numbn(j)), 'color', 'c');
    end
end
%%%--------------------------------------------------------




%%% 3D reconstitution using bwconncomp.
% blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
%        | 6th: matching nuc #(row # in 'nuc')
for i=1:size(thrrna,1)
    if isempty(thrrna{i})
        thrrna{i} = zeros(iminfo(3),iminfo(2));
    end
end
        
tfrna = cat(3,thrrna{:});
temp = bwconncomp(tfrna, 26);
conrna = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
blrna = struct2cell(conrna)';

% calculates intensity of detected blobs from original images.
iint = zeros(1,length(blrna(:,1)));
for i=1:length(blrna(:,1))
    iint(i) = 0;
    rng = round( blrna{i,3}(3) - (blrna{i,3}(6)-1)/2:blrna{i,3}(3) + (blrna{i,3}(6)-1)/2  );
    rng = rng(rng > 0 & rng <= iminfo(4));
    for j=1:length(rng)
        wholeImg = t1{rng(j)};
        rangeY = round(blrna{i,3}(2)):round(blrna{i,3}(2))+blrna{i,3}(5)-1;
        rangeX = round(blrna{i,3}(1)): round(blrna{i,3}(1))+blrna{i,3}(4)-1; 
        cutImg = wholeImg( rangeY  , rangeX  );
%         plot( rangeX(1) , rangeY(1), 'r+');
        iint(i) = iint(i) + sum(sum(blrna{i,4}(:,:,j) .* double(cutImg)));
    end
end

blrna(:,5) = num2cell(iint');
% fprintf('\nsize: %d %d\n', size(blrna,1), size(blrna,2));    

if  isempty(blrna) == 0
    if isempty(blrna{1,1}) == 0
        % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %            | col 6: total intensity | x centroid | y centroid | z centroid | col 10: # ATS in the blob.
        blobrna = zeros(length(blrna(:,1)),1);
        blobrna(:,1) = 1:length(blrna(:,1));
        temp = cell2mat(blrna(:,3));
        blobrna(:,2:4) = temp(:,4:6);
        temp = cell2mat(blrna(:,2));
        blobrna(:,7:9) = temp(:,1:3);

        blobrna(:,6) = cell2mat(blrna(:,5));
        blobrna(:,5) = cell2mat(blrna(:,1));
  
        %%% remove false-positive blobs
        if ~isempty(blobrna)
            % get mean of obj. size (pixels) and intensity (normalized: ind. int. - mean)
            mint = mean(blobrna(:,6));
%             mpix = mean(blobrna(:,5));
        end
        
        % Exclude RNA spots seen less than XX z-planes
        XX = 2;
        if ~isempty(blobrna)
            temp = blobrna(:,4) < XX;
            blobrna(temp,:) = [];
        end

        % Exclude RNA spots that are dim.
        if ~isempty(blobrna)
            temp = blobrna(:,4) < XX & blobrna(:,6) < mint * 0.15 ;
            blobrna(temp,:) = [];
        end

        if ~isempty(blobrna)
            %%% Exclude spots smaller than 2x2 pixels
            temp = blobrna(:,2) < pixr*2 & blobrna(:,3) < pixr*2;
            blobrna(temp ,:) = [];
        end
    else
        blobrna = [];
    end
else
    blobrna = [];
end



%%%------------- Quick visual-----------
% figure, imshow(rZproj*20)
% hold on
% plot(blobrna(:,7), blobrna(:,8), 'ro');
% for i = 1:size(blobrna,1)
%     text(blobrna(i,7), blobrna(i,8), num2str(i), 'color', 'c');
% %     plot(blobrna(:,7), blobrna(:,8), 'r+');
% end








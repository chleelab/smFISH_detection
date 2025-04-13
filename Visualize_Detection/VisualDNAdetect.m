function [] = VisualDNAdetect(af)

%%% ================== DNA (nuclei) display ===========================
for f=1:size(af,1)
%     figure, imshow(af{f,10}*5);
    blankp = ones(size(af{f,10},1),size(af{f,10},2));
%     blankp = ones(size(af{f,5},1),size(af{f,5},2));
    figure, imshow(blankp); hold on
    for i = length(af{f,1}(:,1)):-1:1
        minz = min(af{f,1}(:,3));
        maxz = max(af{f,1}(:,3));
        midz = mean(af{f,1}(:,3));
        diffz = maxz - minz;
        Bval = 1- (af{f,1}(i,3) - minz)/diffz;
        viscircles([af{f,1}(i,1) af{f,1}(i,2)], af{f,1}(i,4)-7, 'EdgeColor', [ 0 0 0 ], 'linewidth', 1);
        viscircles([af{f,1}(i,1) af{f,1}(i,2)], af{f,1}(i,4)-2, 'EdgeColor', [ Bval*1  Bval*1  1-Bval], 'linewidth', 5);
        viscircles([af{f,1}(i,1) af{f,1}(i,2)], af{f,1}(i,4)+3, 'EdgeColor', [ 0 0 0], 'linewidth', 1);
%         viscircles([af{f,1}(i,1) af{f,1}(i,2)], af{f,1}(i,4), 'EdgeColor', [ 0 0 0], 'linewidth', 5);
    end

%     labes = 1:size(af{f,1},1);
%     labes = arrayfun(@num2str, labes, 'uni', 0);
%     text(af{f,1}(:,1), af{f,1}(:,2),labes, 'color', 'r', 'fontsize', 20)
    fprintf('\n%d-th image', f); pause; close all
end
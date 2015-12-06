function [outputImage, colorfulLabels] = VideoGraphCutTextures


    function output = alphaMask(background, overlayedImage, offset, alphaMask)
        xMin = offset(1);
        yMin = offset(2);
        zMin = offset(3);
        xMax = offset(1)+size(overlayedImage, 1)-1;
        yMax = offset(2)+size(overlayedImage, 2)-1;
        if size(size(overlayedImage),2) == 4
            zMax = offset(3)+size(overlayedImage, 4)-1;
        else
            zMax = offset(3)+size(overlayedImage, 3)-1;
        end
        realMask = zeros(size(alphaMask, 1), size(alphaMask, 2), 3, size(alphaMask, 3));
        
        realMask(:,:,1,:) = alphaMask;
        realMask(:,:,2,:) = alphaMask;
        realMask(:,:,3,:) = alphaMask;
        %realMask = repmat(alphaMask, [1, 1, size(background, 3), 1]);
        %alphaMask
        
        size(overlayedImage)
        size(alphaMask)
        size(realMask)
        %size(background(xMin:xMax, yMin:yMax, :, zMin:zMax))
        if size(size(overlayedImage),2) == 4
            background(xMin:xMax, yMin:yMax, :, zMin:zMax) ...
                    = background(xMin:xMax, yMin:yMax, :, zMin:zMax).*(1-realMask) ...
                     + overlayedImage.*realMask;
        else
            background(xMin:xMax, yMin:yMax, zMin:zMax) ...
                    = background(xMin:xMax, yMin:yMax, zMin:zMax).*(1-alphaMask) ...
                     + overlayedImage.*alphaMask;
        end
        output = background;
    end
    function output = getOverlap(outputImagePatchLabels, overlayedImage, offset)
        xMin = offset(1);
        yMin = offset(2);
        zMin = offset(3);
        xMax = offset(1)+size(overlayedImage, 1)-1;
        yMax = offset(2)+size(overlayedImage, 2)-1;
        zMax = offset(3)+size(overlayedImage, 3)-1;
        output = zeros(size(outputImagePatchLabels));
        output(xMin:xMax, yMin:yMax, zMin:zMax) = outputImagePatchLabels(xMin:xMax, yMin:yMax, zMin:zMax)>0;
    end
    function mask = findSeam(patchOnBackground, outputImage)
        % compute derivatives (used for the graphcut cost)
        xDerivativeFilter = [-1/2, 0, 1/2];
        yDerivativeFilter = xDerivativeFilter';
        zDerivativeFilter = zeros(1, 1, 1, 3);
        zDerivativeFilter(1, 1, 1, 1) = -1/2;
        zDerivativeFilter(1, 1, 1, 3) = 1/2;

        patchOnBackgroundXDerivative = imfilter(patchOnBackground, xDerivativeFilter);
        patchOnBackgroundYDerivative = imfilter(patchOnBackground, yDerivativeFilter);
        patchOnBackgroundZDerivative = imfilter(patchOnBackground, zDerivativeFilter);
        outputImageXDerivative = imfilter(outputImage, xDerivativeFilter);
        outputImageYDerivative = imfilter(outputImage, yDerivativeFilter);
        outputImageZDerivative = imfilter(outputImage, zDerivativeFilter);
        
        size(patchOnBackgroundMask)
        size(outputImagePatchLabels)
        overlap = patchOnBackgroundMask & (outputImagePatchLabels>0);%getOverlap(outputImagePatchLabels, patch, offset);
        %figure;
        %subplot(1, 5, 1);
        %imshow(overlap);

        % Find the edges of the overlap region. We'll set graphcut constraints on
        % these edges
        innerEdgeFilter = zeros(3, 3, 3);
        innerEdgeFilter(:,:,2) = [0 -1 0; -1 6 -1; 0 -1 0];
        innerEdgeFilter(2, 2, 1) = -1;
        innerEdgeFilter(2, 2, 3) = -1;
        outerEdgeFilter = -innerEdgeFilter;
        overlapOuterEdges = imfilter(overlap, outerEdgeFilter)>0;
        constraintA = (outputImagePatchLabels>0) & overlapOuterEdges; % The pixels that are required to pertain to image A (the old patch)
        constraintB = patchOnBackgroundMask & overlapOuterEdges; % The pixels that are required to pertain to image B (the new patch)

        pixelsToInputToCutAlgorithmMask = constraintA | constraintB | overlap;

        %outputImageOuterEdges = imfilter(outputImagePatchLabels>0, outerEdgeFilter)>0;
        %figure;
        
        
        
        %subplot(1, 5, 2);
        %imshow(overlapOuterEdges);
        
        %subplot(1, 5, 3);
        %imshow(pixelsToInputToCutAlgorithmMask);

        %subplot(1, 5, 4);
        %imshow(constraintA);

        %subplot(1, 5, 5);
        %imshow(constraintB);

        graphCutInputCount = 0;
        graphCutIndexTable = zeros(size(outputImage)); % This table will map a position in the output image to the corresponding index in the graphcut array
        for i = 1:size(outputImage, 1)
            for j = 1:size(outputImage, 2)
                for k = 1:size(outputImage, 4)
                    if pixelsToInputToCutAlgorithmMask(i,j,k)
                        graphCutInputCount = graphCutInputCount + 1;
                        graphCutIndexTable(i, j, k) = graphCutInputCount;
                    end
                end
            end
        end
        graphCutInputCount2 = 0;
        nodeCosts = zeros(2, graphCutInputCount);
        for i = 1:size(outputImage, 1)
            for j = 1:size(outputImage, 2)
                for k = 1:size(outputImage, 4)
                    if pixelsToInputToCutAlgorithmMask(i,j,k)
                        graphCutInputCount2 = graphCutInputCount2 + 1;
                        newColumn = [0; 0];
                        if (constraintA(i,j,k)==1)
                            newColumn = [0;inf];
                        elseif (constraintB(i,j,k)==1)
                            newColumn = [inf;0];
                        end
                        nodeCosts(:, graphCutInputCount2) = newColumn;
                    end
                end
            end
        end
        
        loopIndex = 0;
        
        %disp('tik')
        %numGraphCutEdges = sum(sum(rightNeighboring)) + sum(sum(bottomNeighboring));
        
        rightNeighboring = pixelsToInputToCutAlgorithmMask & imfilter(pixelsToInputToCutAlgorithmMask, [0; 0; 1]);
        bottomNeighboring = pixelsToInputToCutAlgorithmMask & imfilter(pixelsToInputToCutAlgorithmMask, [0 0 1]);
        behindNeighboring = pixelsToInputToCutAlgorithmMask & imfilter(pixelsToInputToCutAlgorithmMask, reshape([0 0 1], [1 1 3])); % the reshape make sure that the [0 0 1] goes in the 3d dimension 
        
        distanceBetweenOutputImageAndPatchOnBackground = sqrt(sum((outputImage-patchOnBackground).^2, 3));
        
        vertFilterOfDistanceBetweenOutputImageAndPatchOnBackground = imfilter(distanceBetweenOutputImageAndPatchOnBackground, [0;1;1]);
        horFilterOfDistanceBetweenOutputImageAndPatchOnBackground = imfilter(distanceBetweenOutputImageAndPatchOnBackground, [0 1 1]);
        zDirectionFilterOfDistanceBetweenOutputImageAndPatchOnBackground = imfilter(distanceBetweenOutputImageAndPatchOnBackground, reshape([0 1 1], [1 1 3]));
        
        normOutputImageYDerivative = sqrt(sum(outputImageYDerivative.^2, 3));
        normPatchOnBackgroundYDerivative = sqrt(sum(patchOnBackgroundYDerivative.^2, 3));
        verticalDenominator = (1/10000)+imfilter(normOutputImageYDerivative + normPatchOnBackgroundYDerivative, [0; 1; 1]); % The (1/10000) is simply to avoid dividing by zero (not mentioned in the paper)
        verticalResult = vertFilterOfDistanceBetweenOutputImageAndPatchOnBackground ./ verticalDenominator;

        normOutputImageXDerivative = sqrt(sum(outputImageXDerivative.^2, 3));
        normPatchOnBackgroundXDerivative = sqrt(sum(patchOnBackgroundXDerivative.^2, 3));
        horizontalDenominator = (1/10000)+imfilter(normOutputImageXDerivative + normPatchOnBackgroundXDerivative, [0 1 1]);
        horizontalResult = horFilterOfDistanceBetweenOutputImageAndPatchOnBackground ./ horizontalDenominator;
        
        normOutputImageZDerivative = sqrt(sum(outputImageZDerivative.^2, 3));
        normPatchOnBackgroundZDerivative = sqrt(sum(patchOnBackgroundZDerivative.^2, 3));
        zDirectionDenominator = (1/10000)+imfilter(normOutputImageZDerivative + normPatchOnBackgroundZDerivative, reshape([0 1 1], [1 1 3]));
        zDirectionResult = zDirectionFilterOfDistanceBetweenOutputImageAndPatchOnBackground ./ zDirectionDenominator;

        indices = find(rightNeighboring)';

        graphCutEdgeIndicesI = graphCutIndexTable(indices);
        graphCutEdgeIndicesJ = graphCutIndexTable(indices+1);
        graphCutEdgeValues = verticalResult(indices);
        
        indices = find(bottomNeighboring)';
        % The following indices represent the sparsity pattern of our edge matrix.
        % The values will be stored inside 'graphCutEdgeValues'
        graphCutEdgeIndicesI = [graphCutEdgeIndicesI graphCutIndexTable(indices)];
        graphCutEdgeIndicesJ = [graphCutEdgeIndicesJ graphCutIndexTable(indices+size(pixelsToInputToCutAlgorithmMask, 1))];
        graphCutEdgeValues = [graphCutEdgeValues horizontalResult(indices)];
        
        
        indices = find(behindNeighboring)';
        % The following indices represent the sparsity pattern of our edge matrix.
        % The values will be stored inside 'graphCutEdgeValues'
        graphCutEdgeIndicesI = [graphCutEdgeIndicesI graphCutIndexTable(indices)];
        graphCutEdgeIndicesJ = [graphCutEdgeIndicesJ graphCutIndexTable(indices+size(pixelsToInputToCutAlgorithmMask, 1)*size(pixelsToInputToCutAlgorithmMask, 2))];
        graphCutEdgeValues = [graphCutEdgeValues zDirectionResult(indices)];
        % for IND = indices
        %     [i, j] = ind2sub(size(pixelsToInputToCutAlgorithmMask),IND);
        %     loopIndex = loopIndex + 1;
        %     graphCutEdgeValues(loopIndex) = ...
        %         (sqrt(sum((outputImage(i, j, :)-patchOnBackground(i, j, :)).^2)) ...
        %         + sqrt(sum((outputImage(i, j+1, :)-patchOnBackground(i, j+1, :)).^2))) ...
        %          / (sqrt(sum(outputImageXDerivative(i, j, :).^2)) ...
        %            + sqrt(sum(outputImageXDerivative(i, j+1, :).^2)) ...
        %            + sqrt(sum(patchOnBackgroundXDerivative(i, j, :).^2)) ...
        %            + sqrt(sum(patchOnBackgroundXDerivative(i, j+1, :).^2)));
        % end
        

        %disp('tok')
        graphCutEdges = sparse(graphCutEdgeIndicesI, graphCutEdgeIndicesJ, graphCutEdgeValues, graphCutInputCount, graphCutInputCount);
        %disp('tok')

        % Run Graphcut algorithm
        %BK_BuildLib;
        BK_LoadLib;

        size(graphCutEdges)
        size(nodeCosts)

        disp('Computing graph cut...');
        hinc = BK_Create(graphCutInputCount,2*graphCutInputCount);
        BK_SetNeighbors(hinc,graphCutEdges);
        BK_SetUnary(hinc,nodeCosts);
        e_inc = BK_Minimize(hinc);
        graphCutLabels = BK_GetLabeling(hinc)-1;
        disp('Computing graph cut...done.');
        
        

        mask = ones(size(patch, 1), size(patch, 2), size(patch, 4));
        for i = 1:size(patch, 1)
            for j = 1:size(patch, 2)
                for k = 1:size(patch, 4)
                    if pixelsToInputToCutAlgorithmMask(i+offset(1)-1,j+offset(2)-1,k+offset(3)-1)
                        graphCutIndex = graphCutIndexTable(i+offset(1)-1,j+offset(2)-1,k+offset(3)-1);
                        mask(i, j, k) = graphCutLabels(graphCutIndex);
                    end
                end
            end
        end
        
    end


%
% -------------------
%  Reading an Image!
% -------------------
% GIF
%[X, map] = imread('in/basket_flipped.gif');
%patch = ind2rgb(X, map);
%[X, map] = imread('in/jelly.gif');
%patch = ind2rgb(X, map);
%[X, map] = imread('in/nuts6.gif');
%patch = ind2rgb(X, map);

% PNG

% patch1 = im2double(imread('in/basket_png.png'));
% patch = zeros(size(patch1, 1), size(patch1, 2), size(patch1, 3), 2);
% patch(:,:,:,1) = patch1;
% patch(:,:,:,2) = patch1;
%patch2 = im2double(imread('in/basket_png_flipped.png'));

% -------------------
%  Reading a Video!
% -------------------
patchFile = VideoReader('in/flame_original.mpg'); % 240x320. 89 Frames
vidWidth = patchFile.Width;
vidHeight = patchFile.Height;
disp('Reading video!');
patch = zeros(vidHeight,vidWidth,3,0);
k = 1;
while hasFrame(patchFile)
    patch(:,:,:,k) = double(readFrame(patchFile))/255;
    k = k+1;
end
disp('Done reading video!');
class(patch)

%implay(immovie(patch), patchFile.FrameRate);

implay(immovie(patch));

size(patch)
%imshow(patch);

close all;
PLOT_X_NUM = 4;
PLOT_Y_NUM = 3;
plotIndex = 1;

OUT_IMG_WIDTH = 240;
OUT_IMG_HEIGHT = 320;
OUT_IMG_TIME_LENGTH = 150;

% TODO: remove the following line
% Extrudes the 2d image into 3d
%patch = repmat(patch, [1, 1, 1, 20]);
%patch2 = repmat(patch2, [1, 1, 1, 20]);

outputImage = zeros(OUT_IMG_WIDTH, OUT_IMG_HEIGHT, 3, OUT_IMG_TIME_LENGTH);

outputImagePatchLabels = zeros(OUT_IMG_WIDTH, OUT_IMG_HEIGHT, OUT_IMG_TIME_LENGTH);

%% Place the first patch
outputImage(1:size(patch, 1), 1:size(patch, 2), :, 1:size(patch,4)) = patch;
outputImagePatchLabels(1:size(patch, 1), 1:size(patch, 2), 1:size(patch,4)) = ones(size(patch, 1), size(patch, 2), size(patch,4));

%subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
%plotIndex = plotIndex + 1;



implay(immovie(outputImage))

%figure
%subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
%plotIndex = plotIndex + 1;

%imagesc(outputImagePatchLabels)


%patchBackup = patch;

%patch = patch2;
label = 2;
%for label=2:100
    %% Find the offset for the next patch. By offset, we just mean the coordinates of the new patch.
    % known bug: when the offset is [2 2], graphcut algorithm goes into
    % infinite loop
    %offset = [10, 100];
    
    figure;
    implay(immovie(outputImage));
    %offset = uint16(ginput(1));
    %offset = [offset(2), offset(1), 1];
    
    offset = [1, 1, 46]
    
    close all hidden;
    
    

    %% plot the offset second patch
    patchOnBackground = zeros(size(outputImage));
    patchOnBackgroundMask = zeros(size(outputImage, 1), ...
                                  size(outputImage, 2), ...
                                  size(outputImage, 4));
    size(patchOnBackground)
    offset(1)+size(patch, 1)-1
    offset(2)+size(patch, 2)-1
    offset(3)+size(patch, 4)-1
    patchOnBackground(offset(1):offset(1)+size(patch, 1)-1, ...
                      offset(2):offset(2)+size(patch, 2)-1, ...
                      :, ...
                      offset(3):offset(3)+size(patch, 4)-1) = patch;
    patchOnBackgroundMask(offset(1):offset(1)+size(patch, 1)-1, ...
                          offset(2):offset(2)+size(patch, 2)-1, ...
                          offset(3):offset(3)+size(patch, 4)-1) = 1;
    %subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
    %plotIndex = plotIndex + 1;
    %imshow(patchOnBackground);
    implay(immovie(patchOnBackground));
    
    

    %% Find the seam (i.e. compute the mask) for the second patch
    % Find the overlapping region
    disp('Finding the seam...');

    %mask = ones(size(patch, 1), size(patch, 2));
    mask = findSeam(patchOnBackground, outputImage);

    disp('Finding the seam...done.');
    
    


    %% Place the new patch with the calculated offset and seam

    %mask = round(rand(size(patch, 1), size(patch, 2)));
    outputImage = alphaMask(outputImage, patch, offset, mask);
    
    %outputImage = outputImage(:, :, :, 45:89);
    
    implay(immovie(outputImage));
    
    writerObj = VideoWriter('out/movieout.avi');
    open(writerObj);
    writeVideo(writerObj, immovie(outputImage));
    
    
    outputImagePatchLabels = alphaMask(outputImagePatchLabels, mask*label, offset, mask);
    colorfulLabels = reshape(outputImagePatchLabels, [size(outputImagePatchLabels, 1), size(outputImagePatchLabels, 2), 1, size(outputImagePatchLabels, 3)]);
    colorfulLabels = repmat(colorfulLabels, [1, 1, 3, 1]);
    colorfulLabels = colorfulLabels * 0.20;
    
    %colorfulLabels = colorfulLabels(:, :, :, 45:89);
    
    implay(immovie(colorfulLabels));
    
    writerObj2 = VideoWriter('out/seamout.avi');
    open(writerObj2);
    writeVideo(writerObj2, immovie(colorfulLabels));
    
    
    
    %{
    
    %maskLogical = mask == ones(size(patch, 1), size(patch, 2));
    %outputImage((1+offset(1)):(size(patch, 1)+offset(1)), (1+offset(2)):(size(patch, 2)+offset(2)), :) = patch .* repmat(mask, 1, 1, 3);

    %outputImagePatchLabels((1+offset(1)):(size(patch, 1)+offset(1)), (1+offset(2)):(size(patch, 2)+offset(2))) = 2*mask;

    figure
    imshow(outputImage)
    figure
    imagesc(outputImagePatchLabels)
end
%}
end


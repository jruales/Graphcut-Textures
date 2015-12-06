function GraphCutTextures

    %% Helper functions. The actual program begins after these functions

    % Alpha blends the overlayedImage on top of the background image using
    % the given alphaMask. a 1 on the alpha mask represents full opacity of
    % the overlayedImage, and a 0 represents full transparency
    function output = alphaMask(background, overlayedImage, offset, alphaMask)
        xMin = offset(1);
        yMin = offset(2);
        xMax = offset(1)+size(overlayedImage, 1)-1;
        yMax = offset(2)+size(overlayedImage, 2)-1;
        realMask = repmat(alphaMask, 1, 1, size(background, 3));
        %alphaMask
        background(xMin:xMax, yMin:yMax, :) ...
                = background(xMin:xMax, yMin:yMax, :).*(1-realMask) ...
                 + overlayedImage.*realMask;
        output = background;
    end

    %{
    function output = getOverlap(outputImagePatchLabels, overlayedImage, offset)
        xMin = offset(1);
        yMin = offset(2);
        xMax = offset(1)+size(overlayedImage, 1)-1;
        yMax = offset(2)+size(overlayedImage, 2)-1;
        output = zeros(size(outputImagePatchLabels));
        output(xMin:xMax, yMin:yMax) = outputImagePatchLabels(xMin:xMax, yMin:yMax)>0;
    end
    %}

    % Uses graphcuts to find the seam to apply to a new patch given the
    % offset. The input "patchOnBackground" should be an image the size of
    % the output image with the patch image as a subimage at the desired
    % offset
    function mask = findSeam(patchOnBackground, patchOnBackgroundMask, outputImage)
        % compute derivatives (used for the graphcut cost)
        xDerivativeFilter = [-1/2, 0, 1/2];
        yDerivativeFilter = xDerivativeFilter';

        patchOnBackgroundXDerivative = imfilter(patchOnBackground, xDerivativeFilter);
        patchOnBackgroundYDerivative = imfilter(patchOnBackground, yDerivativeFilter);
        outputImageXDerivative = imfilter(outputImage, xDerivativeFilter);
        outputImageYDerivative = imfilter(outputImage, yDerivativeFilter);
        
        overlap = patchOnBackgroundMask & (outputImagePatchLabels>0);%getOverlap(outputImagePatchLabels, patch, offset);
        figure;
        subplot(1, 5, 1);
        imshow(overlap);

        % Find the edges of the overlap region. We'll set graphcut constraints on
        % these edges
        innerEdgeFilter = [0 -1 0; -1 4 -1; 0 -1 0];
        outerEdgeFilter = -innerEdgeFilter;
        overlapOuterEdges = imfilter(overlap, outerEdgeFilter)>0;
        constraintA = (outputImagePatchLabels>0) & overlapOuterEdges; % The pixels that are required to pertain to image A (the old patch)
        constraintB = patchOnBackgroundMask & overlapOuterEdges; % The pixels that are required to pertain to image B (the new patch)

        pixelsToInputToCutAlgorithmMask = constraintA | constraintB | overlap;

        %outputImageOuterEdges = imfilter(outputImagePatchLabels>0, outerEdgeFilter)>0;
        %figure;
        subplot(1, 5, 2);
        imshow(overlapOuterEdges);
        %figure;
        subplot(1, 5, 3);
        imshow(pixelsToInputToCutAlgorithmMask);

        subplot(1, 5, 4);
        imshow(constraintA);

        subplot(1, 5, 5);
        imshow(constraintB);

        % Set up table that maps pixels in the output image to node indices
        % in the graphcut algorithm
        graphCutInputCount = 0;
        graphCutIndexTable = zeros(size(outputImage)); % This table will map a position in the output image to the corresponding index in the graphcut array
        for i = 1:size(outputImage, 1)
            for j = 1:size(outputImage, 2)
                if pixelsToInputToCutAlgorithmMask(i,j)
                    graphCutInputCount = graphCutInputCount + 1;
                    graphCutIndexTable(i, j) = graphCutInputCount;
                end
            end
        end
        % Set graphcut infinity constraints at the non-overlapping portions
        % of the images for the graphcut algorithm
        graphCutInputCount2 = 0;
        nodeCosts = zeros(2, graphCutInputCount);
        for i = 1:size(outputImage, 1)
            for j = 1:size(outputImage, 2)
                if pixelsToInputToCutAlgorithmMask(i,j)
                    graphCutInputCount2 = graphCutInputCount2 + 1;
                    newColumn = [0; 0];
                    if (constraintA(i,j)==1)
                        newColumn = [0;inf];
                    elseif (constraintB(i,j)==1)
                        newColumn = [inf;0];
                    end
                    nodeCosts(:, graphCutInputCount2) = newColumn;
                end
            end
        end
        
        loopIndex = 0;
        %disp('tik')
        %numGraphCutEdges = sum(sum(rightNeighboring)) + sum(sum(bottomNeighboring));
        
        % boolean arrays denoting whether or not a pixel to be used in the
        % graphcut graph has a right or bottom neighboring pixel that will also be part of
        % the graphcut graph
        rightNeighboring = pixelsToInputToCutAlgorithmMask & imfilter(pixelsToInputToCutAlgorithmMask, [0; 0; 1]);
        bottomNeighboring = pixelsToInputToCutAlgorithmMask & imfilter(pixelsToInputToCutAlgorithmMask, [0 0 1]);
        
        % This distance will be important when calculating the edge costs
        distanceBetweenOutputImageAndPatchOnBackground = sqrt(sum((outputImage-patchOnBackground).^2, 3));
        
        % This sum of adjacent distances will be important when calculating
        % the edge costs
        vertFilterOfDistanceBetweenOutputImageAndPatchOnBackground = imfilter(distanceBetweenOutputImageAndPatchOnBackground, [0;1;1]);
        horFilterOfDistanceBetweenOutputImageAndPatchOnBackground = imfilter(distanceBetweenOutputImageAndPatchOnBackground, [0 1 1]);
        
        % These derivatives will be important when calculating
        % the edge costs
        normOutputImageYDerivative = sqrt(sum(outputImageYDerivative.^2, 3));
        normPatchOnBackgroundYDerivative = sqrt(sum(patchOnBackgroundYDerivative.^2, 3));
        normOutputImageXDerivative = sqrt(sum(outputImageXDerivative.^2, 3));
        normPatchOnBackgroundXDerivative = sqrt(sum(patchOnBackgroundXDerivative.^2, 3));
        
        % These denominators use the derivatives above and will be important when calculating
        % the edge costs. In order to avoid divide by zero issues, we add a
        % small value to these denominators, although the paper doesn't
        % address this issue.
        verticalDenominator = (1/10000)+imfilter(normOutputImageYDerivative + normPatchOnBackgroundYDerivative, [0; 1; 1]);
        verticalResult = vertFilterOfDistanceBetweenOutputImageAndPatchOnBackground ./ verticalDenominator;
        horizontalDenominator = (1/10000)+imfilter(normOutputImageXDerivative + normPatchOnBackgroundXDerivative, [0 1 1]);
        horizontalResult = horFilterOfDistanceBetweenOutputImageAndPatchOnBackground ./ horizontalDenominator;

        % We now construct the sparse edge matrix
        
        % Process horizontal edges
        indices = find(rightNeighboring)';
        graphCutEdgeIndicesI = graphCutIndexTable(indices);
        graphCutEdgeIndicesJ = graphCutIndexTable(indices+1);
        graphCutEdgeValues = verticalResult(indices);
        
        % Add vertical edges
        indices = find(bottomNeighboring)';
        graphCutEdgeIndicesI = [graphCutEdgeIndicesI graphCutIndexTable(indices)];
        graphCutEdgeIndicesJ = [graphCutEdgeIndicesJ graphCutIndexTable(indices+size(pixelsToInputToCutAlgorithmMask, 1))];
        graphCutEdgeValues = [graphCutEdgeValues horizontalResult(indices)];
        
        % This is a legacy 'for' loop that I had to replace with the
        % vectorized version for speed
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
        % construct the sparse edges from the weights and connections
        % specified above
        graphCutEdges = sparse(graphCutEdgeIndicesI, graphCutEdgeIndicesJ, graphCutEdgeValues, graphCutInputCount, graphCutInputCount);
        %disp('tok')

        % Run Graphcut algorithm
        %BK_BuildLib;
        % Load Graphcut library
        BK_LoadLib;

        size(graphCutEdges)
        size(nodeCosts)

        % Compute the graph cut using the graph cut library
        disp('Computing graph cut...');
        hinc = BK_Create(graphCutInputCount,2*graphCutInputCount);
        BK_SetNeighbors(hinc,graphCutEdges);
        BK_SetUnary(hinc,nodeCosts);
        e_inc = BK_Minimize(hinc);
        graphCutLabels = BK_GetLabeling(hinc)-1;
        disp('Computing graph cut...done.');

        % Based on the computed graph cut, create and output an appropriate
        % mask (this will represent the seam)
        mask = ones(size(patch, 1), size(patch, 2));
        for i = 1:size(patch, 1)
            for j = 1:size(patch, 2)
                if pixelsToInputToCutAlgorithmMask(i+offset(1)-1,j+offset(2)-1)
                    graphCutIndex = graphCutIndexTable(i+offset(1)-1,j+offset(2)-1);
                    mask(i, j) = graphCutLabels(graphCutIndex);
                end
            end
        end
    end

%% Import the input texture...
% ...in GIF format...
%[X, map] = imread('in/basket_flipped.gif');
%patch = ind2rgb(X, map);
%[X, map] = imread('in/jelly.gif');
%patch = ind2rgb(X, map);
%[X, map] = imread('in/nuts6.gif');
%patch = ind2rgb(X, map);

% ...or in PNG format
patch = im2double(imread('in/basket_png.png'));
%patch = im2double(imread('in/basket_png_flipped.png'));

size(patch)
%imshow(patch);

close all;


PLOT_X_NUM = 4;
PLOT_Y_NUM = 3;
plotIndex = 1;

% These are the dimensions of the output image. They should be at least as
% big as the dimensions of the input texture.
OUT_IMG_WIDTH = 800;
OUT_IMG_HEIGHT = 800;

% This will be the actual output image (in double RGB format)
outputImage = zeros(OUT_IMG_HEIGHT, OUT_IMG_WIDTH, 3);

% This will be an image that indicates what parts of the output image
% belong to what copy of the input image. A value of zero indicates that
% a pixel is empty, and no patch has been placed on it yet
outputImagePatchLabels = zeros(OUT_IMG_HEIGHT, OUT_IMG_WIDTH);

%% Place the first patch
% simply copy the entire patch onto the output image
outputImage(1:size(patch, 1), 1:size(patch, 2), :) = patch;
% update the label image so that the updated values have a label of "1"
outputImagePatchLabels(1:size(patch, 1), 1:size(patch, 2)) = ones(size(patch, 1), size(patch, 2));

subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
plotIndex = plotIndex + 1;

imshow(outputImage)
%figure
subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
plotIndex = plotIndex + 1;

imagesc(outputImagePatchLabels)


patchBackup = patch;

% Repeatedly let the user place new patches. After
% being done, the user can save the output image through the MATLAB imshow interface and quit the program using CTRL-C.
for label=2:1000
    %% Find the offset for the next patch. By offset, we just mean the coordinates of the new patch, so an offset of [1, 1] represents placing a patch on the top left corner.
    % known bug: when the offset is [2 2], graphcut algorithm goes into
    % infinite loop
    %offset = [50, 100];
    
    figure;
    imshow(outputImage);
    % Let the user select the location where to put the next image patch
    offset = int32(ginput(1))
    % the above gives [x y], but we want [y x], so we swap the entries
    offset = [offset(2) offset(1)];
    
    
    patch = patchBackup;
    
    % If the offset causes the image to get off-screen, clip the patch so
    % that is fits entirely on-screen
    patch = patch(max(1,offset(1))-offset(1)+1 : min(size(outputImage,1),offset(1)+size(patch, 1)-1)-offset(1)+1, ...
                  max(1,offset(2))-offset(2)+1 : min(size(outputImage,2),offset(2)+size(patch, 2)-1)-offset(2)+1, ...
                  :);
    
    % adjust the offset for the clipped patch, if necessary
    offset = [max(1, offset(1)), max(1, offset(2))];
    
    %offset = [100 1]
    
    close all;
    
    plotIndex = 1;

    %% plot the offset second patch
    patchOnBackground = zeros(size(outputImage));
    patchOnBackgroundMask = zeros(size(outputImage, 1), size(outputImage, 2));
    patchOnBackground(offset(1):offset(1)+size(patch, 1)-1, ...
                      offset(2):offset(2)+size(patch, 2)-1, ...
                      :) ...
                      = patch;
    patchOnBackgroundMask(offset(1):offset(1)+size(patch, 1)-1, ...
                          offset(2):offset(2)+size(patch, 2)-1, ...
                          :) = 1;
    subplot(PLOT_X_NUM, PLOT_X_NUM, plotIndex);
    plotIndex = plotIndex + 1;
    imshow(patchOnBackground);

    %% Find the seam (i.e. compute the mask) for the second patch
    % Find the overlapping region
    disp('Finding the seam...');

    %mask = ones(size(patch, 1), size(patch, 2));
    % The output mask is a binary 2d image the size of the output image
    % that uses a 1 to inticate that a pixel should come from the new,
    % offset texture, and a 0 to indicate that a pixel should come from the
    % old current output texture
    mask = findSeam(patchOnBackground, patchOnBackgroundMask, outputImage);

    disp('Finding the seam...done.');


    %% Place the new patch with the calculated offset and seam

    %mask = round(rand(size(patch, 1), size(patch, 2)));
    
    % The alphaMask function uses the mask computed above as a transparency
    % mask for the new patch, combining the new patch with the current background
    % image
    outputImage = alphaMask(outputImage, patch, offset, mask);
    outputImagePatchLabels = alphaMask(outputImagePatchLabels, mask*label, offset, mask);

    %maskLogical = mask == ones(size(patch, 1), size(patch, 2));
    %outputImage((1+offset(1)):(size(patch, 1)+offset(1)), (1+offset(2)):(size(patch, 2)+offset(2)), :) = patch .* repmat(mask, 1, 1, 3);

    %outputImagePatchLabels((1+offset(1)):(size(patch, 1)+offset(1)), (1+offset(2)):(size(patch, 2)+offset(2))) = 2*mask;

    % Show the output image and the seam so that we can save them, or add a
    % new patch.
    figure
    imshow(outputImage)
    figure
    imagesc(outputImagePatchLabels)
end
end


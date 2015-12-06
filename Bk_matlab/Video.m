function Video()

patchFile = VideoReader('in/fountain_original.mpg');

vidWidth = patchFile.Width
vidHeight = patchFile.Height

disp('Reading video!');

patch = zeros(vidHeight,vidWidth,3,0,'uint8');
k = 1;
while hasFrame(patchFile)
    patch(:,:,:,k) = readFrame(patchFile);
    k = k+1;
end

disp('Done reading video!');

implay(immovie(patch), patchFile.FrameRate);

writerObj = VideoWriter('out/movieout.avi');
open(writerObj);
writeVideo(writerObj, immovie(patch));

end
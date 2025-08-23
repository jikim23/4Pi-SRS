function stack = readTiffStack(filename)
    info = imfinfo(filename);
    nFrames = numel(info);
    stack = zeros(info(1).Height, info(1).Width, nFrames, 'double');
    for k = 1:nFrames
        stack(:,:,k) = im2double(imread(filename, k));
    end
end

function writeStackToTiff(vol, fname)
% WRITESTACKTOTIFF  Save a 3D array as a multi-page 16-bit TIFF
%   vol   – 3D array [Y × X × Z], can be double, single, uint8, etc.
%   fname – output filename, e.g. 'myVolume.tif'

    %--- 1) Convert to uint16 if needed (ImageJ reads 8/16-bit natively) ---
    if ~isa(vol,'uint8') && ~isa(vol,'uint16')
        % scale double or single into [0,65535]
        vol = uint16( 65535 * mat2gray(vol) );
    end

    %--- 2) Prepare the tag structure ---
    tagstruct.ImageLength          = size(vol,1);
    tagstruct.ImageWidth           = size(vol,2);
    tagstruct.Photometric          = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample        = 16;
    tagstruct.SamplesPerPixel      = 1;
    tagstruct.RowsPerStrip         = 16;              % or size(vol,1)
    tagstruct.PlanarConfiguration  = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression          = Tiff.Compression.None;
    tagstruct.SampleFormat         = Tiff.SampleFormat.UInt;
    tagstruct.Software             = 'MATLAB';

    %--- 3) Open the file and write each slice ---
    t = Tiff(fname,'w');
    for k = 1:size(vol,3)
        t.setTag(tagstruct);
        t.write(vol(:,:,k));
        if k < size(vol,3)
            t.writeDirectory();   % create next IFD
        end
    end
    t.close();
end

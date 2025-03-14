clear all;
%%
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';

%%
saveData = false;

D = gpuDevice;

%%
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/Calibration60X/mito/OMP/241211';
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/Aliasing_OPM/U2OS/241211';
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/LabMembers/Conor/omniOPM/Reto/Aliasing_OPM/U2OS/241211';
dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/U2OS/ER-Henne-dish2/241217';

cellStr = 'Cell1';

chanToken = 'CH00';

tiffList = dir(fullfile(dataPath, cellStr, '*.tif*'));

%% microscope params

% omniOPM
dsFactor = 3;
xyPixelSize = 0.147;
dz = 0.207;
skewAngle = 45.0;

%% reconstruction params

zCrop = 128;

%% set up save folder

condition = sprintf("recon_dz=%dnm", uint16((dz*dsFactor)*1000));
saveDir = fullfile(dataPath, cellStr, condition);

if ~exist(saveDir, "dir")
    mkdir(saveDir);
end

%% loop through timepoints

maxNumTimePoints = 0;

if maxNumTimePoints <= 0
    maxNumTimePoints = length(tiffList);
end

for t = 1:maxNumTimePoints

    timePoint = tiffList(t);

    fprintf("Processing t = %d : %s...\n", t, timePoint.name);

    display = (t == 1) || (t == maxNumTimePoints);

    if ~contains(timePoint.name, chanToken)
        continue;
    end

    imPath = fullfile(timePoint.folder, timePoint.name);

    %% load, deskew and downsample image
    
    im = permute(readtiff(imPath), [2,1,3]);
    
    % push to GPU
    im = gpuArray(single(im));

    % deskew
    im_dsk = deskewFrame3D( ...
        im, ...
        skewAngle, ...
        dz * dsFactor, ...
        xyPixelSize, ...
        'reverse', true ...
        );
    
    outSize = [size(im_dsk, [1,2]), zCrop];
    
    im_dsk(im_dsk(:) == 0) = median(im_dsk(im_dsk(:) > 0));

    im_dsk_rot = rotateFrame3D( ...
        im_dsk, ...
        skewAngle, ...
        dsFactor, ...
        'reverse', true, ...
        'Crop', true, ...
        'outSize', outSize ...
        );
    im_dsk_rot = norm_u16(im_dsk_rot);
    
    if saveData
        writetiff(im_dsk_rot, fullfile(saveDir, timePoint.name));
    end
    
    if display
        figure(t); clf;
        subplot(3,3,[2,5,8]);
        set(gcf, 'color', [1,1,1]);
        % imagesc(mip(im_dsk_rot, 1));
        imagesc([squeeze(im_dsk_rot(:,:,48)); squeeze(im_dsk_rot(400,:,:))']);
        colormap(gca, "parula");
        axis image;
    end

    %% fft of interp
    G_interp = fftn(im_dsk_rot);
    
    if display
        subplot(3,3,7);
        imagesc(pwr(G_interp, 1));
        colormap(gca, "hot");
        axis image;        
    end

    %% "upsampling" downsampled stack: G
    
    % comb upsampling trick
    im_comb = zeros(size(im_dsk) .* [1 1 dsFactor], 'like', im);
    im_comb(:, :, 1:dsFactor:end) = im_dsk;

    % fft
    G = fftn(im_comb);
    
    if display
        subplot(3,3,1);
        imagesc(pwr(G, 1));
        colormap(gca, "hot");
        axis image;
    end

    %% create mask and apply to G
    % fast mask
    mask = single(im_comb);
    mask(:) = 0;

    [~, ~, nz] = size(mask);
    hz = uint16(nz/2);
    wz = uint16((xyPixelSize/dz)*nz/dsFactor/2);
    
    mask(:, :, (hz-wz):(hz+wz)) = 1;
    mask = rotateFrame3D( ...
            mask, ...
            skewAngle * dz/xyPixelSize/dsFactor, ...
            dz/xyPixelSize, ...
            'Crop', true, ...
            'outSize', size(mask) ...
        );
    mask(mask < 1) = 0;
    mask = fftshift(mask);    

    % Gaussian blur, if you're in the mood
    blurSize = 0.0;
    if blurSize > 0
        mask = imgaussfilt3(double(mask), blurSize/dsFactor);
    end
    
    G_mask = mask .* G;
     
    if display
        subplot(3,3,4);
        imagesc(pwr(G_mask, 1));
        colormap(gca, "hot");
        axis image;
    end

    %% ifft on G
    g_recon = ifftn(G_mask);
    g_recon = real(g_recon);

    %% rotate recon
    im_recon_rot = rotateFrame3D( ...
        g_recon, ...
        skewAngle, ...
        1.0, ...
        'reverse', true, ...
        'Crop', true, ...
        'outSize', outSize ...
        );
    im_recon_rot = norm_u16(im_recon_rot);

    if display
        subplot(3,3,[3,6,9]);
        % imagesc(mip(im_recon_rot, 1));
        imagesc([squeeze(im_recon_rot(:,:,48)); squeeze(im_recon_rot(400,:,:))']);
        colormap(gca, "parula");
        axis image;
    end

    %% save
    if saveData
        writetiff(im_recon_rot, fullfile(saveDir, ['recon_', timePoint.name]));
    end
end

reset(D);

%% functions

function [ out ] = norm_u16( in )
    out = double(in);
    out = out - median(out(:));
    out(out < 0) = 0;
    out = out ./ max(out(:));
    out = uint16(65535 .* out);
end
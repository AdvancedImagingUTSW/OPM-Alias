clear all;

%% helper functions
mip = @(x,dim) squeeze(max(rescale(x),[],dim))';
pwr = @(x,dim) squeeze(log10(max(abs(fftshift(x)),[],dim)))';

%% Load TIF from folder. All TIF will be loaded if time series.

dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/Manuscripts/OPM-ALIAS/DataToShare/mesoOPM'; experimentName = 'mesoOPM';
% dataPath = '/archive/bioinformatics/Danuser_lab/Fiolka/MicroscopeDevelopment/omniOPM/U2OS/ER-Henne-dish2/241217/Cell1'; experimentName = 'omniOPM';

tiffList = dir(fullfile(dataPath, '*.tif*'));

pathList = {};
for i = 1:length(tiffList)
    pathList{end+1} = fullfile(tiffList(i).folder, tiffList(i).name);
end

%% parameter setup: use experimentName

% defaults
baseThick = 0;
zCrop = 0;

switch experimentName
    case 'omniOPM'
        osFactor = 1;
        dsFactor = 3;
        xyPixelSize = 0.147;
        dz = 0.207;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 5.0;
        outSize = [512 256 128];
    
    case 'highOPM2'
        osFactor = 2;
        dsFactor = 2;
        xyPixelSize = 0.147;
        dz = 0.207;
        skewAngle = 45.0;
        fillMethod = 'crop';
        angleBump = 3.0;
        outSize = [16,256,96];
        baseThick = 0.5;
        zCrop = 56;
    
    case 'mesoOPM'
        osFactor = 1;
        dsFactor = 4;
        xyPixelSize = 1.15;
        dz = 1.60;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 0.0;
        outSize = [256,128,64];

    otherwise
        osFactor = 1;
        dsFactor = 1;
        xyPixelSize = 0.147;
        dz = 0.210;
        skewAngle = 45.0;
        fillMethod = 'median';
        angleBump = 0.0;   
        outSize = [256,256,64];
end

dz = dz * osFactor;

%% colormap definition

cmap = hsv;
cmap = cmap(1:16:end, :);

colorMap = containers.Map;
colorMap('deskewFrame3D')       = cmap(1,:);
colorMap('rotateFrame3D')       = cmap(3,:);
colorMap('deskewRotateFrame3D') = cmap(2,:);
colorMap('imwarp(z)')           = cmap(4,:);
colorMap('combUpsample')        = cmap(9,:);
colorMap('fftn')                = cmap(6,:);
colorMap('ifftn')               = cmap(8,:);
colorMap('otfMask')             = cmap(10,:);
colorMap('gpuArray')            = cmap(12,:);
colorMap('gpuArray-1')          = cmap(12,:);
colorMap('gpuArray-2')          = cmap(12,:);
colorMap('createMask')          = cmap(15,:);

% plot colors
figure(10); clf;
subplot(1,1,1);
hold on;
dx = 1; dy = 1;
x = 0; y = 0;
for i = 1:length(cmap)
    p = [x y; x y+dy; x+dx y+dy; x+dx y];
    fill(p(:,1), p(:,2), cmap(i,:));
    text(x, y + dy/2, int2str(i));
    y = y + dy;
end
hold off;
axis off;

%% initialize figs

figure(1); clf;
set(gcf, 'color', [1,1,1]);
tiledlayout(1,1, 'Padding', 'loose', 'TileSpacing', 'none');

% GPU
D = gpuDevice;

% number of trials
nTimePts = length(pathList);
nRuns = 5;

if nTimePts > 1
    nRuns = 1;
end

figure(2); clf;
set(gcf, 'color', [1,1,1]);

%% TEST 1: Simply load data and deskewRotateFrame3D

test1_Map = containers.Map;
test1_Name = 'PetaKit5D: fast DSR';

for j = 1:nRuns
    for i = 1:nTimePts
        imPath = pathList{i};
        
        % load image
        im_raw = single(permute(readtiff(imPath), [2,1,3]));
    
        % deskew-rotate
        tic;
        im_dsr = deskewRotateFrame3D( ...
                im_raw, ...
                skewAngle, ...
                dz * dsFactor, ...
                xyPixelSize, ...
                'reverse', true, ...
                'Crop', true ...
        );
        test1_Map = addTime(test1_Map, 'deskewRotateFrame3D', toc);
    end
end

% plot
subplot(2,3,1);
im_crop = ctrcrop(im_dsr, outSize);
imagesc(mip(im_crop, 1));
axis image;
title(test1_Name);

% release GPU
reset(D);

%% TEST 2: Load -> Comb Upsample -> deskewRotateFrame3D -> FT -> Mask -> IFT

test2_Map = containers.Map;
test2_Name = 'Aliasing method: Top-Hat filter after DSR';

for j = 1:nRuns
    for i = 1:nTimePts
        imPath = pathList{i};
        
        % load image
        im_raw = single(permute(readtiff(imPath), [2,1,3]));
    
        % create gpuArrays
        tic;
        [sx, sy, sz] = size(im_raw);
        im_raw = gpuArray(im_raw);
        im_comb = zeros(sx, sy, sz*dsFactor, 'like', im_raw);
        test2_Map = addTime(test2_Map, 'gpuArray-1', toc);
    
        % comb upsampling
        tic;
        im_comb(:, :, 1:dsFactor:end) = im_raw;
        test2_Map = addTime(test2_Map, 'combUpsample', toc);
    
        % deskew-rotate
        tic;
        im_dsr = deskewRotateFrame3D( ...
                im_comb, ...
                skewAngle + angleBump, ...
                dz, ...
                xyPixelSize, ...
                'reverse', true, ...
                'Crop', true ...
            );
        test2_Map = addTime(test2_Map, 'deskewRotateFrame3D', toc);
    
        % fftn
        tic;
        G = fftn(im_dsr);
        test2_Map = addTime(test2_Map, 'fftn', toc);
    
        G_full = gather(G);
    
        tic;
        [~, ~, nz] = size(G);
        wz = uint16((xyPixelSize/dz)*nz/dsFactor/2);
        G(:, :, wz:(end-wz)) = 0;
        test2_Map = addTime(test2_Map, 'otfMask', toc);
    
        tic;
        im_recon = real(ifftn(G));
        test2_Map = addTime(test2_Map, 'ifftn', toc);
    end
end

% plot
subplot(2,3,5);
imagesc(pwr(G_full, 1));
hold on;
hz = uint16(nz/2);
yline(hz-wz, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
yline(hz+wz, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
hold off;
set(gca, 'colormap', hot);
title('Masked OTF');

subplot(2,3,2);
im_crop = ctrcrop(im_recon, outSize);
imagesc(mip(im_crop,1));
axis image;
set(gca, 'colormap', parula);
title(test2_Name);

% release GPU
reset(D);

%% TEST 4: Load -> Comb Upsample -> deskewFrame3D -> FT -> Mask -> IFT -> rotateFrame3D

test4_Map = containers.Map;
tiltedMask_Map = containers.Map;
test4_Name = 'Aliasing method: Top-Hat filter after deskew, before rotation';

for j = 1:nRuns
    for i = 1:nTimePts
        imPath = pathList{i};
        
        % load image
        im_raw = single(permute(readtiff(imPath), [2,1,3]));
    
        % create gpuArrays
        tic;
        [sx, sy, sz] = size(im_raw);
        im_raw = gpuArray(im_raw);
        im_comb = zeros(sx, sy, sz*dsFactor, 'like', im_raw);
        test4_Map = addTime(test4_Map, 'gpuArray-1', toc);
    
        % comb upsampling
        tic;
        im_comb(:, :, 1:dsFactor:end) = im_raw;
        test4_Map = addTime(test4_Map, 'combUpsample', toc);
    
        % deskew
        tic;
        im_dsk = deskewFrame3D( ...
                im_comb, ...
                skewAngle, ...
                dz, ...
                xyPixelSize, ...
                'reverse', true ...
            );    
        test4_Map = addTime(test4_Map, 'deskewFrame3D', toc);
        
        % mask creation
        tic;
    
        mask = single(im_dsk);
        mask(:) = 0;
    
        [~, ~, nz] = size(mask);
        hz = uint16(nz/2);
        wz = uint16((xyPixelSize/dz)*nz/dsFactor/2);
        
        mask(:, :, (hz-wz):(hz+wz)) = 1;
        mask = rotateFrame3D( ...
                mask, ...
                skewAngle * xyPixelSize/dz - angleBump, ...
                0.9*dz/xyPixelSize, ...
                'Crop', true, ...
                'outSize', size(mask) ...
            );
        mask = fftshift(mask);
        tiltedMask_Map = addTime(tiltedMask_Map, 'createMask', toc);
    
        % fftn
        tic;
        G = fftn(im_dsk);
        test4_Map = addTime(test4_Map, 'fftn', toc);
    
        G_full = gather(G);
    
        % masking
        tic;
        G = G .* mask;
        test4_Map = addTime(test4_Map, 'otfMask', toc);
    
        % ifftn
        tic;
        im_recon = real(ifftn(G));
        test4_Map = addTime(test4_Map, 'ifftn', toc);
    
        % rotate
        tic;
        im_dsr = rotateFrame3D( ...
                im_recon, ...
                skewAngle, ...
                1.0, ...
                'reverse', true, ...
                'Crop', true, ...
                'outSize', outSize ...
            );
        test4_Map = addTime(test4_Map, 'rotateFrame3D', toc);
    end
end

% plot
subplot(2,3,6);
imagesc(pwr(G_full, 1));
hold on;
mask = mip(fftshift(mask), 1);
bw = bwperim(mask);
[x,y] = meshgrid(1:size(mask,2), 1:size(mask,1));
x = x(bw(:));
y = y(bw(:));
scatter(x, y, '.');
hold off;
set(gca, 'colormap', hot);
title('Masked OTF');

subplot(2,3,3);
imagesc(mip(im_dsr,1));
axis image;
set(gca, 'colormap', parula);
title(test4_Name);

% release GPU
reset(D);

%% TEST 5: deskew -> interpolate-z -> rotate

test5_Map = containers.Map;
test5_Name = 'PetaKit5D: deskew, interpolate, rotate';

for j = 1:nRuns
    for i = 1:nTimePts
        imPath = pathList{i};
        
        % load image
        im_raw = single(permute(readtiff(imPath), [2,1,3]));
        
        % create gpuArrays
        tic;
        im_raw = gpuArray(im_raw);
        test5_Map = addTime(test5_Map, 'gpuArray-1', toc);
    
        % deskew
        tic;
        im_dsk = deskewFrame3D( ...
                im_raw, ...
                skewAngle, ...
                dz * dsFactor, ...
                xyPixelSize, ...
                'reverse', true ...
        );
        test5_Map = addTime(test5_Map, 'deskewFrame3D', toc);
    
        % manual interpolation
        tic;
        S = [1 0 0 0
            0 1 0 0
            0 0 dsFactor 0
            0 0 0 1];
        
        im_interp = imwarp(im_dsk, affine3d(S));
        test5_Map = addTime(test5_Map, 'imwarp(z)', toc);
    
        tic;
        im_dsr = rotateFrame3D( ...
                im_interp, ...
                skewAngle, ...
                1.0, ...
                'reverse', true, ...
                'Crop', true, ...
                'outSize', outSize ...
            );
         test5_Map = addTime(test5_Map, 'rotateFrame3D', toc);   
    end
end

subplot(2,3,4);
imagesc(mip(im_dsr, 1));
axis image;
title(test5_Name);

% release GPU
reset(D);

%%
fig = figure(1);
nexttile(1);
ax = gca; cla;

hold on;
[tf1, kl1] = timingBoxesPlot(test1_Map,      ax, 0.10, 0.05, test1_Name, colorMap);
[tf2, kl2] = timingBoxesPlot(test5_Map,      ax, 0.60, 0.05, test5_Name, colorMap);
[tf3, kl3] = timingBoxesPlot(test2_Map,      ax, 0.35, 0.05, test2_Name, colorMap);
[tf4, kl4] = timingBoxesPlot(tiltedMask_Map, ax, 0.90, 0.05, '',       colorMap);
[tf5, kl5] = timingBoxesPlot(test4_Map,      ax, 0.85, 0.05, test4_Name, colorMap);
hold off;

kl = unique([kl1, kl2, kl3, kl4, kl5]);

xlim([-0.001 Inf]);
ylim([0 1.05]);

box on;
set(ax, 'FontSize', 18, 'LineWidth', 2.0);

hl = [];
for key = kl
    key = key{1};
    hl = [hl; patch(NaN, NaN, colorMap(key))];
end

legend(hl, kl, 'location', 'southeast', 'fontsize', 18);

yticks([]);
xlabel('execution time / seconds');

%% functions

function [ out ] = addTime( map, key, t )

    try
        val = map(key);
        val.times = [val.times t];
    catch
        val = struct();
        val.times = t;
        
        keys = {};
        try
            keys = map('keys'); % maintain key order
            keys{end+1} = key;
        catch
            keys = {key};
        end

        map('keys') = keys;
    end
    
    map(key) = val;
    out = map;
end

function [ t_f, keyList ] = timingBoxesPlot( map, ax, y, dy, rName, color_map )
    
    bump = 0;
    x = bump;
    i = 1;

    text(ax, 0, y-dy/2, rName, ...
        'rotation', 0, 'FontSize', 22, 'Color', [0,0,0]);

    keyList = {};
    

    for key = map('keys')
        key = key{1};

        temp = strsplit(key, '-');
        temp = temp{1};
        keyList{end+1} = temp;

        t = median(map(key).times);
        % t = mean(map(key).times);
        color = color_map(key);

        p = [x y; x y+dy; x+t y+dy; x+t y];
        fill(ax, p(:,1), p(:,2), color);

        i = i + 1;
        x = x + t;
    end

    t_f = x - bump;
end

function [ out ] = ctrcrop( in, outSize )
    [sx, sy, sz] = size(in);
    hx = uint16(sx/2);
    hy = uint16(sy/2);
    hz = uint16(sz/2);
    
    dx = uint16(outSize(1)/2);
    dy = uint16(outSize(2)/2);
    dz = uint16(outSize(3)/2);

    out = in(...
            (hx-dx):(hx+dx-1), ...
            (hy-dy):(hy+dy-1), ...
            (hz-dz):(hz+dz-1) ...
        );
end
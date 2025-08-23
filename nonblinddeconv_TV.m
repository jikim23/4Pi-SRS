% 1) Add paths (if needed) so your functions are on MATLAB’s search path
% Change fileparts to location of example data and the averagePSF
folder = fileparts('C:\Users\jonik\OneDrive\Desktop\FuLab\CodeResources\4PiArbitraryPhase\FIgure3_Deconv\ROI1_bigCrop\');
addpath(genpath(folder));

% 2) Load your data
rawFile      = fullfile(folder,'roi1fordeconv.tif');
g_obs        = readTiffStack(rawFile);        % size [Nx×Ny×Nz] or [Nx×Ny×Nz×K]

psfFile       = fullfile(folder,'average_psf.tif');
psf        = readTiffStack(psfFile);
psf = psf / sum(psf(:));


%% Downsampling psfDan from zoom 6 to zoom 4 
% voxel dimensions for figure 3 A549 

% Original voxel sizes
vx_orig = 0.049;
vy_orig = 0.049;
vz_orig = 0.03;

% Target voxel sizes
vx_new = 0.072;
vy_new = 0.072;
vz_new = 0.03;

% Original array size
[nx, ny, nz] = size(psf);

% New grid size (downsampled)
nx_new = round((vx_orig / vx_new) * nx);
ny_new = round((vy_orig / vy_new) * ny);
nz_new = nz;  % unchanged

% Create original coordinate grid
[x, y, z] = ndgrid(1:nx, 1:ny, 1:nz);

% Create new coordinate grid
xq = linspace(1, nx, nx_new);
yq = linspace(1, ny, ny_new);
zq = 1:nz;  % no change

[xqg, yqg, zqg] = ndgrid(xq, yq, zq);

% Interpolate to new grid
psfzoom4 = interpn(x, y, z, psf, xqg, yqg, zqg, 'linear');

% Optional: Normalize PSF
psfzoom4 = psfzoom4 / sum(psfzoom4(:));


%%
f = double(g_obs);
f = f / max(f(:));
% --- your existing pre-step (optional RL warm start) ---
f = deconvlucy(f, psfzoom4, 10);

% --- params (add stop rules; keep your weights) ---
params.iter        = 5;        % hard cap (e.g., 30–50)
params.mu          = 20;
params.lambda_tv   = 0.005;
params.lambda_l0   = 0.03;

% optional background mask for sigma (logical same size as g)
% params.bgMask = bgMask;   

% stopping rules (all optional; these are sensible defaults)
params.stop.rel_change_tol = 1.e-4;   % ||Δo||/||o|| < tol
params.stop.window         = 3;      % …for M consecutive iters
params.stop.chi2_tol       = 0.11;   % |χ²/N - 1| < 5% (Gaussian)
params.stop.minIter        = 5;      % don0t stop too early
params.stop.data_delta_tol = 1e-4;   % |Δ data residual| < tol

% --- run (now returns diagnostics too) ---
[f, diag] = deconv_residual(f, psfzoom4, params);

% --- quick convergence plots ---
figure('Color','w','Position',[80 80 1000 320]);
subplot(1,3,1); semilogy(diag.iter, diag.data_rel,'-o','LineWidth',2); grid on; box on;
xlabel('iter'); ylabel('||H*o-g||/||g||'); title('data residual');

subplot(1,3,2); plot(diag.iter, diag.chi2,'-o','LineWidth',2); hold on; yline(1,'--k'); grid on; box on;
xlabel('iter'); ylabel('\chi^2 / N'); title('discrepancy');

subplot(1,3,3); semilogy(diag.iter, diag.rel_change,'-o','LineWidth',2); hold on;
yline(params.stop.rel_change_tol,'--'); grid on; box on; xlabel('iter'); ylabel('||Δo||/||o||'); title('relative change');





toc
beep
%%  Fast plotting deconvolution output to check quality

[nx, ny, nz] = size(f);
cx = round(nx/2);
cy = round(ny/2);
cz = round(nz/2);

figure('Position',[100 100 800 800]);

% 1) XY slice at z = center
ax1 = subplot(2,2,1);
imagesc(f(:,:,cz));
axis image off;
title(sprintf('XY @ z = %d', 48));
colormap(ax1,'jet');
colorbar(ax1,'Location','eastoutside');

% 2) XZ slice at y = center
ax2 = subplot(2,2,2);
imagesc(squeeze(f(:,cy,:)));
axis image off;
xlabel('X'); ylabel('Z');
title(sprintf('XZ @ y = %d', cy));
colormap(ax2,'jet');
colorbar(ax2,'Location','eastoutside');

% 3) YZ slice at x = center
ax3 = subplot(2,2,3);
imagesc(squeeze(f(cx,:,:))');
axis image off;
xlabel('Y'); ylabel('Z');
title(sprintf('YZ @ x = %d', cx));
colormap(ax3,'jet');
colorbar(ax3,'Location','eastoutside');

% 4) Axial (Z) profile at (cx,cy)
subplot(2,2,4);
plot((1:nz)*30, squeeze(f(cx,cy,:)), '-o','MarkerSize',4);
xlabel('Z position / nm');
ylabel('Intensity / A.U.');
title('Axial profile @ center');

sgt = sgtitle('Orthogonal Views & Center Line Profile');
sgt.FontSize = 14;



%% Write image back to tiff
writeStackToTiff(f , 'decvonvolutedtest5iter.tif');


function y = tv_solver_eq39(b, eta, maxIter)
%TV_SOLVER_EQ39  Solves y = argmin_y (eta/2)||b - y||^2 + TV(y)
%   using the alternating minimization algorithm of Wang et al. (2008)
%   References eq. (39) in the target paper.
%
%   y = tv_solver_eq39(b, eta, maxIter)
%     b       - input blurred image (2D or 3D array)
%     eta     - fidelity weight in (eta/2)||b-y||^2 + TV(y)
%     maxIter - number of AM iterations
%
%   Output:
%     y       - denoised/deconvolved image

% Ensure double
b = double(b);
dimsB = ndims(b);

% TV weight = alpha in 1/2||b-y||^2 + alpha*TV(y)
alpha = 1/eta;

% Initialize u and dual fields w = D u
y = b;  % initial u^0 = b

% Initialize w = zeros of same size but dims+1 for each gradient dir
if dimsB==2
    [ny,nx] = size(b);
    wx = zeros(ny,nx);
    wy = zeros(ny,nx);
elseif dimsB==3
    [ny,nx,nz] = size(b);
    wx = zeros(ny,nx,nz);
    wy = zeros(ny,nx,nz);
    wz = zeros(ny,nx,nz);
else
    error('Only 2D/3D supported');
end

% AM parameters
beta = 2*alpha;  % initial penalty
tau  = 2;        % beta multiplier per iteration

% Precompute Laplacian spectrum for FFT solve
if dimsB==2
    [wy_freq, wx_freq] = ndgrid(2*pi*(0:ny-1)/ny, 2*pi*(0:nx-1)/nx);
    Lap = (2-2*cos(wy_freq)) + (2-2*cos(wx_freq));
elseif dimsB==3
    [wy_freq, wx_freq, wz_freq] = ndgrid(2*pi*(0:ny-1)/ny, 2*pi*(0:nx-1)/nx, 2*pi*(0:nz-1)/nz);
    Lap = (2-2*cos(wy_freq)) + (2-2*cos(wx_freq)) + (2-2*cos(wz_freq));
end



for k = 1:maxIter
    %--- 1) u-update via FFT solve of (I - beta*Lap) u = b - beta*div(w)
    % Compute divergence of w
    if dimsB==2
        div_w = (wx - circshift(wx,[0,1])) + (wy - circshift(wy,[1,0]));
    else
        div_w = (wx - circshift(wx,[0,1,0])) + ...
                (wy - circshift(wy,[1,0,0])) + ...
                (wz - circshift(wz,[0,0,1]));
    end
    rhs = b - beta * div_w;
    % FFT solve
    denom = 1 + beta * Lap;
    y = real(ifftn(fftn(rhs) ./ denom));

    %--- 2) w-update via shrinkage: w = shrink(Dy, alpha/beta)
    if dimsB==2
        dyx = circshift(y,[0,-1]) - y;  % forward diff x
        dyy = circshift(y,[-1,0]) - y;  % forward diff y
        mag = sqrt(dyx.^2 + dyy.^2);
        shrink = max(mag - alpha/beta, 0) ./ (mag + eps);
        wx = shrink .* dyx;
        wy = shrink .* dyy;
    else
        dyx = circshift(y,[0,-1,0]) - y;
        dyy = circshift(y,[-1,0,0]) - y;
        dyz = circshift(y,[0,0,-1]) - y;
        mag = sqrt(dyx.^2 + dyy.^2 + dyz.^2);
        shrink = max(mag - alpha/beta, 0) ./ (mag + eps);
        wx = shrink .* dyx;
        wy = shrink .* dyy;
        wz = shrink .* dyz;
    end

    %--- 3) update beta
    beta = tau * beta;
end

% Return final y
end

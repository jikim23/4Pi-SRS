function [o, diag] = deconv_residual(g, h, params)
% deconv_hybrid_tvl0  Hybrid TV + L0‐gradient non‐blind deconvolution
% Usage:
%   [o, diag] = deconv_hybrid_tvl0(g, h, params)
% Inputs:
%   g                 observed 3D volume (Y×X×Z)
%   h                 PSF (same size as g or smaller)
%   params.iter       hard cap on outer iterations (e.g., 30–50)
%   params.mu         penalty weight
%   params.lambda_tv  TV weight
%   params.lambda_l0  L0‐gradient threshold weight
% Optional (for stopping rules):
%   params.bgMask               logical mask for background (same size as g)
%   params.stop.rel_change_tol  default 1e-3
%   params.stop.window          default 3
%   params.stop.chi2_tol        default 0.05
%   params.stop.minIter         default 5
%   params.stop.data_delta_tol  default 1e-4
% Outputs:
%   o                 deconvolved volume
%   diag              struct with per-iteration diagnostics

% -- ensure double & basic setup
g = double(g);  h = double(h);
imgSize = size(g);
H  = fftn(h, imgSize);        Hc = conj(H);
G  = fftn(g);

mu     = params.mu;
lam_tv = params.lambda_tv;
lam_l0 = params.lambda_l0;
Kmax   = params.iter;

% --- stopping defaults
st.rel_change_tol = getfielddef(params,'stop','rel_change_tol',1e-3);
st.window         = getfielddef(params,'stop','window',3);
st.chi2_tol       = getfielddef(params,'stop','chi2_tol',0.05);
st.minIter        = getfielddef(params,'stop','minIter',5);
st.data_delta_tol = getfielddef(params,'stop','data_delta_tol',1e-4);

% --- noise sigma for χ² (Gaussian); prefer bgMask if provided
if isfield(params,'bgMask') && ~isempty(params.bgMask)
    sigma = std(double(g(params.bgMask)));
else
    sigma = 1.4826 * mad(double(g(:)),1);   % robust MAD
end
sigma = max(sigma, eps);

% init
o = g;                 % warm start (or use your RL output before calling)
o_prev = [];
belowCnt = 0;

diag.iter       = [];
diag.data_rel   = [];  % ||H*o - g|| / ||g||
diag.chi2       = [];  % mean( (r/sigma).^2 )
diag.rel_change = [];  % ||o - o_prev|| / ||o_prev||

for k = 1:Kmax
  %---- 1) Data step: Wiener-ish inversion ----%
  O = fftn(o);
  Z = (Hc .* G + mu * O) ./ (abs(H).^2 + mu);
  o = real(ifftn(Z));

  %---- 2) TV denoise step ----%
  o = tv_solver_eq39(o, mu/lam_tv, 20);

  %---- 3) L0-gradient step ----%
  [dx, dy, dz] = gradient(o);
  mag   = sqrt(dx.^2 + dy.^2 + dz.^2);
  thr   = lam_l0 / mu;
  keep  = mag > thr;
  dx(~keep) = 0;  dy(~keep) = 0;  dz(~keep) = 0;

  % divergence (backward diffs via diff + padding)
  div_x = cat(1, dx(1,:,:),    diff(dx,1,1));
  div_y = cat(2, dy(:,1,:),    diff(dy,1,2));
  div_z = cat(3, dz(:,:,1),    diff(dz,1,3));
  divg  = div_x + div_y + div_z;

  %---- 4) Gradient-based correction ----%
  o = o - (lam_l0/mu) * divg;

  % non-negativity
  o(o<0) = 0;

  %================ Diagnostics & stopping =================
  % residuals
  g_pred   = real(ifftn(H .* fftn(o)));
  r        = g_pred - g;

  data_rel = norm(r(:)) / max(norm(g(:)), eps);
  chi2     = mean( (r(:)./sigma).^2 );

  if isempty(o_prev)
      rel_change = NaN;
  else
      rel_change = norm(o(:)-o_prev(:)) / max(norm(o_prev(:)), eps);
  end

  % store
  diag.iter(end+1)       = k;
  diag.data_rel(end+1)   = data_rel;
  diag.chi2(end+1)       = chi2;
  diag.rel_change(end+1) = rel_change;

  % update rolling counter for rel_change
  if ~isnan(rel_change) && rel_change < st.rel_change_tol
      belowCnt = belowCnt + 1;
  else
      belowCnt = 0;
  end

  % optional: also require data residual to be stalling
  data_delta_ok = (numel(diag.data_rel) < 2) || ...
                  (abs(diag.data_rel(end) - diag.data_rel(end-1)) < st.data_delta_tol);

  % stopping logic
  hit_rel = (belowCnt >= st.window) && (k >= st.minIter);
  hit_chi = (k >= st.minIter) && (abs(chi2 - 1) <= st.chi2_tol);

  if hit_rel || (hit_chi && data_delta_ok)
      % fprintf('Stopped at iter %d (rel-change/chi2 criteria).\n', k);
      break;
  end

  o_prev = o;
end
end

% --------- tiny helper ---------
function out = getfielddef(S,sub,field,def)
    if isfield(S,sub) && isfield(S.(sub),field) && ~isempty(S.(sub).(field))
        out = S.(sub).(field);
    else
        out = def;
    end
end

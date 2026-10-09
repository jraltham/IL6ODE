function feats = toxicityFeatures(t,x)

% -----------------------------
% Extract toxicity
% -----------------------------
X = x(:,6);

% -----------------------------
% Check for zero toxicity
% -----------------------------
if max(X) <= 0
    % No toxicity developed
    feats.Tonset = NaN;
    feats.Tres   = NaN;
    feats.Xmax   = 0;
    return
end

% -----------------------------
% Normalize toxicity
% -----------------------------
Xnorm = X / max(X);

% -----------------------------
% Onset: first time Xnorm > 50% of peak
% -----------------------------
idx_on = find(Xnorm > 0.5, 1, 'first');
if isempty(idx_on)
    feats.Tonset = t(end); % set to end of simulation if never reaches threshold
else
    feats.Tonset = t(idx_on);
end

% -----------------------------
% Peak toxicity
% -----------------------------
[~, iPeak] = max(Xnorm);
tPeak = t(iPeak);

% -----------------------------
% Resolution: first time after peak Xnorm drops below 50%
% -----------------------------
% Make sure we only search **after peak**
idx_postPeak = iPeak+1 : length(t);
idx_res = find(Xnorm(idx_postPeak) < 0.5, 1, 'first');

if isempty(idx_res)
    feats.Tres = t(end) - tPeak;  % resolution not reached → use end of simulation
else
    feats.Tres = t(idx_postPeak(idx_res)) - tPeak;
end

% -----------------------------
% Peak toxicity (un-normalized)
% -----------------------------
feats.Xmax = max(X);

end
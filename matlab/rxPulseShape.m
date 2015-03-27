function [ yr ] = rxPulseShape( yc )

inputSig = yc(:);

% do pulse-shaping & upsample
filtSpan = 30;       % Filter span in symbol durations
beta = 0.35;        % Roll-off factor
sps = 2;

% generate the coefficients for the filter (same exact code as the
% txPulseShape_80211b
t = (-filtSpan:filtSpan)/2;
c = -4.*beta./sps .* ( cos((1+beta).*pi.*t) + ...
    sin((1-beta).*pi.*t) ./ (4.*beta.*t) ) ...
    ./ (pi .* ((4.*beta.*t).^2 - 1));
% set the middle tap properly
c(filtSpan+1) = -1 ./ (pi.*sps) .* (pi.*(beta-1) - 4.*beta );
% normalize power
c = c / sqrt(sum(c.^2));

inputSig = [inputSig; zeros(filtSpan*sps/2,1)];

yr = filter(c,1,inputSig);
yr = yr(filtSpan+1:end);
% convert back to row vector b/c i like those
yr = conj(yr');


end


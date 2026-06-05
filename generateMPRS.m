function [u_seq, u_fun, t_vec, seed] = generateMPRS( ...
    t0, tf, step, levels, Td_min, Td_max, seed)

if nargin < 7 || isempty(seed)
    rng('shuffle');             
    s = rng;
    seed = s.Seed;
else
    rng(seed);                   
end

t_vec = t0:step:tf;
N = length(t_vec);
u_seq = zeros(2, N);

k = 1;
while k <= N

    % Random dwell time
    Td = Td_min + (Td_max - Td_min) * rand;
    Nk = max(1, round(Td / step));

    % Random levels for two input channels
    idx = randi(numel(levels), [2, 1]);
    u_level = levels(idx);
    u_level = u_level(:);   % FORCE column vector (2×1)

    % Assign over dwell interval
    k_end = min(N, k + Nk - 1);
    u_seq(:, k:k_end) = repmat(u_level, 1, k_end - k + 1);

    k = k_end + 1;
end

u_fun = @(t) interp1(t_vec', u_seq', t, 'previous', 'extrap')';

end

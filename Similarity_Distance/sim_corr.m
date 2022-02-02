function s = sim_corr(r_a, r_b)

s = sum(r_a.*r_b) / sqrt(sum(r_a .* r_a)) / sqrt(sum(r_b .* r_b));

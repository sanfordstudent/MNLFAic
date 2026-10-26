# Small example (a real analysis would use more persons, items and draws);
# item 4 has residual DIF for English learners
sim <- simulate_mnlfa_ic(
  person_data = data.frame(el = rep(0:1, 150)), person = ~ el,
  item_data = data.frame(domain = rep(c("number", "geometry", "data"), length.out = 16)),
  item = ~ domain, slopes = rep(1.2, 16), intercepts = seq(-1, 1, length.out = 16),
  impact = -0.3, residual_dif = matrix(replace(numeric(16), 4, 1)), seed = 6)
fit <- mnlfa_ic(sim$responses, sim$person_data, ~ el,
                item_data = sim$item_data, item = ~ domain,
                control = mnlfa_control(M = 200, n_pilot = 100, seed = 1))

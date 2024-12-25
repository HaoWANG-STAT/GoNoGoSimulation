SD <- 3084.06
mean <- 4356.67
CV <- SD/mean
CV
sigma <- sqrt(log(CV^2+1))
sigma


CV <- sqrt(exp(sigma^2)-1)
CV
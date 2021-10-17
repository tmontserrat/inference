# Useful statistical inference functions
# One-sample
# Function to calculate a confidence interval using a z-critical value
create_conf_interval_normal <- function(mean.sample, sd.sample, sample.size, percent) {
  critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
  lwr.bound <- mean.sample-(critical.z*(sd.sample/sqrt(sample.size)))
  upr.bound <- mean.sample+(critical.z*(sd.sample/sqrt(sample.size)))
  return(c(lwr.bound, upr.bound))
}

# Function to calculate a confidence interval using a t-critical value
create_conf_interval_student <- function(mean.sample, sd.sample, sample.size, 
                                         percent, alternative="two.sided") {
  if (alternative=="two.sided") {
    df <- sample.size - 1
    critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
    lwr.bound <- mean.sample-(critical.t*(sd.sample/sqrt(sample.size)))
    upr.bound <- mean.sample+(critical.t*(sd.sample/sqrt(sample.size)))
    return(c(lwr.bound, upr.bound))
  }
  else if (alternative=="lower.one.sided") {
    df <- sample.size - 1
    critical.t <- qt(percent/100, df)
    lwr.bound <- mean.sample-(critical.t*(sd.sample/sqrt(sample.size)))
    upr.bound <- mean.sample+(critical.t*(sd.sample/sqrt(sample.size)))
    return(c(Inf, upr.bound))
  }
  else if (alternative=="upper.one.sided") {
    df <- sample.size - 1
    critical.t <- qt(percent/100, df)
    lwr.bound <- mean.sample-(critical.t*(sd.sample/sqrt(sample.size)))
    upr.bound <- mean.sample+(critical.t*(sd.sample/sqrt(sample.size)))
    return(c(lwr.bound, Inf))
  }
  
}

# Function that simulate samples from a given variable from a dataset and return
# a table with the proportions of those samples that capture the correct mean
check_capture_mean <- function(replicates, sample.size, mu, percent, dataset, variable) {
  # Empty vectors to store results
  sample.means <- c()
  m <- c() 
  contains.mu <- c()
  
  # Store results in a vector
  set.seed(2017)
  
  for (k in 1:replicates) {
    # Take a random sample
    samples <- dataset[sample(nrow(dataset), sample.size), variable]
    sample.means[k] <- mean(as.numeric(samples), na.rm=TRUE)
    
    # Calculate margin of error
    m[k] <- qnorm(percent/100 + (1-(percent/100))/2)*(sd(samples, na.rm=TRUE)/sqrt(sample.size))
    
    # Check if the confidence interval captures the population mean
    upr.bound <- sample.means[k] + m[k] 
    lwr.bound <- sample.means[k] - m[k]
    if (mu <= upr.bound & mu >= lwr.bound) {
      contains.mu[k] <- TRUE
    } else {
      contains.mu[k] <- FALSE
    }
  }
  
  return(prop.table(table(contains.mu)))
}


# Function that simulate samples from a given variable from a dataset and return
# a table with the proportions of those samples that it's t-statistic (for the mean)
# is in the rejection zone
check_in_rejection_zone <- function(sample_mean,
                                    mu,
                                    sample_sd,
                                    sample_size,
                                    alpha,
                                    replicates, 
                                    data, 
                                    variable) {
  # Empty vector to store results
  t_in_rejection <- c("logical", replicates)
  for (k in 1:replicates) {
    # Complete that for the variable of interest
    data.complete <- data[!is.na(data[, variable]), ]
    # Take a sample of size = 100
    # 100 random rows
    random.rows <- sample(nrow(data.complete), sample_size)
    data.sample <- data.complete[random.rows, ]
    sample_mean <- mean(data.sample[, variable])
    sample_sd <- sd(data.sample[, variable])
    
    # t-statistic
    t_st <- t_statistic(sample_mean, mu, sample_sd, sample_size)
    # Critical t
    t_critical <- abs(qt(alpha/2, df))
    
    # Check if t-statistic is in the rejection zone
    t_in_rejection[k] <- abs(t_st) > t_critical
    
  }
  return(prop.table(table(t_in_rejection)))
}

# Return the t statistic given a sample mean, population mean,
# sample standard deviation and sample size
t_statistic <- function(sample_mean, population_mean, sd_sample, n){
  return((sample_mean - population_mean)/(sd_sample/sqrt(n)))
} 


# Adapted from Ref: https://gist.github.com/jrnold/6799152
# Function that generates a normal curve and color the area of interest
normal_prob_area_plot <- function(lb, ub, mean = 0, sd = 1, 
                                  limits = c(mean - 3 * sd, mean + 3 * sd),
                                  type = "lower_upper_bounds") {
  if(type == "lower_upper_bounds") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax <- seq(xmin, xmax, length.out = 100)
    area <- data.frame(x = areax, ymin = 0, ymax = dnorm(areax, mean = mean, sd = sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y))
           + geom_area(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y), fill = "steelblue")
           + geom_ribbon(data = area, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "red")
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "lower_bound") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax <- seq(limits[1], xmin, length.out = 100)
    area <- data.frame(x = areax, ymin = 0, ymax = dnorm(areax, mean = mean, sd = sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y))
           + geom_area(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y), fill = "steelblue")
           + geom_ribbon(data = area, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "red")
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "upper_bound") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax <- seq(xmax, limits[2], length.out = 100)
    area <- data.frame(x = areax, ymin = 0, ymax = dnorm(areax, mean = mean, sd = sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y))
           + geom_area(data.frame(x = x, y = dnorm(x, mean = mean, sd = sd)),
                       mapping = aes(x = x, y = y), fill = "steelblue")
           + geom_ribbon(data = area, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "red")
           + scale_x_continuous(limits = limits) +
             theme_bw())
  }
}

# Function that generate a t student distribution curve and color
# the interested area. Useful for visualizing p-values
student_prob_area_plot <- function(lb_t_st, ub_t_up, df, 
                                   limits = c(-4, 4),
                                   type = "lower_upper_bounds") {
  if (type == "lower_upper_bounds") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 100)
    areax2 <- seq(limits[2], xmax, length.out = 100)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dt(areax1, df))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dt(areax2, df))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dt(x, df)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area1, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "lower_bound") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 100)
    areax2 <- seq(limits[2], xmax, length.out = 100)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dt(areax1, df))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dt(areax2, df))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dt(x, df)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area1, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "upper_bound") {
    x <- seq(limits[1], limits[2], length.out = 100)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 100)
    areax2 <- seq(limits[2], xmax, length.out = 100)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dt(areax1, df))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dt(areax2, df))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dt(x, df)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  }
}



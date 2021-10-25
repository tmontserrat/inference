# Copyright 2021 Tomàs Montserrat Ayuso
# Useful statistical inference functions

# Load helper functions
source("helper_inference_functions.R")

# Function to calculate a standard error
calculate_standard_error <- function(sd, sample.size) {
  return(sigma/sqrt(sample_size))
}

# Function to calculate a confidence interval
create_conf_interval <- function(mean.sample, sd=1, sample.size, 
                                  percent, alternative="two.sided",
                                  type="student", proportion=FALSE) {
  # Calculate the confidence interval using a t-statistic
  if (type == "student") {
    # For means
    if (proportion == FALSE) {
      if (alternative == "two.sided") {
        conf_interval_two_sided(mean.sample, sd, sample.size, percent)
      }
      else if (alternative == "lower.one.sided") {
        conf_interval_one_sided_left(mean.sample, sd, sample.size, percent)
      }
      else if (alternative=="upper.one.sided") {
        conf_interval_one_sided_right(mean.sample, sd, sample.size, percent)
      }
    } 
    # For proportions
    else if (proportion == TRUE) {
      if (alternative=="two.sided") {
        conf_interval_two_sided_prop(mean.sample, sample.size, percent)
      } 
      else if (alternative=="lower.one.sided") {
        conf_interval_one_sided_left_prop(mean.sample, sample.size, percent)
      } 
      else if (alternative=="upper.one.sided") {
        conf_interval_one_sided_right_prop(mean.sample, sample.size, percent)
      }
    }
  } # Calculate the confidence interval using a z-statistic 
  else if (type == "normal") {
    # For means
    if (proportion == FALSE) {
      if (alternative == "two.sided") {
        conf_interval_two_sided_normal(mean.sample, sd, sample.size, percent)
      } 
      else if (alternative == "lower.one.sided") {
        conf_interval_one_sided_left_normal(mean.sample, sd, sample.size, percent)
      } 
      else if (alternative == "upper.one.sided") {
        conf_interval_one_sided_right_normal(mean.sample, sd, sample.size, percent)
      }
    }
    # For proportions
    else if (proportion == TRUE) {
      if (alternative == "two.sided") {
        conf_interval_two_sided_prop(mean.sample, sample.size, percent)
      }
      else if (alternative == "lower.one.sided") {
        conf_interval_one_sided_left_normal(mean.sample, sd, sample.size, percent)
      }
      else if (alternative == "upper.one.sided") {
        conf_interval_one_sided_right_normal(mean.sample, sample.size, percent)
      }
    }
  }
}

# Function that returns the sample size needed for a given precision. 
# Needs an initial proportion and sample size as reference. For proportions
get_sample_size <- function(prop, precision, percent) {
  critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
  return(critical.z^2 * prop*(1-prop)/(precision)^2)
}

# Function to calculate a margin of error
calculate_margin_error <- function(sd, sample.size, percent, type="student") {
  if (type == "student") {
    df <- sample.size - 1
    critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
    margin_error <- critical.t*(sd/sqrt(sample.size))
    return(margin_error)
  } else if (type == "normal") {
    critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
    margin_error <- critical.z*(sd/sqrt(sample.size))
    return(margin_error)
  }
}

# Function that calculates the degrees of freedom for two samples
calculate_df_two_sample <- function(s1, s2, n1, n2) {
  df.nom <- (((s1^2)/n1) + ((s2^2)/n2))^2
  df.denom.1 <- (((s1^2)/n1)^2)/(n1-1)
  df.denom.2 <- (((s2^2)/n2)^2)/(n2-1)
  return(df.nom / (df.denom.1+df.denom.2))
}

# Function to calculate a confidence interval using a t-critical value. 
create_conf_interval_two_sample <- function(x1.bar, x2.bar, s1, s2, n1, n2, percent) {
  df <- calculate_df_two_sample(s1, s2, n1, n2)
  critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
  lwr.bound <- (x1.bar-x2.bar) - (critical.t*(sqrt( ((s1^2)/n1) + ((s2^2)/n2)) ))
  upr.bound <- (x1.bar-x2.bar) + (critical.t*(sqrt( ((s1^2)/n1) + ((s2^2)/n2)) ))
  return(c(lwr.bound, upr.bound))
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

check_p_values <- function(control.n, control.mean, control.sigma, 
                           treatment.n, treatment.mean, treatment.sigma,
                           alpha, replicates) {
  p.values <- vector("numeric", replicates)
  for (k in 1:replicates) {
    control <- rnorm(n=control.n, mean=control.mean, sd=control.sigma)
    treatment <- rnorm(n=treatment.n, mean=treatment.mean, sd=treatment.sigma)
    p.values[k] <- t.test(control, treatment, 
                          conf.level=(1-alpha), paired=FALSE,
                          alternative="two.sided", mu=0)$p.value
  }
  reject <- p.values <= alpha
  return(prop.table(table(reject)))
}

# Function that performs t.tests, between two simulated normal samples and
# check how many of them are significant. Useful to study the type II error
check_typeii_error <- function(control.n, control.mean, control.sigma,
                               treatment.n, treatment.mean, treatment.sigma,
                               alpha, replicates) {
  p.values <- vector("numeric", replicates)
  for (i in 1:replicates) {
    control <- rnorm(control.n, mean=control.mean, sd=control.sigma)
    treatment <- rnorm(treatment.n, mean=treatment.mean, sd=treatment.sigma)
    p.values[i] <- t.test(control, treatment, 
                          alternative="two.sided", mu=0, conf.level=1-alpha)$p.value
  }
  reject <- p.values <= alpha
  return(prop.table(table(reject)))
}

# Return the t statistic given a sample mean, population mean,
# sample standard deviation and sample size
t_statistic <- function(sample_mean, population_mean, sd_sample, n){
  return((sample_mean - population_mean)/(sd_sample/sqrt(n)))
} 

# Return the t statistic given in a independent two-sample test.
t_statistic_two_samples <- function(x1.bar, x2.bar, s1, s2, n1, n2){
  return((x1.bar - x2.bar)/(sqrt( ((s1^2)/n1) + ((s2^2)/n2)) ))
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
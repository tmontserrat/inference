# Copyright 2021 Tom?s Montserrat Ayuso
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
        conf_interval_two_sided_prop_normal(mean.sample, sample.size, percent)
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
# Function to calculate a confidence interval for two sample. 
create_conf_interval_two_sample <- function(x1.bar, x2.bar, s1=1, s2=1, n1, n2, 
                                            percent, proportions=FALSE) {
  if (proportions == FALSE) {
    df <- calculate_df_two_sample(s1, s2, n1, n2)
    critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
    se <- sqrt( ((s1^2)/n1) + ((s2^2)/n2))
    lwr.bound <- (x1.bar-x2.bar) - (critical.t*se)
    upr.bound <- (x1.bar-x2.bar) + (critical.t*(sqrt( ((s1^2)/n1) + ((s2^2)/n2)) ))
    return(c(lwr.bound, upr.bound))
  }
  else if (proportions == TRUE) {
    diff <- x1.bar - x2.bar
    radicand.1 <- (x1.bar*(1-x1.bar))/n1
    radicand.2 <- (x2.bar*(1-x2.bar))/n2
    se <- sqrt(radicand.1 + radicand.2)
    critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
    lwr.bound <- diff - (critical.z*se)
    upr.bound <- diff + (critical.z*se)
    return(c(lwr.bound, upr.bound))
  }
}

calculate_p_hat_pooled <- function(x1, x2, n1, n2) {
  return((x1+x2) / (n1+n2))
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

# Return the standardize test statistic given a sample mean, population mean,
# sample standard deviation and sample size
standardize_mean <- function(sample_mean, population_mean, sd_sample, n){
  return((sample_mean - population_mean)/(sd_sample/sqrt(n)))
} 

# Return the t statistic given in a independent two-sample test.
standardize_two_samples <- function(x1.bar, x2.bar, s1, s2, n1, n2){
  return((x1.bar - x2.bar)/(sqrt( ((s1^2)/n1) + ((s2^2)/n2)) ))
} 

# Return z statistic for a one sample proportion.
z_statistic_one_sample_prop <- function(p.hat, p.zero, alpha) {
  nom <- p.hat-p.zero
  denom <- sqrt((p.zero*(1-p.zero))/n)
  z.star <- nom/denom
  return(z.star)
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

normal_prob_area_plot_p_val <- function(lb, ub, mean = 0, sd = 1, 
                                        limits = c(mean - 3 * sd, mean + 3 * sd),
                                        type = "lower_upper_bounds") {
  if (type == "lower_upper_bounds") {
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dnorm(areax1, mean, sd))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dnorm(areax2, mean, sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean, sd)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area1, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "lower_bound") {
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dnorm(areax1, mean, sd))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dnorm(areax2, mean, sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean, sd)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area1, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "upper_bound") {
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb, limits[1])
    xmax <- min(ub, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dnorm(areax1, mean, sd))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dnorm(areax2, mean, sd))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dnorm(x, mean, sd)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
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
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
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
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
    area1 <- data.frame(x = areax1, ymin = 0, ymax = dt(areax1, df))
    area2 <- data.frame(x = areax2, ymin = 0, ymax = dt(areax2, df))
    return(ggplot()
           + geom_line(data.frame(x = x, y = dt(x, df)),
                       mapping = aes(x = x, y = y))
           + geom_ribbon(data = area1, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
           + scale_x_continuous(limits = limits) +
             theme_bw())
  } else if (type == "upper_bound") {
    x <- seq(limits[1], limits[2], length.out = 1000)
    xmin <- max(lb_t_st, limits[1])
    xmax <- min(ub_t_up, limits[2])
    areax1 <- seq(limits[1], xmin, length.out = 1000)
    areax2 <- seq(limits[2], xmax, length.out = 1000)
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

f_prob_area_plot <- function(f, df1, df2, 
                             limits = c(0, 6)) {
  x <- seq(limits[1], limits[2], length.out = 1000)
  xmax <- min(f, limits[2])
  areax2 <- seq(limits[2], xmax, length.out = 1000)
  area2 <- data.frame(x = areax2, ymin = 0, ymax = df(areax2, df1, df2))
  return(ggplot()
         + geom_line(data.frame(x = x, y = df(x, df1, df2)),
                     mapping = aes(x = x, y = y))
         + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
         + scale_x_continuous(limits = limits) +
           theme_bw())
}

chi_squared_prob_area_plot <- function(q, df, 
                                       limits = c(0, df*4)) {
  x <- seq(limits[1], limits[2], length.out = 1000)
  xmax <- min(q, limits[2])
  areax2 <- seq(limits[2], xmax, length.out = 1000)
  area2 <- data.frame(x = areax2, ymin = 0, ymax = dchisq(areax2, df))
  return(ggplot()
         + geom_line(data.frame(x = x, y = dchisq(x, df)),
                     mapping = aes(x = x, y = y))
         + geom_ribbon(data = area2, mapping = aes(x = x, ymin = ymin, ymax = ymax), fill = "dodgerblue", alpha=0.4)
         + scale_x_continuous(limits = limits) +
           theme_bw())
}

# Useful function to shadow a specific region in a binomial distribution
bar_plot_binomial_prob <- function(min.prob, 
                                   max.prob,
                                   n, 
                                   prob, 
                                   x.labels.by=1 , 
                                   color.success="coral", 
                                   color.failure="dodgerblue") {
  # Create data frame
  df = data.frame(x=0:n, y=dbinom(0:n, n, prob))
  df$group <- df$x
  
  # Create the borders for the shadowed region
  prob_X <- c(min.prob:max.prob)
  
  # Create two groups for each observation
  df$group[df$x %in% prob_X] = "color"
  df$group[!(df$x %in% prob_X)] = "no_color"
  
  # Create the plot
  bar_plot_binomial = ggplot(df, aes(x, y, fill = group)) +
    geom_bar(stat = "identity",
             col = "black") +
    scale_fill_manual(values = c("no_color" = "dodgerblue", "color" = "coral"), guide = "none") +
    labs(
      x = "X",
      y = "Probability"
    ) +
    scale_x_continuous(breaks = round(seq(min(df$x), 
                                          max(df$x), 
                                          by = x.labels.by),1)) +
    theme_bw()
  
  return(bar_plot_binomial)
}

# Function to calculate the t-statistic to evaluate the linear relationship
t_linear_cor <- function(cor, n) {
  return(cor/(sqrt((1-cor^2)/(n-2))))
}
  
# Function to calculate the f-statistic to evaluate ANOVA
f_statistic <- function(ssb, ssw, k, n) {
  return( (ssb/(k-1)) / (ssw/(n-k)) )
}

plot_of_means <- function(dataframe, factor, response, percent=95, labels=c("Factor", "Means")) {
  # Split the dataframe by the levels of the factor
  factor_levels <- levels(dataframe[, factor])
  response_split <- split(dataframe[, response], dataframe[, factor])
  response_means <- unlist(lapply(response_split, mean))
  response_sds <- unlist(lapply(response_split, sd))
  response_ns <- unlist(lapply(response_split, length))
  
  # Empty vectors to store the data
  lower_bound_conf_int <- c()
  upper_bound_conf_int <- c()
  # Iterate over the differents levels of the factor
  for (i in 1:length(factor_levels)) {
    # Store the information of interest
    lower_bound_conf_int[i] <- create_conf_interval(response_means[i],
                                                    response_sds[i],
                                                    response_ns[i],
                                                    percent)[1]
    upper_bound_conf_int[i] <- create_conf_interval(response_means[i],
                                                    response_sds[i],
                                                    response_ns[i],
                                                    percent)[2]
  }
  # Create the dataframe with the data
  d <- data.frame(type=factor_levels, 
                  mean=response_means, 
                  lower_bound=lower_bound_conf_int,
                  upper_bound=upper_bound_conf_int)
  
  plot <- ggplot(data=d, aes(x=type, y=mean, group=1)) +
    geom_point(shape=16, size=3) + 
    geom_line() +
    geom_errorbar(aes(ymin=lower_bound, ymax=upper_bound),
                  width=0.2, alpha=0.4, linetype="dashed") +
    labs(x=labels[1], y=labels[2]) +
    theme_bw()
  
  return(plot)
}

compare_means_dot_plot <- function(dataframe, factor, response) {
  # Retrieve the levels
  factor_levels <- levels(dataframe[, factor])
  # Split the response variable by the levels
  response_split <- split(dataframe[, response], dataframe[, factor])
  # Store the grouped means and median in vectors
  response_means <- unlist(lapply(response_split, mean))
  response_medians <- unlist(lapply(response_split, median))
  # Store the global mean and medians in vectors
  response_mean <- mean(dataframe[, response])
  response_median <- median(dataframe[, response])
  # Generate the visualization
  plot <- ggplot(dataframe, aes_string(factor, response)) + 
    geom_point(alpha=0.4) +
    annotate("segment", y = response_mean, 
             yend = response_mean,
             x = 0.5, xend = length(factor_levels)+0.5,
             colour = "blue", size=0.5, linetype="dashed") +
    annotate("segment", y = response_median, 
             yend = response_median,
             x = 0.5, xend = length(factor_levels)+0.5,
             colour = "green", size=0.5, linetype="dashed")
  for (i in 1:length(response_means)) {
    plot <- plot +
      annotate("segment", y = response_means[i], 
               yend = response_means[i],
               x = 0.7+i-1, xend = 1.3+i-1,
               colour = "blue", size=0.5) +
      annotate("segment", y = response_medians[i],
               yend = response_medians[i],
               x = 0.7+i-1, xend = 1.3+i-1,
               colour = "green", size=0.5)
  }
  
  plot <- plot + theme_bw()
  return(plot)
}

median_conf_int_bootstrap <- function(data, iterations=100, bins=25, percent=95) {
  alpha <- 1 - (percent/100)
  medians <- replicate(iterations, {
    # Creamos una nueva muestra
    ggt_sample <- sample(data, length(data), replace=TRUE)
    
    # Calculamos la mediana
    ggt_median_sample <- median(ggt_sample)
  })
  
  # Calculamos la desviación estándar de la distribución de las medianas
  sd_medians <- sd(medians)
  
  # Tomamos los percentiles 25 y 97.5 para calcular el intervalo de confianza
  conf_interval <- quantile(medians, c(alpha/2, 1-(alpha/2)))
  
  # Creamos el histograma de la distribución
  medians_df <- data.frame(medians=medians)
  histogram <- ggplot(medians_df, aes(x=medians)) +
    geom_histogram(bins=bins, color="black", fill="steelblue") +
    labs(y="Frequency", x="Medians") +
    theme_bw()
  
  return(list("medians"=medians, 
              "sd"=sd_medians, 
              "conf_int"=conf_interval, 
              "histogram"=histogram))
}

median_conf_int_conventional <- function(data, percent=95) {
  # Sample size
  n <- length(data)
  
  # Quantile of interest
  q <- 0.5
  
  # Critical Z
  z_critico <- calculate_critical_z(percent=percent)
  
  # Indices for the confidence interval
  j <- round(n*q - z_critico*sqrt(n*q*(1-q)))
  k <- round(n*q + z_critico*sqrt(n*q*(1-q)))
  
  # Sort values
  sorted_data <- sort(data)
  
  # Build the confidence interval
  median_conf_int <- c(sorted_data[j], sorted_data[k])
  
  return(median_conf_int)
}
  
  
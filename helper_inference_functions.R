# Copyright 2021 Tomàs Montserrat Ayuso
# Helper functions for statistical_inference_functions.R

# Confidence intervals for means
conf_interval_two_sided <- function(mean.sample, sd, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
  lwr.bound <- mean.sample-(critical.t*(sd/sqrt(sample.size)))
  upr.bound <- mean.sample+(critical.t*(sd/sqrt(sample.size)))
  return(c(lwr.bound, upr.bound))
}

conf_interval_two_sided_normal <- function(mean.sample, sd, sample.size, percent) {
  critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
  lwr.bound <- mean.sample-(critical.z*(sd/sqrt(sample.size)))
  upr.bound <- mean.sample+(critical.z*(sd/sqrt(sample.size)))
  return(c(lwr.bound, upr.bound))
}

conf_interval_one_sided_left <- function(mean.sample, sd, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100, df)
  upr.bound <- mean.sample+(critical.t*(sd/sqrt(sample.size)))
  return(c(Inf, upr.bound))
}

conf_interval_one_sided_left_normal <- function(mean.sample, sd, sample.size, percent) {
  critical.z <- qnorm(percent/100)
  upr.bound <- mean.sample+(critical.z*(sd/sqrt(sample.size)))
  return(c(Inf, upr.bound))
}

conf_interval_one_sided_right <- function(mean.sample, sd, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100, df)
  lwr.bound <- mean.sample-(critical.t*(sd/sqrt(sample.size)))
  return(c(lwr.bound, Inf))
}

conf_interval_one_sided_right_normal <- function(mean.sample, sd, sample.size, percent) {
  critical.z <- qnorm(percent/100)
  lwr.bound <- mean.sample-(critical.z*(sd/sqrt(sample.size)))
  return(c(lwr.bound, Inf))
}

# Confidence intervals for proportions
conf_interval_two_sided_prop <- function(mean.sample, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100 + (1-(percent/100))/2, df)
  lwr.bound <- mean.sample-(critical.t*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  upr.bound <- mean.sample+(critical.t*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(lwr.bound, upr.bound))
}

conf_interval_two_sided_prop_normal <- function(mean.sample, sample.size, percent) {
  critical.z <- qnorm(percent/100 + (1-(percent/100))/2)
  lwr.bound <- mean.sample-(critical.z*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  upr.bound <- mean.sample+(critical.z*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(lwr.bound, upr.bound))
}

conf_interval_one_sided_left_prop <- function(mean.sample, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100, df)
  upr.bound <- mean.sample+(critical.t*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(Inf, upr.bound))
}

conf_interval_one_sided_left_prop_normal <- function(mean.sample, sample.size, percent) {
  critical.z <- qnorm(percent/100)
  upr.bound <- mean.sample+(critical.z*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(Inf, upr.bound))
}

conf_interval_one_sided_right_prop <- function(mean.sample, sample.size, percent) {
  df <- sample.size - 1
  critical.t <- qt(percent/100, df)
  lwr.bound <- mean.sample-(critical.t*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(lwr.bound, Inf))
}

conf_interval_one_sided_right_prop_normal <- function(mean.sample, sample.size, percent) {
  critical.t <- qnorm(percent/100)
  lwr.bound <- mean.sample-(critical.t*(sqrt((mean.sample*(1-mean.sample))/sample.size)))
  return(c(lwr.bound, Inf))
}
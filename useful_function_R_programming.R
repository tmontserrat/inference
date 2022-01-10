# Useful functions in R programming

#' To factor
#' Converts to factor a variable in a dataframe
#' @param dataframe The whole dataframe
#' @param variables_to_convert A vector with the column names to convert

to_factor <- function(dataframe, variables_to_convert) {
    for (variable in variables_to_convert) {
      dataframe[, variable] <-  factor(dataframe[[variable]])
    }
    return(dataframe)
}

#' Evaluate normality
#' Create a histogram and a qq-plot
#' @param dataframe The dataframe with the data
#' @param variable The variable to analyze
#' @param bins The number of bins the histogram will have
#' @param title Should the plots have titles?
#' @param title_hist The title of the histogram
#' @param title_qq The title of the qq-plot

evaluate_normality <- function(dataframe, 
                               variable, 
                               bins=30,
                               title=FALSE, title_hist='', title_qq='') {
  if (title) {
    histogram <- ggplot(dataframe, aes_string(x=variable)) +
      geom_histogram(bins=bins, color="black", fill="steelblue") +
      labs(y="Frecuencia", title=title_hist) +
      theme_bw()
    
    qq_plot <- ggplot(dataframe, aes_string(sample=variable)) + 
      stat_qq(color="coral") +
      geom_qq_line() +
      labs(y="Cuantiles de la muestra", x="Cuantiles teóricos", title=title_qq) +
      theme_bw()
    
    return(grid.arrange(histogram, qq_plot, nrow=1))
  }
  else {
    histogram <- ggplot(dataframe, aes_string(x=variable)) +
      geom_histogram(bins=bins, color="black", fill="steelblue") +
      labs(y="Frecuencia") +
      theme_bw()
    
    qq_plot <- ggplot(dataframe, aes_string(sample=variable)) + 
      stat_qq(color="coral") +
      geom_qq_line() +
      labs(y="Cuantiles de la muestra", x="Cuantiles teóricos") +
      theme_bw()
    
    return(grid.arrange(histogram, qq_plot, nrow=1))
  }
  
}

#' Calculate variance ratio
#' Calculates the variance of a variable grouped by the levels of a factor
#' and give the ratio between the maximum and minimum
#' @param response The vector of a dataframe of the response variable
#' @param factor The vector of a dataframe of the factor variable

calculate_variance_ratio <- function(response, factor) {
  # Calculamos la varianza de albumina1 por nivel de fib4_3
  var_fib <- tapply(response, factor, var)
  
  # Calculamos la ratio entre la máxima varianza y la mínima
  max_var <- max(var_fib)
  min_var <- min(var_fib)
  var_ratio <- max_var/min_var
  names(var_ratio) <- c("ratio")
  results <- list(var_fib, var_ratio)
  names(results) <- c("variances", "ratio")
  return(results)
}
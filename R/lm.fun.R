lm.fun <- function(y, x, varNames = c()) {
  
  
  # check and generate variable names
  
  hasNames <- FALSE
  yName <- "Y"
  xNames <- c("(Intercept)")
  
  if (length(varNames) == length(x)+1) {
    hasNames <- TRUE
    yName <- varNames[1]
    xNames  <- c(xNames, varNames[2:length(varNames)])
  } else {
    warning("not all or no variable names are provided, default is used.")
    for (i in 1:length(x)) {
      xNames <- c(xNames, paste0("X",i))
    }
  }
  
  
  # checks for y
  
  if (!(is.vector(y)) | is.list(y)) {
    stop("y must be of type vector()")
  } else if (length(y) < 2) {
    stop("y is empty or contains only one value")
  } else {
    for (i in 1:length(y)) {
      if (!(is.numeric(y[i]))) {
        stop("elements of y must be of type 'numeric'")
      }
    }
  }
  
  
  # containers to save information about variables of type factor
  
  ref.categories <- c()  # saves reference categories for factorial variables
  fac.vars.num <- c()  # saves variable numbers of factor variables
  fac.vars.verb <- c()  # saves names of factor variables
  fac.data <- data.frame(dummyrow = rep(666, length(y))) # saves observations of factor variables as numeric
  cnames.fac <- c()  # saves the column names for the data frame 'fac.data'
  
  
  # checks for x
  
  if (!(is.list(x)) | length(x) < 1) {
    
    stop("x must be a list of vectors (one vector at least)")
    
  } else {
    
    for (i in 1:length(x)) {
      
      if (!(is.vector(x[i]))) {
        
        stop("x must be a list of vectors")
      }
      
      if (length(y) != length(x[[i]])) {
        
        stop("vector y and the vectors of list x must be of same length")
      }
      
      if (!(is.numeric(x[[i]]))) {
        
        if (!(is.factor(x[[i]]))) {
          
          stop("vectors of x must be of type 'numeric' or 'factor'")
          
        } else {  # if variable is of type factor
          
          if (length(levels(x[[i]])) < 2) {
            
            stop("factor variables must consist of more than one level")
            
          } else {  # transform a factor variable with z levels to z-1 numeric variables
            
            fac.vars.num <- c(fac.vars.num, i)
            fac.vars.verb <- c(fac.vars.verb, xNames[i+1])
            ref.categories <- c(ref.categories, levels(x[[i]])[1])
            
            for (j in 1:(length(levels(x[[i]]))-1)) {
              
              currLevel <- levels(x[[i]])[j+1]
              dataVar <- rep(0, length(x[[i]]))
              dataVar[which(x[[i]] == currLevel)] <- 1
              fac.data <- cbind(fac.data, dataVar)
              cnames.fac <- c(cnames.fac, paste0(xNames[i+1], "_", currLevel))
            }
          }
        }
      }
    }
  }
  
  
  # generate formula string
  
  formString <- paste0(yName, " ~ ")
  for (i in 2:length(xNames)) {
    if (i == 2) {
      formString <- paste0(formString, xNames[i])
    } else {
      formString <- paste0(formString, " + ", xNames[i])
    }
  }
  
  
  # generate matrices
  
  Y <- as.matrix(y)
  X <- data.frame(rep(1, nrow(Y)))
  for (i in 1:length(x)) {
    X <- cbind(X, x[i])
  }
  X <- cbind(X, fac.data[,-1])  # add the dummy variables for the factor variables
  X <- X[,-(fac.vars.num+1)]  # delete the factor variables
  X <- as.matrix(X)
  xNames <- c(xNames, cnames.fac)
  xNames <- xNames[-(fac.vars.num+1)]
  colnames(X) <- xNames
  
  
  
  # define n and k
  n <- length(y)
  k <- ncol(X)
  
  ## first we estimate the betas
  beta.hat = solve(t(X) %*% X) %*% t(X) %*% Y
  
  # we store it into a data frame
  resCoeff = as.data.frame(beta.hat)
  colnames(resCoeff) = c("Estimate")
  
  ## next we calculate the standard errors
  # calculate predicted values
  y.hat = X %*% beta.hat
  
  # calculate residuals
  epsilon = y - y.hat
  
  # calculate variance-covariance-matrix
  vcov = 1 / (n-k) * as.numeric(t(epsilon) %*% epsilon) * solve(t(X) %*% X)
  
  # calculate the standard errors
  s.e. = sqrt(diag(vcov)) # the variances are the diagonal elements
  
  # add stand errors to result
  resCoeff$'Std. Error' = s.e.
  
  # calculate t value and add them to the result
  t = beta.hat / s.e.
  resCoeff$'t value' = t
  
  # calculate p value and add them to the result
  p = 2*pt(abs(t), df=n-k,lower.tail= FALSE)
  resCoeff$'Pr(>|t|)' = p
  
  # add stars
  stars = ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,"."," "))))
  resCoeff$' ' = stars
  
  # calculate residual quartils
  res.quartils = quantile(epsilon,p=seq(0,1,length.out=5))
  names(res.quartils) <- c("Min", "1Q", "Median", "3Q", "Max")
  
  res.error = sqrt(sum((y.hat-Y)^2)/(n-k))
  
  # calculate multiple rsquared
  mul.rsq = sum((y.hat - mean(Y))^2) / sum((Y-mean(Y))^2)
  
  # calculate adjusted rsquared
  adj.rsq = 1 - (1 - mul.rsq)*(n-1) / (n-k)
  
  # calculate F-statistics and its p-value
  f =  mul.rsq/(k - 1) / ((1 - mul.rsq) / (n-k))
  f.p = pf(f, k-1, n-k, lower.tail=F)
  
  # round numbers
  resCoeff$Estimate <- round(resCoeff$Estimate, 6)
  resCoeff$`Std. Error` <- round(resCoeff$`Std. Error`, 6)
  resCoeff$`t value` <- round(resCoeff$`t value`, 3)
  resCoeff$`Pr(>|t|)` <- ifelse(resCoeff$`Pr(>|t|)`<2e-16, "< 2e-16",
                                format(resCoeff$`Pr(>|t|)`, scientific = TRUE, digits = 3))
  res.quartils <- round(res.quartils, 4)
  res.error <- round(res.error, 3)
  mul.rsq <- round(mul.rsq, 5)
  adj.rsq <- round(adj.rsq, 5)
  f <- round(f, 2)
  f.p <- ifelse(f.p<2.2e-16, "< 2.2e-16",
                format(resCoeff$`Pr(>|t|)`, scientific = TRUE, digits = 4))
  
  
  # print result to console
  
  writeLines(paste0("\nFormula:\n", formString))
  
  if (length(fac.vars.verb) > 0 & length(fac.vars.verb) == length(ref.categories)) {
    writeLines("\nDummies:")
    for (i in 1:length(fac.vars.verb)) {
      writeLines(paste0(fac.vars.verb[i], ": reference category = '", ref.categories[i], "'"))
    }
  }
  writeLines("\nResiduals:")
  print(res.quartils)
  writeLines("\nCoefficients:")
  print(resCoeff)
  writeLines(paste0("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n",
                    "Residual standard error: ", res.error, " on ", n-k, " degrees of freedom", "\n",
                    "Multiple R-squared:  ", mul.rsq, "    Adjusted R-squared:  ", adj.rsq, "\n",
                    "F-statistic: ", f, " on ", k-1, " and ", n-k, " DF,  p-value: ", f.p, "\n"))
}
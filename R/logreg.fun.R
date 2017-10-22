l.link <- function(x){
  
  p <- exp(x) / (1 + exp(x))
  p <- as.matrix(p)
  
  return(p)
}

loglike <- function(theta, Y, X) {
  
  betas   <- theta
  betas <- as.matrix(betas)
  
  # Log lik function
  mu <- X %*% betas
  logl <- sum(t(Y)%*%log(l.link(mu)) + t((1 - Y))%*%log(1-l.link(mu)))
  
  return(logl)
}

logreg.fun <- function(formula, data) {
  
  data <- as.data.frame(data)
  
  Y = as.data.frame(data[,as.character(formula[2])])
  names = trimws(unlist(strsplit(as.character(formula[3]),"\\+")))
  
  allVars <- FALSE
  
  # select covariates
  if (length(names) == 1) {
    
    if (names == ".") {
      
      allVariables <- TRUE
      
      pos <- 0
      for (i in 1:length(colnames(data))) {
        if (as.character(formula[2]) == colnames(data)[i]) {
          pos <- i
          break
        }
      }
      if (pos > 0) {
        X = data[,-pos]
        names <- colnames(X)
      } else {
        stop("problems with variable name finding")
      }
    }
  }
  
  if (!allVars) {
    X = data[,names]
  }
  
  X = data[,names]
  X = cbind(rep(1,nrow(X)),X)
  colnames(X)[1] <- "(Intercept)"
  
  newDat <- cbind(Y,X)
  rowsBefore <- nrow(newDat)
  newDat <- na.omit(newDat)
  numNAs <- rowsBefore  - nrow(newDat)
  
  Y <- newDat[,1]
  X <- newDat[,-1]
  
  # containers to save information about variables of type factor
  ref.categories <- c()  # saves reference categories for factor variables
  fac.vars.num <- c()  # saves variable numbers of factor variables
  fac.vars.verb <- c()  # saves names of factor variables
  fac.data <- data.frame(dummyrow = rep(666, length(Y))) # saves observations of factor variables as numeric
  cnames.fac <- c()  # saves the column names for the data frame 'fac.data'
  
  # check covariates
  for (i in 1:ncol(X)) {
    
    if (!is.numeric(X[,i])) {
      
      if (!is.factor(X[,i])) {
        
        stop("used variables in 'data' need to be of type numeric or factor")
        
      } else {
        
        if (length(levels(X[,i])) < 2) {
          
          stop("factor variables must consist of more than one level")
          
        } else {  # transform a factor variable with z levels to z-1 numeric variables
          
          fac.vars.num <- c(fac.vars.num, i)
          fac.vars.verb <- c(fac.vars.verb, colnames(X)[i])
          ref.categories <- c(ref.categories, levels(X[,i])[1])
          
          for (j in 1:(length(levels(X[,i]))-1)) {
            
            currLevel <- levels(X[,i])[j+1]
            dataVar <- rep(0, length(X[,i]))
            dataVar[which(X[,i] == currLevel)] <- 1
            fac.data <- cbind(fac.data, dataVar)
            cnames.fac <- c(cnames.fac, paste0(colnames(X)[i], "_", currLevel))
          }
        }
      }
    }
  }
  
  xNames <- colnames(X)
  X <- cbind(X, fac.data[,-1])  # add the dummy variables for the factor variables
  X <- X[,-(fac.vars.num)]  # delete the factor variables
  X <- as.matrix(X)
  xNames <- c(xNames, cnames.fac)
  xNames <- xNames[-(fac.vars.num)]
  colnames(X) <- xNames
  
  # check dependent variable
  if (!is.numeric(Y)) {
    
    if (!is.factor(Y)) {
    
      stop("dependent variable must be of type numeric or factor")
      
    } else {
      
      if (length(levels(Y)) != 2) {
        
        stop("dependent variable must be dichotomous")
        
      } else {
        
        Y <- as.numeric(Y)-1
      }
    }
  } else {
    
    if (length(unique(Y)) != 2) {
      
      stop("dependent variable must be dichotomous")
    }
  }
  
  Y <- as.matrix(Y)
  
  n <- length(Y)
  k <- ncol(X)
  
  # optimize our model
  res <- optim(rep(0, ncol(X)),
               loglike,
               Y=(Y),X=(X),
               control=list(fnscale=-1),
               hessian = T)
  
  # optimize null model (to get 'null deviance')
  res0 <- optimize(f=loglike, interval = c(0,1), maximum = TRUE, Y=(Y), X=(X[,1]))
  
  # get residuals and their deviance
  beta.hat <- as.matrix(res$par)
  e <- Y - X %*% beta.hat
  res.quartils = quantile(e,p=seq(0,1,length.out=5))
  names(res.quartils) <- c("Min", "1Q", "Median", "3Q", "Max")
  
  # get coefficients and add them to the result
  result <- res$par
  result <- as.data.frame(result)
  colnames(result) <- c("Estimate")
  rownames(result) <- colnames(X)
  
  # get standard errors and add them to the result
  fisher_info <- solve(-res$hessian)
  se <- sqrt(diag(fisher_info))
  result$'Std. Error' <- se
  
  # calculate z values and add them to the result
  z = res$par / se
  result$'z value' = z
  
  # calculate p value and add them to the result
  p = 2*pt(abs(z), df=n-k,lower.tail= FALSE)
  result$'Pr(>|z|)' = p
  
  # add stars
  stars = ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,"."," "))))
  result$' ' = stars
  
  # get null deviance
  nulldev <- - 2 * res0$objective
  
  # get residual deviance
  resdev <- -2 * res$value
  
  # calculate AIC
  aic <- -2*res$value + 2*k
  
  # round all the values
  res.quartils <- round(res.quartils, 4)
  result$Estimate <- round(result$Estimate, 6)
  result$`Std. Error` <- round(result$`Std. Error`, 6)
  result$`z value` <- round(result$`z value`, 3)
  result$`Pr(>|z|)` <- ifelse(result$`Pr(>|z|)`<2e-16, "< 2e-16", ifelse(result$`Pr(>|z|)`<0.00001, format(result$`Pr(>|z|)`, scientific = TRUE, digits = 3), round(result$`Pr(>|z|)`,6)))
  nulldev <- round(nulldev, 1)
  resdev <- round(resdev, 1)
  aic <- round(aic, 1)
  
  
  # write results to console
  writeLines(paste0("\nFormula:\n", formula[2], " ", formula[1], " ", formula[3]))
  if (length(fac.vars.verb) > 0 & length(fac.vars.verb) == length(ref.categories)) {
    writeLines("\nDummies:")
    for (i in 1:length(fac.vars.verb)) {
      writeLines(paste0(fac.vars.verb[i], ": reference category = '", ref.categories[i], "'"))
    }
  }
  writeLines("\nDeviance Residuals:")
  print(res.quartils)
  writeLines("\nCoefficients:")
  print(result)
  writeLines(paste0("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n",
                    "    Null deviance: ", nulldev," on ", n-1, " degrees of freedom\n",
                    "Residual deviance: ", resdev, " on ", n-k, " degrees of freedom\n",
                    "  (", numNAs, " observations deleted due to missingness)\n",
                    "AIC: ", aic, "\n"))
}

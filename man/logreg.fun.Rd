\name{logreg.fun}
\alias{logreg.fun}

\title{
logistic regression
}
\description{
The function calculates a logistic regression and prints the results to the console.
}
\usage{
logreg.fun(formula, data)
}
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
}
  \item{data}{
a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
}
}
\details{
The function prints the logistic regression results to the console; including information about the coefficients, standard errors, z- and p-values and the AIC. However, the function does not return. The printed output is similar to the one that is provided by summary() for glm() models of family 'binomial'.
}
\value{
The function does not return a value, but prints the logistic regression results to the console.
}

\author{
  Saro Gibilisco, \email{sarogiba@bluewin.ch}
  Rolf Badat, \email{rolf.badat@uzh.ch}
}

\examples{
logreg.fun(W.Hnd ~ Sex + Clap + Fold, MASS::survey)
}

\keyword{ logistic regression }

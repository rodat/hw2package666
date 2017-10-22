\name{lm.fun}
\alias{lm.fun}

\title{
ordinary least squares (OLS) regression
}
\description{
The function calculates an OLS regression and prints the results to the console.
}
\usage{
lm.fun(y, x, varNames = c())
}
\arguments{
  \item{y}{
dependent variable which must be of type numeric
}
  \item{x}{
list of independent variables which must be vectors and of type numeric or factor.
}
  \item{varNames}{
OPTIONAL a vector of all the variable names: the first variable name in the vector must be the one of the dependent variable, followed by the names of the independent variables. default: c()
}
}
\details{
The function prints the OLS regression results to the console; including information about the coefficients, standard errors, t- and p-values, residuals, R squared and the F-Statistics. However, the function does not return. The printed output is similar to the one that is provided by summary() for lm() models.
}
\value{
  The function does not return a value, but prints the OLS regression results to the console.
}

\author{
  Saro Gibilisco, \email{sarogiba@bluewin.ch}
  Rolf Badat, \email{rolf.badat@uzh.ch}
}

\examples{
dat <- MASS::survey
dat <- na.omit(dat)
lm.fun(dat$Pulse, list(dat$Age, dat$Sex, dat$Smoke), c("Pulse", "Age", "Sex", "Smoke"))
}

\keyword{ lm }

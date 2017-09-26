# This is a script demonstrating how to compute prediction interval from a GAM.
# The predict.lm() of R has support for such intervals, but predict.gam() does not.
# predict.gam() does optionally compute a standard error of the individual fits
# and this can be used to construct the intervals, but the means is not immediately
# apparent.
#
#
# We will start by loading the package for GAM support.
library(mgcv)
# Next, set a sample size large enough that we should see strong similarities when the
# smoothing of a GAM is turned on to avoid vagaries of sampling.
n <- 1000
# Now create some data.
X <- rnorm(n, mean=100, sd=20) # simulate some X
Y <- 0.6*X + rnorm(n, sd=5) # simulate some Y
plot(X,Y) # show strong linear relation

# STEP 1: Start from a know location. Fit a line, extract the prediction interval as if
# the sample were the new data being used for prediction.
LM <- lm(Y~X) #  # fit the model, conventional OLS
PL <- predict(LM, interval = c("prediction"), se.fit=TRUE)
PL <- as.data.frame(PL) # this casts all of the results in columnar format.
print(head(PL))

# Next compute the quantile of the t-distribution give the 0.95 prediction interval
QT <- abs(qt(0.05/2, (n-2))) # quantile of t-distribution
# A review of Helsel and Hirsch (2002) shows that we are basically interested in the
# computation  SIGMA*QT*sqrt(1+HAT) where SIGMA is the residual standard error, QT as
# described already, and sqrt(1+HAT) for the component unique to a data points position
# in parameter space. The "1+" is the feature for prediction interval and not confidence
# interval.

# STEP 2: Now that we have the fit from OLS, let us compute the GAM without smoothing
# for the same data. This is an MLE solution though to my understanding and not OLS.
GM <- gam(Y~X) # fit the model, maximum likelihood (MLE).
PG <- predict(GM, se.fit=TRUE) # we only get the standard error of fit, not prediction
# intervals
PG <- as.data.frame(PG) # this casts all of the results in columnar format.
GM$sigma <- sqrt(summary(GM)$scale) # the SIGMA of a GAM is a variance, though it is
# labeled as a "scale"
PG$residual.scale <- GM$sigma # just assign the SIGMA to the residual.scale to mimic the
# PL nomenclature
print(head(PG))

# Next, let us extract the degrees of freedom of the GAM. The value is 2, which matches
# that of the linear model because the GAM is fitting a linear model: slope + intercept.
Gdf <- sum(GM$edf) # value is 2 (!!!!! no smoothing !!!!!)
# Next, let us compute the equivalent QT.
GQT <- abs(qt(0.05/2, (n-Gdf))) # quantile of t-distribution
# Inspection show that QT and GQT are the same as they must be.

# Now heuristic manipulation of the prediction intervals from lm() against the
# se.fit from gam() yield that we can isolate the sqrt(1+HAT) term in the prediction
# interval like so:
sqrt1phat_GAM <- sqrt(1+(PG$se.fit/PG$residual.scale)^2) # LOOK HERE!
# And the user can compare lm() prediction intervals to those now developed for gam().
PG$lwr <- PG$fit - GQT*PG$residual.scale*sqrt1phat_GAM
PG$upr <- PG$fit + GQT*PG$residual.scale*sqrt1phat_GAM
summary(PG$lwr - PL$fit.lwr) # basically 1 part in 1E-12
summary(PG$upr - PL$fit.upr) # The values are all zero though not quite zero because
# MLE and OLS are numerically distinct pathways to the fit but result in effectively
# the same fit.

# This the equivalent of a hat value (leverage) is the ratio of the standard error
# variance to the residual variance. Let us check that for the linear model.
head((PL$se.fit/PL$residual.scale)^2)
head(hatvalues(LM))
# We see that the are the same. Thus, the "leverage" of a GAM is
head((PG$se.fit/PG$residual.scale)^2)

# STEP 2: Now that we have shown that OLS and GAM are the same and how to access what we
# need, let us make the jump to smoothing term. The data are very linear so the smooth
# itself should have limited impact.
GMs <- gam(Y~s(X)) # fit the model, maximum likelihood (MLE) and cross validiation on the
# smooth
PGs <- predict(GMs, se.fit=TRUE) # we only get the standard error of fit, not prediction
# intervals
PGs <- as.data.frame(PGs) # this casts all of the results in columnar format.
GMs$sigma <- sqrt(summary(GMs)$scale) # the SIGMA of a GAM is a variance, though it is
# labeled as a "scale"
PGs$residual.scale <- GMs$sigma # just assign the SIGMA to the residual.scale to mimic the
#  PL nomenclature
print(head(PGs))

# Next, let us extract the degrees of freedom of the GAM. The value is 2, which matches
# that of the linear model because the GAM is fitting a linear model: slope + intercept.
Gdfs <- sum(GMs$edf) # value is 5.912325 (!!!!! smoothing !!!!!). Note that the df have
# gone up because of the "knot points" in the smoothing. It is possible that for very high
# correlation and n that a test run produces a Gdfs exactly equal to 2---the smooth is a
# line. Try rerunning.
GQTs <- abs(qt(0.05/2, (n-Gdfs))) # quantile of t-distribution

# We proceed with computation of the prediction interval and compare to those from the
# linear model. We see that they are indeed different, but the mean is about zero.
sqrt1phat_GAMs <- sqrt(1+(PGs$se.fit/PGs$residual.scale)^2)
PGs$lwr <- PGs$fit - GQTs*PGs$residual.scale*sqrt1phat_GAMs
PGs$upr <- PGs$fit + GQTs*PGs$residual.scale*sqrt1phat_GAMs

summary(PGs$lwr - PL$fit.lwr)
summary(PGs$upr - PL$fit.upr) # Below my test case
#LOWER:    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#      -0.79110 -0.21290 -0.07269  0.01100  0.05497  2.68800
#
#UPPER:    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#      -0.36900 -0.24810 -0.10780 -0.01100  0.01448  4.19100

# Let this compare the 1:1 relation (as maybe) between the linear model and smoothed GAM
# fit along with the lower and upper prediction limits. Things look pretty good.
plot(PL$fit.fit, PGs$fit, log="yx")
points(PL$fit.lwr, PGs$lwr, col=2)
points(PL$fit.upr, PGs$upr, col=4)


# Finally, this helper function could be used for prediction limit computation.

"gamIntervals" <-
function(gam_predicts_with_se.fit, gam=NULL,
         interval=c("none", "confidence", "prediction"), level=0.95, ...) {
   # Demo: library(mgcv)
   #       X <- 2*pi*(1:360)/360 # simulate some X
   #       Y <- 1.6*sin(X) + 40*cos(X) + rnorm(length(X), sd=12)
   #     GAM <- gam(Y~s(X)); PGAM <- predict(GAM, se.fit=TRUE)
   #     PGAM <- gamIntervals(PGAM, gam=GAM)
   #     print(head(PGAM))
   #     print(head(PGAM$leverage)); print(head(GAM$hat)) # see they are the value
   # plot(GAM$hat, (PGAM$se.fit/PGAM$residual.scale)^2)
   # Compare what the GAM says its leverage values are to back computed.
   # The plot() only works because predict() called back on the actuall model.
   if(class(gam)[1] != "gam") {
      warning("need the actual GAM model too via the 'gam' argument")
      return()
   }
   z <- as.data.frame(gam_predicts_with_se.fit)
   if(! any(names(z) == "se.fit")) {
      warning("need gam predictions with se.fit=TRUE passed for 'gam_predicts_with_se.fit'")
      return()
   }
   interval <- match.arg(interval)
   sum.gam <- summary(gam); n <- sum.gam$n # summary.gam() and the sample size
   z$residual.scale <- sigma <- sqrt(sum.gam$scale) # residual standard error
   df <- n-sum(gam$edf)           # total degrees of freedom
   QT <- abs(qt((1-level)/2, df)) # will do the +/- separately
   z$leverage <- (z$se.fit/sigma)^2
   if(interval == "none") {
      z$lwr <- z$upr <- NA
   } else {
        one <- ifelse(interval == "confidence", 0, 1)
        tmp <- sqrt(one+z$leverage)
      z$lwr <- z$fit - sigma*QT*tmp
      z$upr <- z$fit + sigma*QT*tmp
   }
   print(head(z))
   attr(z, "interval")                  <- interval
   attr(z, "level")                     <- level
   attr(z, "t-dist_degrees_of_freedom") <- df
   return(z)
}

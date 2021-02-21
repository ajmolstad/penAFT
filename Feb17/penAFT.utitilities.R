penAFT.predict <- function(fit, Xnew, lambda = NULL){
  if(!is.matrix(Xnew)) {
    stop("Xnew must be a matrix of dimension n_new x p")
  }
  if (class(fit)!="penAFT" & class(fit)!="penAFT.cv") {
    stop("Input 'fit' must be a model fit from penAFT or penAFT.cv")
  }
  if (class(fit) == "penAFT") {
    if(is.null(lambda) | !any(fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    } else {
      if(fit$standardize){
        Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))/(rep(1, dim(Xnew)[1])%*%t(fit$X.sd))
      } else {
        if(fit$center){
          Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))
        } else {
          Xpred <- Xnew
        }
      }
      s <- which(fit$lambda == lambda)
      preds <- Xpred%*%as.matrix(fit$beta[,s])
    }
    
  } else {
    if(any(fit$full.fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    }
    if(is.null(lambda)){
      s <- which.min(fit$cv.err.linPred)
    }
    fit <- fit$full.fit
    if(fit$standardize){
      Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))/(rep(1, dim(Xnew)[1])%*%t(fit$X.sd))
    } else {
      if(fit$center){
        Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))
      } else {
        Xpred <- Xnew
      }
    }
    preds <- Xpred%*%as.matrix(fit$beta[,s])
  }
  return(preds)
}



penAFT.coef <- function(fit, lambda = NULL){
  
  if (class(fit)!="penAFT" & class(fit)!="penAFT.cv") {
    stop("Input 'fit' must be a model fit from penAFT or penAFT.cv")
  }
  if (class(fit) == "penAFT") {
    if(is.null(lambda) | !any(fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    } else {
      if(fit$standardize){
        s <- which(fit$lambda == lambda)
        beta.out <- (rep(1, dim(Xnew)[1])%*%t(fit$X.sd))*as.matrix(fit$beta[,s])
        warning("Coefficients are for centered predictors!")
      } else {
        if(fit$center){
          s <- which(fit$lambda == lambda)
          beta.out <- as.matrix(fit$beta[,s])
          warning("Coefficients are for centered predictors!")
        } else {
          s <- which(fit$lambda == lambda)
          beta.out <- as.matrix(fit$beta[,s])
        }
      }
    }
    
  } else {
    if(any(fit$full.fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    }
    if(is.null(lambda)){
      s <- which.min(fit$cv.err.linPred)
    }
    fit <- fit$full.fit
    if(fit$standardize){
      beta.out <- (rep(1, dim(Xnew)[1])%*%t(fit$X.sd))*as.matrix(fit$beta[,s])
      warning("Coefficients are for centered predictors!")
    } else {
      if(fit$center){
        beta.out <- as.matrix(fit$beta[,s])
        warning("Coefficients are for centered predictors!")
      } else {
        beta.out <- as.matrix(fit$beta[,s])
      }
    }
  }
  return(list("beta" = beta.out, "mean.adjustment" = fit$X.mean))
  
}

penAFT.plot <- function(fit){
  
  if(class(fit)!="penAFT.cv"){
    stop("Input 'fit' must be a model fit using penAFT.cv.")
  }
  library(ggplot2)
  data <- data.frame(
    "log10lambda" = log10(fit$full.fit$lambda), 
    "linPred" = fit$cv.err.linPred,
    "ObjErr" = colMeans(fit$cv.err.obj)
  )
  # A few constants
  t1 <- "#69b3a2"
  t2 <- rgb(0.2, 0.6, 0.9, 1)
  
  p1 <- ggplot(data, aes(x=log10lambda)) +
    geom_line( aes(y=linPred/max(linPred), size = 1.2), color=t1) + 
    geom_line( aes(y=ObjErr/max(ObjErr), size = 1.2), color=t2) + 
    scale_y_continuous(name = "Linear predictor CV error", sec.axis = sec_axis(~.*1, name="Average within-fold CV error")) + 
    theme_light() +
    theme(
      axis.title.y = element_text(color = t1, size=13),
      axis.title.y.right = element_text(color = t2, size=13)
    ) + ggtitle("Relative cross-validation errors") + xlab(expression(log[10](lambda)))
  
  return(p1)
}
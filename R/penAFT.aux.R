penAFT.predict <- function(fit, Xnew, lambda = NULL){

  if(!is.matrix(Xnew)) {
    if(length(Xnew) == length(fit$X.mean)){
      Xnew <- matrix(Xnew, nrow=1)
    } else {
      stop("Xnew must be a matrix of dimension n_new x p")
    }
  }

  if (!inherits(fit, "penAFT") & !inherits(fit, "penAFT.cv")) {
    stop("Input 'fit' must be a model fit from penAFT or penAFT.cv")
  }

  if (inherits(fit, "penAFT")) {
    if(is.null(lambda) | !any(fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    } else {
      if(fit$standardize){
        Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))/(rep(1, dim(Xnew)[1])%*%t(fit$X.sd))
      } else {
        Xpred <- Xnew
      }
      s <- which(fit$lambda == lambda)
      preds <- Xpred%*%as.matrix(fit$beta[,s])
    }

  } else {
    
    if(is.null(lambda)){
      s <- min(which(fit$cv.err.linPred == min(fit$cv.err.linPred)))
    } else {
      if(!any(fit$full.fit$lambda == lambda)){
        stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
      } else {
        s <- which(fit$full.fit$lambda == lambda)
      }
    }

    fit <- fit$full.fit
    if(fit$standardize){
      Xpred <- (Xnew - rep(1, dim(Xnew)[1])%*%t(fit$X.mean))/(rep(1, dim(Xnew)[1])%*%t(fit$X.sd))
    } else {
      Xpred <- Xnew
    }
    preds <- Xpred%*%as.matrix(fit$beta[,s])
  }
  return(preds)
}



penAFT.coef <- function(fit, lambda = NULL){

  if (!inherits(fit, "penAFT") & !inherits(fit, "penAFT.cv")) {
    stop("Input 'fit' must be a model fit from penAFT or penAFT.cv")
  }
  if (inherits(fit, "penAFT")) {
    if(is.null(lambda) | !any(fit$lambda == lambda)){
      stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
    } else {
      if (fit$standardize) {
        s <- which(fit$lambda == lambda)
        beta.out <- (1/fit$X.sd)*as.matrix(fit$beta[,s])
      } else {
        s <- which(fit$lambda == lambda)
        beta.out <- as.matrix(fit$beta[,s])
      }
    }

  } else {
    
    if (is.null(lambda)) {
      s <- min(which(fit$cv.err.linPred == min(fit$cv.err.linPred)))
    } else {
      if(!any(fit$full.fit$lambda == lambda)){
        stop("Must supply input 'lambda' equal to element of penAFT$lambda, or use penAFT.cv for model fitting.")
      } else {
        s <- which(fit$full.fit$lambda == lambda)
      }
    }
    fit <- fit$full.fit
    if (fit$standardize) {
      beta.out <- (1/fit$X.sd)*as.matrix(fit$beta[,s])
    } else {
      beta.out <- as.matrix(fit$beta[,s])
    }
  }
  return(list("beta" = beta.out))

}




penAFT.plot <- function(fit){

  if (!inherits(fit, "penAFT.cv")) {
    stop("Input 'fit' must be a model fit using penAFT.cv.")
  }

  dat <- data.frame(
    "log10lambda" = log10(fit$full.fit$lambda),
    "linPred" = fit$cv.err.linPred,
    "ObjErr" = colMeans(fit$cv.err.obj),
    "sesObjErr" = apply(fit$cv.err.obj, 2, sd)/sqrt(dim(fit$cv.err.obj)[1])
  )
  
  # colors for both y-axes
  t1 <- "dodgerblue3"
  t2 <- "black"
  log10lambda <- dat$log10lambda
  ObjErr <- dat$ObjErr
  sesObjErr <- dat$sesObjErr
  linPred <- dat$linPred

  p1 <- ggplot(dat, aes(x=log10lambda)) +
      geom_ribbon(aes(ymin = ObjErr - sesObjErr, ymax = ObjErr  + sesObjErr), fill = "grey70", alpha=0.3) +
    geom_line( aes(y=ObjErr), color=t2,alpha=0.8) +
    geom_line( aes(y=linPred), color=t1,alpha=0.8) +
    scale_y_continuous(name = "Linear predictor CV error", sec.axis = sec_axis(~.*1, name="Within-fold CV error")) +
    theme_bw() +
    theme(axis.title.y = element_text(color = t1, size=13),
      axis.title.y.right = element_text(color = t2, size=13)) + ggtitle("Cross-validation errors") + xlab(expression(log[10](lambda))) + theme(plot.title = element_text(hjust = 0.5))+
    geom_vline(xintercept = log10(fit$full.fit$lambda)[which.min(fit$cv.err.linPred)], linetype="dotted",
                color = t1) +
    geom_vline(xintercept = log10(fit$full.fit$lambda)[min(which(ObjErr <= min(ObjErr + sesObjErr)))], linetype="dotted",
                color = t2)

  return(p1)
}






penAFT.trace <- function(fit, groupNames = NULL){
  
  if(!inherits(fit, "penAFT.cv") & !inherits(fit, "penAFT")){
    stop("Input 'fit' must be a model fit using penAFT.cv or penAFT.")
  }
  
  if (inherits(fit, "penAFT")) {
    fit$full.fit <- fit
  }

  beta <- fit$full.fit$beta
  
  if(is.null(fit$full.fit$groups)){
    
    t1 <- "dodgerblue3"
    t2 <- "black"
    
    if(inherits(fit, "penAFT.cv"){
      out <- as.matrix(fit$full.fit$beta)
        dat2 <- data.frame(
          "log10lambda" = log10(fit$full.fit$lambda),
          "ObjErr" = colMeans(fit$cv.err.obj),
          "sesObjErr" = apply(fit$cv.err.obj, 2, sd)/sqrt(dim(fit$cv.err.obj)[1])
        )
        out <- as.matrix(fit$full.fit$beta)
        dat <- data.frame(
          "values" = c(out[which(rowSums(abs(out)) > 0),]),
          "log10lambda" = rep(log10(fit$full.fit$lambda), each = sum(rowSums(abs(out)) != 0)),
          "predictor" = rep(1:sum(rowSums(abs(out)) != 0), length(fit$full.fit$lambda))
        )
        values <- dat$values
        log10lambda <- dat$log10lambda
        predictor <- dat$predictor
        
        p1 <- ggplot(dat, aes(x=log10lambda, y = values, color=as.factor(predictor))) +
          theme_bw() + theme(legend.position="") + ylab(expression(hat(beta)[j])) +
          geom_line() + geom_vline(xintercept = log10(fit$full.fit$lambda)[which.min(fit$cv.err.linPred)], linetype="dotted",
                                   color = t1) + geom_vline(xintercept = log10(fit$full.fit$lambda)[min(which(dat2$ObjErr <=  min(dat2$ObjErr + dat2$sesObjErr)))], linetype="dotted",
                                                            color = t2) +  ggtitle("Elastic net trace plot") +
          xlab(expression(log[10](lambda))) + theme(plot.title = element_text(hjust = 0.5))
    } else {
      
      out <- as.matrix(fit$full.fit$beta)
      dat <- data.frame(
        "values" = c(out[which(rowSums(abs(out)) > 0),]),
        "log10lambda" = rep(log10(fit$full.fit$lambda), each = sum(rowSums(abs(out)) != 0)),
        "predictor" = rep(1:sum(rowSums(abs(out)) != 0), length(fit$full.fit$lambda))
      )
      values <- dat$values
      log10lambda <- dat$log10lambda
      predictor <- dat$predictor
      
      p1 <- ggplot(dat, aes(x=log10lambda, y = values, color=as.factor(predictor))) +
        theme_bw() + theme(legend.position="") + ylab(expression(hat(beta)[j])) +
        geom_line()  +  ggtitle("Elastic net trace plot") +
        xlab(expression(log[10](lambda))) + theme(plot.title = element_text(hjust = 0.5))
    }

  } else {

  groups <- fit$full.fit$groups
  out <- matrix(0, nrow=length(unique(groups)), ncol=length(fit$full.fit$lambda))

  for(k in 1:length(fit$full.fit$lambda)) {
    for(j in unique(groups)){
      out[j,k] <- sqrt(sum(beta[groups==j,k]^2))
    }
  }

  if (inherits(fit, "penAFT.cv")) {
    dat2 <- data.frame(
      "log10lambda" = log10(fit$full.fit$lambda),
      "linPred" = fit$cv.err.linPred,
      "ObjErr" = colMeans(fit$cv.err.obj),
      "sesObjErr" = apply(fit$cv.err.obj, 2, sd)/sqrt(dim(fit$cv.err.obj)[1])
    )
  
   # A few constants
    t1 <- "dodgerblue3"
    t2 <- "black"
    if(is.null(groupNames)){
      dat <- data.frame(
        "values" = c(t(out)),
        "Groups" = as.factor(rep(unique(groups), each = length(fit$full.fit$lambda))),
        "log10lambda" = rep(log10(fit$full.fit$lambda), length(unique(groups))))
    } else {
       dat <- data.frame(
        "values" = c(t(out)),
        "Groups" = as.factor(rep(groupNames, each = length(fit$full.fit$lambda))),
        "log10lambda" = rep(log10(fit$full.fit$lambda), length(unique(groups))))
    }
  
    Groups <- dat$Groups
    values <- dat$values
    log10lambda <- dat$log10lambda
  
    p1 <- ggplot(dat, aes(x=log10lambda, y = values, color=Groups)) +
      theme_bw() + ylab(expression("||"*hat(beta)[G[g]]*"||"[2])) +
      geom_line() + geom_vline(xintercept = log10(fit$full.fit$lambda)[which.min(fit$cv.err.linPred)], linetype="dotted",
                  color = t1) + geom_vline(xintercept = log10(fit$full.fit$lambda)[min(which(dat2$ObjErr <= min(dat2$ObjErr + dat2$sesObjErr)))], linetype="dotted",
                  color = t2) +  ggtitle("Group lasso trace plot") +
                  xlab(expression(log[10](lambda))) + theme(plot.title = element_text(hjust = 0.5))
  } else {

    # A few constants
    t1 <- "dodgerblue3"
    t2 <- "black"
    if(is.null(groupNames)){
      dat <- data.frame(
        "values" = c(t(out)),
        "Groups" = as.factor(rep(unique(groups), each = length(fit$full.fit$lambda))),
        "log10lambda" = rep(log10(fit$full.fit$lambda), length(unique(groups))))
    } else {
      dat <- data.frame(
        "values" = c(t(out)),
        "Groups" = as.factor(rep(groupNames, each = length(fit$full.fit$lambda))),
        "log10lambda" = rep(log10(fit$full.fit$lambda), length(unique(groups))))
    }
    
    Groups <- dat$Groups
    values <- dat$values
    log10lambda <- dat$log10lambda
    
    p1 <- ggplot(dat, aes(x=log10lambda, y = values, color=Groups)) +
      theme_bw() + ylab(expression("||"*hat(beta)[G[g]]*"||"[2])) +
      geom_line()  +  ggtitle("Group lasso trace plot") +
      xlab(expression(log[10](lambda))) + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  return(p1)
}





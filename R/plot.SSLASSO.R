#' Plot coefficients from a "SSLASSO" object
#'
#' Produces a coefficient profile plot of the SSLASSO solution path.
#'
#' @param x An object of class "SSLASSO", usually, a result of a call to \code{SSLASSO}.
#' @param ... Other arguments to be passed to \code{plot}.
#' 
#' @return No return value, called for side effects (produces a plot of the coefficient paths).
#'  
#' @importFrom graphics matplot matpoints title par text
#' @exportS3Method plot SSLASSO
#'
plot.SSLASSO<-function(x, ...){

  # Save and restore graphical parameters on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
	betas <-t(x$beta)
        
     	v0s <-x$lambda0

        select <-apply(betas,2,function(x){as.numeric(x!=0)})

        matplot(v0s,betas,xlab=expression(lambda[0]),ylab=expression(hat(beta)),lwd=1,col="grey",lty=2,type="l")

        betas[select==F]=NA

        matpoints(v0s,betas*select,xlab=expression(lambda[0]),ylab=expression(hat(beta)),lwd=1,col=4,lty=2,pch=19)

        matpoints(v0s,betas*(1-select),xlab=expression(lambda[0]),ylab=expression(hat(beta)),lwd=1,col=2,lty=2,pch=19)

        title("Spike-and-Slab LASSO")

        par(xpd=T)

        labels=paste("X",1:ncol(betas),sep="")

        labels[select[length(v0s),]==0]<-""

        text(max(v0s)*(1.1),betas[length(v0s),],labels=labels,col=4)
    	


}



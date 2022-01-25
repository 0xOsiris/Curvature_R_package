### Differential Geometry - Curvature of Parametric Curve R Package

### Installation
```markdown
> library(devtools)
> install_github("LeytonTaylor/Curvature_R_package");
> library(CurvatureCalc);
> curvature(c(expression(1+.5*cos(t),1+.1*sin(t),sin(t))), 1)
```

### Usage
```markdown
@param alpha A vector of 3 expressions defining the parametric curve
@param t Evaluation parameter
@return Scalar value Curvature of alpha at time t
@details alpha must be a vector of 3 parametric expressions with param 't'
@examples 
kAlpha <- curvature(c(expression(1+.5*cos(t),1+.1*sin(t),-1)), 1)
dataSet= matrix(c(expression(1+.5*cos(t),1+.1*sin(t),sin(t))), nrow=1,ncol=3,byrow=TRUE)
alphaN=dataSet[1,]
kAlpha<-curvature(alphaN,1)
```

### Components

## Curvature Calculation

```markdown

curvature <- function(alpha,t){
  
  alphaPrimeMatrix <- d_util$alphaDerivativeMat(alpha,t)
  alphaP <- c(eval(alphaPrimeMatrix[2,1]),eval(alphaPrimeMatrix[2,2]),eval(alphaPrimeMatrix[2,3]))
  alphaDP <- c(eval(alphaPrimeMatrix[3,1]),eval(alphaPrimeMatrix[3,2]),eval(alphaPrimeMatrix[3,3]))
  normCross <-xprod_util$xprod(alphaP,alphaDP)
  normAD <-norm(as.matrix(normCross),"F")
  normAP <-norm(as.matrix(alphaP),"F")
  kAlpha = normAD /(normAP^3)
  
  
  return(kAlpha)
  
}


```
## Cross product calculation
```markdown

#' @export
xprod_util = new.env()
xprod_util$xprod <- function(...) {
  args <- list(...)

  if (length(args) == 0) {
    stop("No data supplied")
  }
  len <- unique(sapply(args, FUN=length))
  if (length(len) > 1) {
    stop("All vectors must be the same length")
  }
  if (len != length(args) + 1) {
    stop("Must supply N-1 vectors of length N")
  }

  m <- do.call(rbind, args)
  sapply(seq(len),
         FUN=function(i) {
           det(m[,-i,drop=FALSE]) * (-1)^(i+1)
         })
}


```

## Derivative Matrix
```markdown
#' @export
d_util = new.env()

d_util$alphaDerivativeMat <- function(alpha, t){
  
  alphaP <- c(D(alpha[1],'t'), D(alpha[2],'t'),D(alpha[3],'t'))
  alphaDP <- c(D(D(alpha[1],'t'),'t'),D(D(alpha[2],'t'),'t'),D(D(alpha[3],'t'),'t'))
  
  alphaPrimeMatrix <-matrix(
    c(alpha,alphaP,alphaDP),
    nrow=3,
    ncol=3,
    byrow=TRUE
    
  )
  
  dimnames(alphaPrimeMatrix)=list(c("alpha","alphaP","alphaDP"),c("X","Y","Z"))
  
  return(alphaPrimeMatrix)
  
}

```


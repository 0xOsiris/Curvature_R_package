# Differential Geometry - R Utility to Calculate Curvature of a Parametric Curve

# Installation Instructions
#### > library(devtools)
#### > install_github("LeytonTaylor/Curvature_R_package");
#### > library(CurvatureCalc);
#### > curvature(c(expression(1+.5*cos(t),1+.1*sin(t),sin(t))), 1)
# Usage:
#### @param alpha A vector of expressions defining the curve
#### @param t Evaluation parameter
#### @return Scalar value Curvature of alpha at time t
#### @details alpha must be a vector of parametric expressions with param 't'
#### @examples 
#### kAlpha <- curvature(c(expression(1+.5*cos(t),1+.1*sin(t),-1)), 1)
#### dataSet= matrix(c(expression(1+.5*cos(t),1+.1*sin(t),sin(t))), nrow=1,ncol=3,byrow=TRUE)
#### alphaN=dataSet[1,]
#### kAlpha<-curvature(alphaN,1)

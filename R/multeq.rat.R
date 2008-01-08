`multeq.rat` <-
function(data,grp,resp=NULL,margin.lo=NULL,margin.up=NULL,method="single.step",FWER=0.05) {

if (length(grp) > 1) {
  stop("Specify only one grouping variable")
}
tr.names <- levels(data[,grp])
if (length(tr.names) > 2) {
  stop("Grouping factor must have exactly 2 levels")
}
comp.name <- paste(tr.names[1], tr.names[2], sep = "-")
Resp.X <- subset(data, data[,grp]==tr.names[1])
Resp.Y <- subset(data, data[,grp]==tr.names[2])
if (is.null(resp)) {
  n.ep <- length(names(data))-1                                          # number of endpoints
  Resp.X <- Resp.X[,-which(names(data)%in%grp)]
  Resp.Y <- Resp.Y[,-which(names(data)%in%grp)]
} else {
  n.ep <- length(resp)                                                   # number of endpoints
  Resp.X <- Resp.X[resp]
  Resp.Y <- Resp.Y[resp]
}
if (is.numeric(margin.lo) & length(margin.lo) != n.ep) {
  stop("Length of margin.lo is not equal to the number of response variables")
}
if (is.numeric(margin.up) & length(margin.up) != n.ep) {
  stop("Length of margin.up is not equal to the number of response variables")
}
method <- match.arg(method, choices = c("single.step", "step.up"))

X.n <- nrow(Resp.X); Y.n <- nrow(Resp.Y)                                 # sample sizes for X and Y
test.stat <- numeric(n.ep)

X.mean <- mean(Resp.X); Y.mean <- mean(Resp.Y)                           # mean vectors for X and Y
estimate <- X.mean/Y.mean
cov.mat=((X.n-1)*cov(Resp.X)+(Y.n-1)*cov(Resp.Y))/(X.n+Y.n-2)            # common estimated covariance matrix of the data
                                                                         # just the variances needed
degr.fr <- X.n+Y.n-2                                                     # df

if (is.numeric(margin.lo)) {
  T.up <- (X.mean-margin.lo*Y.mean)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + (margin.lo^2)/Y.n )  )
}
if (is.numeric(margin.up)) {
  T.do <- (X.mean-margin.up*Y.mean)/(  sqrt(diag(cov.mat))*sqrt( 1/X.n + (margin.up^2)/Y.n )  )
}

if ((is.numeric(margin.lo)) & is.numeric(margin.up)) {
  for (i in 1:n.ep) test.stat[i] <- min(T.up[i],-T.do[i])
}
if ((is.numeric(margin.lo)) & is.numeric(margin.up)==FALSE) {
  test.stat <- T.up
}
if ((is.numeric(margin.lo)==FALSE) & is.numeric(margin.up)) {
  test.stat <- -T.do
}

if (method == "step.up")
{
  p.value <- numeric(n.ep)
  g <- numeric(n.ep)                                                     # see Dilba et al. 2004 p. 446
  lower <- numeric(n.ep); upper <- numeric(n.ep)
  for (i in 1:n.ep){
    pos <- which(p.value<FWER)
    if (length(pos)==n.ep+1-i){
      p.value[p.value<FWER]=pt(q=test.stat[pos],df=degr.fr,lower.tail=FALSE)*i
      g[pos]=(qt(p=1-FWER/i,df=degr.fr))^2*diag(cov.mat)[pos]/( Y.n*(Y.mean[pos])^2 )
      if (is.numeric(margin.lo)) {
        lower[pos]=( estimate[pos]-sqrt(g[pos])*sqrt( (estimate[pos])^2+(1-g[pos])*Y.n/X.n ) )/(1-g[pos])
      } else { lower[pos]=-Inf }
      if (is.numeric(margin.up)) {
        upper[pos]=( estimate[pos]+sqrt(g[pos])*sqrt( (estimate[pos])^2+(1-g[pos])*Y.n/X.n ) )/(1-g[pos])
      } else { upper[pos]=Inf }
    }
  p.value[p.value>1]=1
  }
}

if (method == "single.step")
{
  p.value=pt(q=test.stat,df=degr.fr,lower.tail=FALSE)*n.ep; p.value[p.value>1]=1
  g=(qt(p=1-FWER/n.ep,df=degr.fr))^2*diag(cov.mat)/( Y.n*(Y.mean)^2 )
  if (is.numeric(margin.lo)) {
    lower=( estimate-sqrt(g)*sqrt( (estimate)^2+(1-g)*Y.n/X.n ) )/(1-g)
  } else { lower=rep(-Inf,n.ep) }
  if (is.numeric(margin.up)) {
    upper=( estimate+sqrt(g)*sqrt( (estimate)^2+(1-g)*Y.n/X.n ) )/(1-g)
  } else { upper=rep(Inf,n.ep) }
}

if (any(g>1)) lower <- upper <- "NSD"

value <- list(comp.name=comp.name,estimate=estimate,degr.fr=degr.fr,test.stat=test.stat,
              p.value=p.value,lower=lower,upper=upper,margin.lo=margin.lo,margin.up=margin.up,
              method=method,FWER=FWER)

names(value$test.stat) <- resp
names(value$p.value) <- resp
names(value$lower) <- names(value$upper) <- resp
if (is.numeric(margin.lo)) { names(value$margin.lo) <- resp }
if (is.numeric(margin.up)) { names(value$margin.up) <- resp }
class(value) <- "multeq.rat"

return(value)

}


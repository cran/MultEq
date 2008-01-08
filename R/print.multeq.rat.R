`print.multeq.rat` <-
function(x,digits=4,...) {

cat("", "\n")
cat("Alternative hypotheses: ratios ")
if (is.numeric(x$margin.lo)) cat("larger than", x$margin.lo)
if (is.numeric(x$margin.lo) & is.numeric(x$margin.up)) cat(" and ")
if (is.numeric(x$margin.up)) cat ("smaller than", x$margin.up)
cat("", "\n")
cat("Method:", x$method,
    "\n")
if (is.numeric(x$lower) & is.numeric(x$upper)) {
  out <- cbind(x$estimate, x$lower, x$upper, x$p.value)
  colnames(out) <- c("estimate", "lower", "upper", "p.value")
} else {
  out <- cbind(x$estimate, x$p.value)
  colnames(out) <- c("estimate", "p.value")
  cat("                                      ", "\n")
  cat("   The mean in the denominator is not significantly different from zero. ", 
      "\n")
  cat("                                      ", "\n")
}
cat("", "\n")
print(out, digits=digits)
cat("", "\n")
invisible(x)

}


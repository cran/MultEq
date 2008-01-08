`print.multeq.diff` <-
function(x,digits=4,...) {

cat("", "\n")
cat("Alternative hypotheses: differences ")
if (is.numeric(x$margin.lo)) cat("larger than", x$margin.lo)
if (is.numeric(x$margin.lo) & is.numeric(x$margin.up)) cat(" and ")
if (is.numeric(x$margin.up)) cat ("smaller than", x$margin.up)
cat("", "\n")
cat("Method:", x$method,
    "\n")
out <- cbind(x$estimate, x$lower, x$upper, x$p.value)
colnames(out) <- c("estimate", "lower", "upper", "p.value")
cat("", "\n")
print(out, digits=digits)
cat("", "\n")
invisible(x)

}


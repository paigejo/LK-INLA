
library(latex2exp)

phi1 = function(x) {
  2 - 2*x
}
phi2 = function(x) {
  2*x
}
f = function(x, c1=1, c2=1) {
  c1*phi1(x) + c2*phi2(x)
}
varphi1 = function(x, varc1=1) {
  varc1*phi1(x)^2
}
varphi2 = function(x, varc2=1) {
  varc2*phi2(x)^2
}
covTerm = function(x, varc1=1, varc2=varc1, corc1c2=0) {
  2*corc1c2*phi1(x)*phi2(x)
}
varfx = function(x, varc1=1, varc2=varc1, corc1c2=0) {
  varc1*phi1(x)^2 + varc2*phi2(x)^2 + 2*corc1c2*phi1(x)*phi2(x)
}

pdf(file="meshVariance.pdf", height=6, width=6)
par(mfrow=c(2, 2))
xs = seq(0, 1, l=100)
plot(xs, f(xs), type="l", xlab="x", main=TeX("f(x)"), ylim=c(0,4), ylab="")
lines(xs, phi1(xs), col="blue", lty=2)
lines(xs, phi2(xs), col="blue", lty=2)
legend("top", c(TeX("f(x)"), TeX("$c_1 \\phi_1(x)$"), TeX("$c_2 \\phi_2(x)$")), cex=.7, 
       lty=c(1, 2, 2), col=c("black", "blue", "blue"))

cols = rainbow(3)
plot(xs, varfx(xs), type="l", xlab="x", main=TeX("Var(f(x)), Cor($c_1$, $c_2$)=0"), ylim=c(0,7), ylab="")
lines(xs, varphi1(xs), col=cols[1], lty=2)
lines(xs, varphi2(xs), col=cols[2], lty=2)
lines(xs, covTerm(xs), col=cols[3], lty=2)
legend("top", c(TeX("Var(f(x))"), TeX("$c_1 \\phi_1(x)$"), TeX("$c_2 \\phi_2(x)$"), TeX("$2 Cov(c_1, c_2) \\phi_1(x) \\phi_2(x)$")), cex=.7, 
       lty=c(1, 2, 2, 2), col=c("black", cols))

plot(xs, varfx(xs, corc1c2=.5), type="l", xlab="x", main=TeX("Var(f(x)), Cor($c_1$, $c_2)=0.5$"), ylim=c(0,7), ylab="")
lines(xs, varphi1(xs), col=cols[1], lty=2)
lines(xs, varphi2(xs), col=cols[2], lty=2)
lines(xs, covTerm(xs, corc1c2=.5), col=cols[3], lty=2)
legend("top", c(TeX("Var(f(x))"), TeX("$c_1 \\phi_1(x)$"), TeX("$c_2 \\phi_2(x)$"), TeX("$2 Cov(c_1, c_2) \\phi_1(x) \\phi_2(x)$")), cex=.7, 
       lty=c(1, 2, 2, 2), col=c("black", cols))

plot(xs, varfx(xs, corc1c2=1), type="l", xlab="x", main=TeX("Var(f(x)), Cor($c_1$, $c_2$)=1"), ylim=c(0,7), ylab="")
lines(xs, varphi1(xs), col=cols[1], lty=2)
lines(xs, varphi2(xs), col=cols[2], lty=2)
lines(xs, covTerm(xs, corc1c2=1), col=cols[3], lty=2)
legend("top", c(TeX("Var(f(x))"), TeX("$c_1 \\phi_1(x)$"), TeX("$c_2 \\phi_2(x)$"), TeX("$2 Cov(c_1, c_2) \\phi_1(x) \\phi_2(x)$")), cex=.7, 
       lty=c(1, 2, 2, 2), col=c("black", cols))
dev.off()



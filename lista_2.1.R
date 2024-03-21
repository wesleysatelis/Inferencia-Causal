# Consideremos uma amostra de tamanho 100 com 60 tratamentos e 40 placebos em que os resultados potenciais são gerados com efeito individual constante.
# A Scientific Table é então

n <- 100
n1 <- 60
n0 <- 40
y0 <- sort(rexp(n), decreasing = TRUE)
tau <- 1
y1 <- y0 + tau

# Escolhemos 1 amostra e calculamos: τ^, V^, intervalo de confiança e V(τ^).
z <- sample(c(rep(1, n1), rep(0, n0)))
tau_hat <- mean(y1[z == 1]) - mean(y0[z == 0])
tau_hat

s_1 <- sqrt(sum(y1[z==1]-mean(y1[z==1]))^2/(n-1))
s_0 <- sqrt(sum(y0[z==0]-mean(y0[z==0]))^2/(n-1))

V_hat <- var(y1[z == 1]) / n1 + var(y0[z == 0]) / n0

intervalo_confianca <- c(tau_hat - 1.96 * sqrt(V_hat), tau_hat + 1.96 * sqrt(V_hat))
intervalo_confianca

var_hat_tau <- var(y1) / n1 + var(y0) / n0 - var(y1 - y0) / n
var_hat_tau

# -------------------------------------------------------------------------

tau_hat_p <- c()
V_hat_p <- c()
lim_sup <- c()
lim_inf <- c()
est <- c()
for (i in 1:10^4) {
  z_permut <- sample(z)
  tau_hat_p[i] <- mean(y1[z_permut == 1]) - mean(y0[z_permut == 0])
  V_hat_p[i] <- var(y1[z_permut == 1]) / n1 + var(y0[z_permut == 0]) / n0
  lim_sup[i] <- tau_hat_p[i] + 1.96*sqrt(V_hat_p[i])
  lim_inf[i] <- tau_hat_p[i] - 1.96*sqrt(V_hat_p[i])
  est[i] <- (tau_hat_p[i] - tau) / sqrt(V_hat_p[i])
}
mean(V_hat_p)

cobertura <- c()
for (i in 1:10^4) {
  cobertura[i] <- ifelse(lim_inf[i] < tau && tau < lim_sup[i], 1, 0)
}
mean(cobertura)
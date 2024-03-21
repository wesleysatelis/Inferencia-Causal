library(Matching)
data(lalonde)

# QUESTÃO 3 ---------------------------------------------------------------

# t studentizada ----------------------------------------------------------
z <- lalonde$treat
y <- lalonde$re78

n <- nrow(lalonde)
n1 <- sum(z) # quantos receberam tratamento
n0 <- n - n1 # quantos não receberam tratamento
tau <- mean(y[z == 1]) - mean(y[z == 0]) # efeito causal médio
s2 = var(y)

media_y_1 <- mean(y[z == 1])
s2_1 <- sum((y[z == 1] - media_y_1)^2)/(n1 - 1)

media_y_0 <- mean(y[z == 0])
s2_0 <- sum((y[z == 0] - media_y_0)^2)/(n0 - 1)

est_teste_frt <- tau / sqrt(s2_1/n1 + s2_0/n0)

# a t studentizada é assintoticamente N(0, 1), então
pvalor_frt_clt <- 2 - 2*pnorm(est_teste_frt) # Bilateral

# FRT Monte Carlo
mc <- 10^5
est_teste_frt_mc <- rep(0, mc)
for (i in 1:mc) {
  zpermut <- sample(z) # permutando z
  media_y_1 <- mean(y[zpermut == 1])
  s2_1 <- sum((y[zpermut == 1] - media_y_1)^2)/(n1 - 1)
  
  media_y_0 <- mean(y[zpermut == 0])
  s2_0 <- sum((y[zpermut == 0] - media_y_0)^2)/(n0 - 1)
  
  tau_permut <- mean(y[zpermut == 1]) - mean(y[zpermut == 0]) # efeito causal medio da permutação
  est_teste_frt_mc[i] <- tau_permut/sqrt(s2_1/n1 + s2_0/n0) # estatistica do teste da permutação
}

pvalor_frt_mc <- mean(abs(est_teste_frt_mc) >= abs(est_teste_frt)); pvalor_frt_mc
# H0: não existe efeito causal
# se p-valor < 0.05, rejeita H0 a um nivel de significância de 95%.

# FRT clássico teste t para dif medias 
# essa é a versão assintótica?
est_test_frt_classico_t <- t.test(y[z == 1], y[z == 0], var.equal = FALSE)$statistic; est_test_frt_classico_t
pvalor_frt_classico_t <- t.test(y[z == 1], y[z == 0], var.equal = FALSE)$p.value; pvalor_frt_classico_t


# wilcoxon ----------------------------------------------------------------
z <- lalonde$treat
y <- lalonde$re78

W_obs <- sum(rank(y)[z == 1])

# FRT Monte Carlo
mc <- 10^5
W_mc <- rep(0, mc)
for (i in 1:mc) {
  zpermut <- sample(z) # permutando z
  W_mc[i] <- sum(rank(y)[zpermut == 1]) # estatistica do teste da permutação
}

pvalor_W <- mean(W_mc >= W_obs); pvalor_W

wilcox.test(y[z == 1], y[z==0], exact = FALSE)$p.value
wilcox.test(y[z == 1], y[z==0], exact = TRUE)$p.value


# QUESTÃO 4 ---------------------------------------------------------------

teste <- function(mc){
  z <- rbinom(100, 1, 0.7)
  y <- rnorm(100, 0.5, 5)
  
  s2 <- sum((y-mean(y))^2)/(n-1)
  estat <- (mean(y[z==1]) - mean(y[z==0]))/(n*s2/(n1*n0))
  
  # mc <- 10^5
  estat_mc <- numeric(mc)
  for(i in 1:mc){
    zpermut <- sample(z)
    estat_mc[i] <- (mean(y[zpermut==1]) - mean(y[zpermut==0]))/(n*s2/(n1*n0))
  }
  
  # bilateral
  p_valor1 <- (1+sum(estat_mc >= estat))/(1+mc)
  p_valor2 <- sum(estat_mc >= estat)/mc
  
  return(tibble(R=mc, p_chapeu=p_valor1, p_til=p_valor2))
}

teste(100)

permutation = function(n, n1){
  M = choose(n, n1)
  treat.index = combn(n, n1)
  Z_matriz = matrix(0, n, M)
  for(m in 1:M){
    treat = treat.index[, m]
    Z_matriz[treat, m] = 1
  }
  return(Z_matriz)
}

Z_matriz = permutation(n, n1)

## p_FRT
estat_comb = c()
for(i in 1:ncol(Z_matriz)){
  zpermut = Z_matriz[, i]
  estat_comb[i] = (mean(y[zpermut==1]) - mean(y[zpermut==0]))/(n*s2/(n1*n0))
}
mean(abs(estat_comb) >= abs(estat))







## Exercicio 04 ----------------------------------------------------------------
set.seed(123)
n = 20
y0 = rnorm(n)
tau = rnorm(n, -0.5)
y1 = y0 + tau
z = rbinom(n, 1, 0.5)
y = z*y1 + (1-z)*y0
n1 = sum(z)
n0 = n - n1

## efeito causal medio
mean(y[z == 1]) - mean(y[z == 0])

s2 =  sum((y-mean(y))^2)/(n-1)
estat = (mean(y[z==1]) - mean(y[z==0]))/(n*s2/(n1*n0))

permutation = function(n, n1){
  M = choose(n, n1)
  treat.index = combn(n, n1)
  Z_matriz = matrix(0, n, M)
  for(m in 1:M){
    treat = treat.index[, m]
    Z_matriz[treat, m] = 1
  }
  return(Z_matriz)
}

Z_matriz = permutation(n, n1)

## p_FRT
estat_comb = c()
for(i in 1:ncol(Z_matriz)){
  zpermut = Z_matriz[, i]
  estat_comb[i] = (mean(y[zpermut==1]) - mean(y[zpermut==0]))/(n*s2/(n1*n0))
}
mean(abs(estat_comb) >= abs(estat))


## p_FRT chapeu (mc)
R = 10^5
estat_mc = c()
for(i in 1:R){
  zpermut = sample(z)
  estat_mc[i] = (mean(y[zpermut==1]) - mean(y[zpermut==0]))/(n*s2/(n1*n0))
}
mean(abs(estat_mc) >= abs(estat))

## p_FRT tio (finite-sample valid Monte Carlo approximation)
(1+sum(abs(estat_mc) >= abs(estat)))/(R+1)



data(lalonde)

X <- as.matrix(lalonde[, c('age', 'educ', 'black', 'hisp', 'married', 'nodegr', 're74', 're75', 'u74', 'u75')])
Z <- lalonde$treat
Y <- lalonde$re78

# padroniza covariaveis
X_media <- colMeans(X)
X_sd <- apply(X, 2, sd)
X_std <- scale(X, center = X_media, scale = X_sd)

# tratado e controle
X_trat <- X_std[Tr == 1, ]
X_contr <- X_std[Tr == 0, ]

# distancia de Mahalanobis
mahalanobis <- function(x, y, cov_mat) {
  diff <- x - y
  sqrt(t(diff) %*% solve(cov_mat) %*% diff)
}

# matriz de cov das covariaveis
cov_mat <- cov(X_contr)

# mais proximos da unidade tratada
proximos <- integer(nrow(X_trat))
for (i in 1:nrow(X_trat)) {
  dist <- numeric(nrow(X_contr))
  for (j in 1:nrow(X_contr)) {
    dist[j] <- mahalanobis(X_trat[i, ], X_contr[j, ], cov_mat)
  }
  proximos[i] <- which.min(dist)
}

# match exato do outcome controle
match_contr_outc <- Y[Tr == 0][proximos]

# efeito causal
tau <- mean(Y[Tr == 1] - match_contr_outc)
# print(paste("Tau sem correção de viés:", tau))

# correção de viés
vies_correc <- numeric(nrow(X_trat))

for (i in 1:nrow(X_trat)) {
  # i <- 3
  match_control <- X_contr[proximos[i], ]
  dist <- numeric(nrow(X_contr))
  for (j in 1:nrow(X_contr)) {
    dist[j] <- mahalanobis(match_control, X_contr[j, ], cov_mat)
  }
  dis_ordenadas <- order(dist)
  proximos_contr <- X_contr[dis_ordenadas[1:nrow(X_trat)], ]
  vies_correc[i] <- sum(X_trat[i, ] - colMeans(proximos_contr)) / ncol(X_trat)
}

tau_adj <- tau + mean(vies_correc)
print(paste("Tau com correção de viés:", tau_adj))


# com o pacote ------------------------------------------------------------
match <- Match(Y = Y, Tr = Z, X = X, BiasAdjust = TRUE)
summary(match)


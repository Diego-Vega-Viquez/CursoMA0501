library(tidyverse)

neville <- function(nodos, valores, x) {
  stopifnot(is.numeric(nodos),
            is.numeric(valores),
            length(nodos) == length(valores))
  n <- length(nodos)
  Q <- matrix(NA_real_, nrow = n, ncol = n)
  
  # Columna inicial con valores de Y
  Q[, 1] <-  valores
  
  # Construccion de la tabla de Neville
  for (i in 2:n) {
    for (j in 2:i) {
      numerador <- ((x - nodos[i - j + 1]) * Q[i, j - 1] - (x - nodos[i]) * Q[i -
                                                                                1, j - 1])
      denominador <- nodos[i] - nodos[i - j + 1]
      Q[i, j] <-  numerador / denominador
    }
  }
  return(list(valor = Q[n, n], tabla = Q))
}


graficar.polinomio <- function(nodos,
                               a,
                               b,
                               metodo,
                               f = NULL,
                               valores = NULL,
                               df = NULL,
                               derivadas.clamped = NULL,
                               f.referencia = NULL) {
  stopifnot(is.numeric(nodos), length(nodos) >= 2)
  
  if (is.null(f) == is.null(valores)) {
    stop("Debe proveer exactamente uno: 'f' (función) o 'valores' (numérico).")
  }
  
  if (!is.null(f)) {
    stopifnot(is.function(f))
    valores_nodos <- f(nodos)
  } else {
    stopifnot(is.numeric(valores), length(valores) == length(nodos))
    valores_nodos <- valores
  }
  
  derivadas_nodos <- NULL
  if (!is.null(df)) {
    if (is.function(df)) {
      derivadas_nodos <- df(nodos)
    } else if (is.numeric(df)) {
      stopifnot(length(df) == length(nodos))
      derivadas_nodos <- df
    } else {
      stop("`df` debe ser función o vector numérico de derivadas en los nodos.")
    }
  }
  
  if (!is.null(derivadas_nodos)) {
    H <- function(x)
      vapply(
        x,
        function(xx)
          metodo(nodos, valores_nodos, derivadas_nodos, xx)$valor,
        numeric(1)
      )
  } else {
    if (is.null(derivadas.clamped)) {
      H <- function(x)
        vapply(
          x,
          function(xx)
            metodo(nodos, valores_nodos, xx)$valor,
          numeric(1)
        )
    } else {
      H <- function(x)
        vapply(
          x,
          function(xx)
            spline.sujeto(nodos, valores_nodos, derivadas.clamped, xx)$valor,
          numeric(1)
        )
    }
  }
  
  xi <- seq(a, b, length.out = 400)
  
  fx.col <- if (!is.null(f)) f(xi) else NA_real_
  
  fx.ref.col <- if (!is.null(f.referencia)) {
    stopifnot(is.function(f.referencia))
    f.referencia(xi)
  } else {
    NA_real_
  }
  
  df_plot <- data.frame(
    x      = xi,
    Hx     = H(xi),
    fx     = fx.col,
    fx.ref = fx.ref.col
  )
  
  df_nodos <- data.frame(x = nodos, y = valores_nodos)
  
  p <- ggplot(df_plot, aes(x = x))
  
  if (!is.null(f)) {
    p <- p +
      geom_line(
        aes(y = fx,
            color = "Original",
            linetype = "Original"),
        linewidth = 1
      )
  }
  
  p <- p +
    geom_line(
      aes(y = Hx,
          color = "Interpolación",
          linetype = "Interpolación"),
      linewidth = 1
    )
  
  if (!is.null(f.referencia)) {
    p <- p +
      geom_line(
        aes(y = fx.ref,
            color = "Referencia",
            linetype = "Referencia"),
        linewidth = 1
      )
  }
  
  p <- p +
    geom_point(
      data = df_nodos,
      aes(x = x, y = y),
      shape = 21,
      size = 3,
      fill = "white"
    ) +
    scale_color_manual(
      values = c(
        "Original"        = "blue",
        "Interpolación"   = "red",
        "Referencia"      = "darkgreen"
      ),
      name = "Serie"
    ) +
    scale_linetype_manual(
      values = c(
        "Original"        = "solid",
        "Interpolación"   = "dashed",
        "Referencia"      = "solid"
      ),
      guide = "none"
    ) +
    labs(
      title = paste0("Interpolación por ", deparse(substitute(metodo))),
      y = "Valor"
    ) +
    theme_minimal(base_size = 14)
  
  p
}

#############################################

metodo.euler <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    W[i] <- W[i-1] + h*f(T[i-1],W[i-1])
    T[i] <- t + (i - 1)*h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

prueba <- metodo.euler(1, 2, 10, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos <- prueba$t
valores <- prueba$w

graficar.polinomio(nodos, 1, 2, neville, valores = valores, f.referencia = function(x) x^2*(exp(x)-exp(1)))

metodo.euler.predictor.corrector <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    W[i] <- W[i-1] + h*f(T[i-1],W[i-1]) #predictor
    T[i] <- T[i - 1] + h
    W[i] <- W[i-1] + (h/2)*(f(T[i-1], W[i-1]) + f(T[i], W[i])) #corrector
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

prueba.corrector <- metodo.euler.predictor.corrector(1, 2, 10, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.corr <- prueba.corrector$t
valores.corr <- prueba.corrector$w

graficar.polinomio(nodos.corr, 1, 2, neville, valores = valores.corr, f.referencia = function(x) x^2*(exp(x)-exp(1)))

runge.kutta.punto.medio <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    W[i] <- W[i-1] + h*f(T[i-1] + (h/2),W[i-1] + (h/2)*f(T[i-1], W[i-1])) #predictor
    T[i] <- T[i - 1] + h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

pba.rkpm <- runge.kutta.punto.medio(1, 2, 10, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.rkpm <- pba.rkpm$t
valores.rkpm <- pba.rkpm$w

graficar.polinomio(nodos.rkpm, 1, 2, neville, valores = valores.rkpm, f.referencia = function(x) x^2*(exp(x)-exp(1)))

runge.kutta.cuarto.orden <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    k1 <- h*f(T[i-1], W[i-1])
    k2 <- h*f(T[i-1] + h/2, W[i-1] + k1/2)
    k3 <- h*f(T[i-1] + h/2, W[i-1] + k2/2)
    k4 <- h*f(T[i-1] + h, W[i-1] + k3)
    
    W[i] <- W[i-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    T[i] <- T[i - 1] + h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

pba.rk4o <- runge.kutta.cuarto.orden(1, 2, 3, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.rk4o <- pba.rk4o$t
valores.rk4o <- pba.rk4o$w

graficar.polinomio(nodos.rk4o, 1, 2, neville, valores = valores.rk4o, f.referencia = function(x) x^2*(exp(x)-exp(1)))


euler.modificado <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    W[i] <- W[i-1] + h*f(T[i-1] + h/2, W[i-1] + (h/2)*f(T[i-1], W[i-1]))
    T[i] <- T[i - 1] + h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

pba.e.mod <- euler.modificado(1, 2, 3, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.e.mod <- pba.e.mod$t
valores.e.mod <- pba.e.mod$w

graficar.polinomio(nodos.e.mod, 1, 2, neville, valores = valores.e.mod, f.referencia = function(x) x^2*(exp(x)-exp(1)))

tercer.orden.huen <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    k1 <- f(T[i-1], W[i-1])
    k2 <- f(T[i-1] + h/3, W[i-1] + k1*(h/3))
    k3 <- f(T[i-1] + (2/3)*h, W[i-1] + (2*h/3)*k2)

    W[i] <- W[i-1] + (h/4)*(k1 + 3*k3)
    T[i] <- T[i - 1] + h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

pba.huen <- euler.modificado(1, 2, 3, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.huen <- pba.huen$t
valores.huen <- pba.huen$w

graficar.polinomio(nodos.huen, 1, 2, neville, valores = valores.huen, f.referencia = function(x) x^2*(exp(x)-exp(1)))

tercer.orden.runge.kutta <- function(a, b, N, alfa, f){
  h <- (b - a) / N
  t <- a
  w <- alfa
  T <- numeric(N + 1)
  W <- numeric(N + 1)
  T[1] <- t
  W[1] <- w
  for (i in 2:(N+1)) {
    k1 <- f(T[i-1], W[i-1])
    k2 <- f(T[i-1] + h/2, W[i-1] + k1*(h/2))
    k3 <- f(T[i-1] + h, W[i-1] - h*k1 + 2*h*k2)
    
    W[i] <- W[i-1] + (h/6)*(k1 + 4*k2 + k3)
    T[i] <- T[i - 1] + h
  }
  
  tabla <- data.frame(t = T, w = W)
  
  return(list(t = T, w = W, tabla = tabla))
}

pba.runge.kutta.3 <- euler.modificado(1, 2, 3, 0, function(t, w) (2/t)*w + t^2*exp(t))
nodos.runge.kutta.3 <- pba.runge.kutta.3$t
valores.runge.kutta.3 <- pba.runge.kutta.3$w

graficar.polinomio(nodos.runge.kutta.3, 1, 2, neville, valores = valores.runge.kutta.3, f.referencia = function(x) x^2*(exp(x)-exp(1)))


exp_matriz <- function(A, t) {
  eig <- eigen(A)
  C <- eig$vectors
  D <- diag(exp(t * eig$values))
  C %*% D %*% solve(C)
}

A <- matrix(c(3, -1, -2, 2), nrow = 2, byrow = TRUE)
t <- 1
exp_matriz(A, t)

resolver.sistema.diferencial <- function(A, condiciones.iniciales, t.puntos) {
  nombres <- sapply(condiciones.iniciales, function(x) as.character(x[[1]]))
  B <- matrix(as.numeric(sapply(condiciones.iniciales, function(x) x[[2]])), ncol = 1)
  sol.mat <- sapply(t.puntos, function(t) exp_matriz(A, t) %*% B)
  if (is.vector(sol.mat)) sol.mat <- matrix(sol.mat, nrow = nrow(A))
  rownames(sol.mat) <- nombres
  soluciones <- lapply(seq_len(nrow(sol.mat)), function(i) as.numeric(sol.mat[i, ]))
  names(soluciones) <- nombres
  list(tiempos = t.puntos, soluciones = soluciones)
}

A <- matrix(c(3, -1, -2, 2), nrow = 2, byrow = TRUE)
condiciones <- list(list("x1", 90), list("x2", 150))
t.grid <- seq(0, 1, by = 0.25)
res <- resolver.sistema.diferencial(A, condiciones, t.grid)
data.frame(t = res$tiempos, do.call(cbind, res$soluciones))

A <- matrix(c(3,-2,0,-2,3,0,0,0,5), 3, 3, byrow = TRUE)
cond <- list(list("x1",2), list("x2",1), list("x3",3))
ts <- seq(0, 1, by = 0.25)
res <- resolver.sistema.puntos(A, cond, ts)
data.frame(t = res$tiempos, do.call(cbind, res$soluciones))

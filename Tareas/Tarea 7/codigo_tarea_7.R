

punto.fijo.vv <- function(p0, F, tol = 1e-8, N = 100){
  x.prev <- as.numeric(p0)
  m <- length(p0)
  x.next <- x.prev
  for (i in 1:N) {
    for (j in 1:m) {
      x.next[j] <- F[[j]](x.prev)
    }
    if(abs(max(x.next - x.prev)) < tol){
      return(x.next)
    }
    x.prev <- x.next
  }
  return(paste0("Numero maximo de iteraciones excedido: ", N))
}


f1 <- function(x) (1/3)*cos(x[2]*x[3]) + 1/6
f2 <- function(x) (1/9)*sqrt(x[1]^2 + sin(x[3]) + 1.06) - 0.1
f3 <- function(x) -(1/20)*exp(-x[1]*x[2]) - (10*pi - 3)/60

F <- list(f1, f2, f3)

p0 <- c(0, 0, 0)
sol <- punto.fijo.vv(p0, F)
sol

evaluar.componente <- function(gj, x){
  return(gj(x))
}

evaluar.funcion.vectorial <- function(G, x){
  vapply(G, evaluar.componente, numeric(1), x = x)
}

calcular.jacobiano <- function(G, x, h = sqrt(.Machine$double.eps)){
  x <- as.numeric(x)
  m <- length(x)
  J <- matrix(0, nrow = m, ncol = m)
  for (j in 1:m) {
    ej <- rep(0, m)
    ej[j] <- 1
    fp <- evaluar.funcion.vectorial(G, x + h * ej)
    fm <- evaluar.funcion.vectorial(G, x - h * ej)
    J[, j] <- (fp - fm)/(2*h)
  }
  return(J)
}

newton.raphson.vv <- function(p0, F, tol = 1e-8, N = 50){
  x.prev <- as.numeric(p0)
  for (i in 1:N) {
    J <- calcular.jacobiano(F, x.prev)
    Fx <- evaluar.funcion.vectorial(F, x.prev)
    delta <- try(solve(J, Fx), silent = TRUE)
    if(inherits(delta, "try-error") || any(!is.finite(delta))) {
      stop("Falla al resolver el sistema lineal: jacobiana singular o mal condicionada")
    }
    x.next <- x.prev - as.numeric(delta)
    
    if(abs(max(x.next - x.prev)) < tol){
      return(x.next)
    }
    x.prev <- x.next
  }
  return("Numero maximo de iteraciones excedido")
}

r1 <- function(x) x[1]^2 + x[2]^2 + 0.6*x[2] - 0.16
r2 <- function(x) x[1]^2 - x[2]^2 + x[1] -1.6*x[2] - 0.14

R <- list(r1, r2)
p0 <- c(0.6, 0.25)

sol <- newton.raphson.vv(p0, R)
sol

newton.kantorovich.vv <- function(p0, F, tol = 1e-8, N = 50){
  x.prev <- as.numeric(p0)
  J <- calcular.jacobiano(F, p0)
  for (i in 1:N) {
    Fx <- evaluar.funcion.vectorial(F, x.prev)
    delta <- try(solve(J, Fx), silent = TRUE)
    if(inherits(delta, "try-error") || any(!is.finite(delta))) {
      stop("Falla al resolver el sistema lineal: jacobiana singular o mal condicionada")
    }
    x.next <- x.prev - as.numeric(delta)
    
    if(abs(max(x.next - x.prev)) < tol){
      return(x.next)
    }
    x.prev <- x.next
  }
  return("Numero maximo de iteraciones excedido")
}

r1 <- function(x) x[1]^2 + x[2]^2 + 0.6*x[2] - 0.16
r2 <- function(x) x[1]^2 - x[2]^2 + x[1] -1.6*x[2] - 0.14

R <- list(r1, r2)
p0 <- c(0.6, 0.25)

sol <- newton.kantorovich.vv(p0, R)
sol

# Algoritmos de derivacion numerica vistos en clase

Der2Puntos <- function(f, x_0, h) {
  return((f(x_0 + h) - f(x_0)) / h)
}

Der3PuntosA <- function(f, x_0, h) {
  return(((-3 / 2) * f(x_0) + 2 * f(x_0 + h) - (1 / 2) * f(x_0 + 2 * h)) /
           h)
}

Der3PuntosB <- function(f, x_0, h){
  return(((-1/2)*f(x_0 - h) + (1/2)*f(x_0 + h))/h)
}

Der5PuntosA <- function(f, x_0, h){
  return(( (-25/12)*f(x_0) + 4*f(x_0 + h) - 3*f(x_0 + 2*h) +
             (4/3)*f(x_0 + 3*h) - (1/4)*f(x_0 + 4*h) ) / h)
}

Der5PuntosB <- function(f, x_0, h){
  return(( f(x_0 - 2*h) - 8*f(x_0 - h) + 8*f(x_0 + h) - f(x_0 + 2*h) ) / (12*h))
}

Richardson <- function(N, x0, h, n){
  Q <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in 1:n) {
    Q[i,1] <- N(x0, h/(2^(i - 1)))
  }
  for (i in 2:n) {
    for (j in 2:i) {
      Q[i,j] <- (4^(j-1)*Q[i,j-1] - Q[i-1,j-1])/(4^(j-1)-1)
    }
  }
  return(Q[n,n])
}

Trapecio <- function(f, a, b){
  return((b-a)*(f(a) + f(b))/2)
}

Simpson <- function(f, a, b){
  return(((b-a)/6)*(f(a) + 4*f((a+b)/2)+f(b)))
}

SimpsonCompuesto <- function(f, a, b, m){
  h <- (b-a)/(2 * m)
  xi0 <- f(a)+f(b)
  xi1 <- 0 #suma de nodos impares
  xi2 <- 0 #suma de nodos pares
  for(j in 1:(2 * m - 1)){
    x <- a+j*h
    if(j%%2 == 0){
      xi2 <- xi2 + f(x)
    } else{
      xi1 <- xi1 + f(x)
    }
  }
  return((h/3)*(xi0 + 2*xi2 + 4*xi1))
}

Romberg <- function(f, a, b, n){
  Q <- matrix(NA_real_, nrow = n, ncol = n)
  h <- b-a
  Q[1,1] <- (h/2)*(f(a) + f(b))
  for (i in 2:n){
    suma <- 0
    for (k in 1:(2^(i-2))){
      suma <- suma + f(a + (k - 1/2)*h)
    }
    Q[i,1] <- 0.5*(Q[i-1,1] + h*suma)
    h <- h/2
  }
  
  for (i in 2:n) {
    for (j in 2:i) {
      Q[i,j] <- (4^(j-1)*Q[i,j-1] - Q[i-1,j-1])/(4^(j-1)-1)
    }
  }
  return(Q[n,n])
}

IntegralDobleSCT <- function(f, g1, g2, a, b, n, m){
  h <- (b-a)/(2*n)
  xi0 <- 0 # inicial + final
  xi1 <- 0 # impares
  xi2 <- 0 # pares
  for (i in 0:(2*n)) {
    xw <- a + i*h
    yi <- SimpsonCompuesto(function(y) f(xw, y), g1(xw), g2(xw), m)
    if(i == 0 || i == 2*n){
      xi0 <- xi0 + yi
    } else if(i %% 2 == 0){
      xi2 <- xi2 + yi
    } else {
      xi1 <- xi1 + yi
    }
  }
  return((h/3)*(xi0+4*xi1+2*xi2))
}

TrapecioExtendido <- function(f, a, b, n){
  h <- (b-a)/n
  suma <- 0
  for(i in 1:(n-1)){
    suma <- suma + f(a + i*h)
  }
  
  return((h/2)*(f(a) + f(b) + 2*suma))
}

coef.fourier <- function(f, k, tipo = "coseno") {

  integrando <- switch(
    tipo,
    "constante" = function(x) f(x),
    "coseno"    = function(x) f(x) * cos(k * x),
    "seno"      = function(x) f(x) * sin(k * x),
    stop("El parámetro 'tipo' debe ser 'constante', 'coseno' o 'seno'.")
  )
  
  # Cálculo del coeficiente según el tipo
  if (tipo == "constante") {
    return((1 / (2 * pi)) * SimpsonCompuesto(integrando, -pi, pi, 200))
  } else {
    return((1 / pi) * SimpsonCompuesto(integrando, -pi, pi, 200))
  }
}

serie.fourier <- function(f, n, x) {
  a0 <- coef.fourier(f, 0, tipo = "constante")
  suma <- a0
  for (k in 1:n) {
    a_k <- coef.fourier(f, k, tipo = "coseno")
    b_k <- coef.fourier(f, k, tipo = "seno")
    suma <- suma + a_k * cos(k * x) + b_k * sin(k * x)
  }
  return(suma)
}


# serie.fourier <- function(f, n, x){
#   serie <- 0
#   for (k in 0:(2*n-1)) {
#     if(k == 0){
#       polinomio_k <- 1/sqrt(2*pi)
#     } else{
#       polinomio_k <- (1/sqrt(pi)) * cos(k*x)
#     }
#     
#     serie <- serie + coef.fourier(f, k) * polinomio_k
#   }
#   return(serie)
# }

funcion.prueba <- function(x) {
  if(-pi <= x & x < 0){
    return (-1)
  } else if (0 <= x & x <= pi){
    return(1)
  } else{
    return(NULL)
  }
}

serie.fourier(funcion.prueba, 3, 6)




TrapecioExtendido(f, 0, 2, 1000)



Simpson(f, 0, 2)
Trapecio(f, 0, 2)
SimpsonCompuesto(f, 0, 2, 100)
Romberg(f, 0, 2, 5)
IntegralDobleSCT(function(x, y) exp(-1*(x+y)), function(x) x^2, function(x) sqrt(x), 0, 1, 10, 10)

f <- function(x) x^4

Richardson(function(x, h) Der2Puntos(prueba, x, h), 2.19, 0.01, 10)

prueba <- function(x) (x^3)*exp(x^2) - sin(x)

Der5PuntosA(prueba, 69, 1)
Der5PuntosB(prueba, 69, 1)

# Regla del Trapecio Extendida para integrales dobles tipo I
trapecio_doble <- function(f, a, b, c, d, nx, ny) {
  hx <- (b - a) / nx
  hy <- (d - c) / ny
  suma <- 0
  
  for (i in 0:nx) {
    x <- a + i * hx
    for (j in 0:ny) {
      y <- c + j * hy
      peso <- 1
      
      if (i == 0 || i == nx) peso <- peso / 2
      if (j == 0 || j == ny) peso <- peso / 2
      
      suma <- suma + peso * f(x, y)
    }
  }
  
  I <- hx * hy * suma
  return(I)
}

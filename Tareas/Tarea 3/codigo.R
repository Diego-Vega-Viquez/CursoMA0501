#funciones

############
#Punto fijo#
############

punto.fijo <- function(p0, tol, n, g) {
  i <- 1
  p0_temp <- p0
  while (i <= n) {
    p <- g(p0_temp)
    if (abs(p - p0_temp) < tol) {
      return(p)
    }
    i <- i + 1
    p0_temp <- p
  }
  retun(NULL)
}

punto.fijo.rec <- function(p0, tol, n, g) {
  p1 <- g(p0)
  if (abs(p1 - p0) < tol || n < 1) {
    if (n >= 1) {
      return(p1)
    }
    else{
      return(Inf)
    }
  }
  else{
    return(punto.fijo.rec(p1, tol, n - 1, g))
  }
}

pba <- function(x)
  sqrt(1 + 1 / x)

###########
#Biseccion#
###########

biseccion <- function(a, b, tol, n, g) {
  i <- 1
  a1 <- a
  b1 <- b
  if (g(a) * g(b) > 0) {
    return("No se cumplen las hipotesis")
  }
  else{
    while (i <= n) {
      x <- (a1 + b1) / 2
      if (g(a1) * g(x) > 0) {
        a1 <- x
      }
      else{
        b1 <- x
      }
      if (abs(b1 - a1) < tol) {
        return(x)
      }
      i <- i + 1
    }
    return(NULL)
  }
}

biseccion.rec <- function(a, b, tol, n, g) {
  a1 <- a
  b1 <- b
  x <- (a + b) / 2
  if (abs(b1 - a1) < tol || n < 1) {
    if (n >= 1) {
      return(x)
    }
    else{
      return(x)
    }
  }
  else {
    if (g(a) * g(x) > 0) {
      return(biseccion.rec(x, b1, tol, n - 1, g))
    }
    else {
      return(biseccion.rec(a1, x, tol, n - 1, g))
    }
  }
}

f <- function(x) (cos(x))^2 - 2 * sin(x)
df <- function(x) - 2 * cos(x) * sin(x) - 2 * sin(x)

biseccion(0, 1, 1e-8, 100, f)
biseccion.rec(0, 1, 1e-8, 100, f)

################
#newton-raphson#
################

newton.raphson <- function(x0, tol, n, f, df) {
  i <- 1
  x0_temp <- x0
  while (i <= n) {
    x <- x0_temp - f(x0_temp) / df(x0_temp)
    if (abs(x - x0_temp) < tol) {
      return(x)
    }
    i <- i + 1
    x0_temp <- x
  }
  return(NULL)
}

newton.raphson.rec <- function(x0, tol, n, f, df) {
  x0_temp <- x0
  x <- x0_temp - f(x0_temp) / df(x0_temp)
  if (abs(x - x0) < tol || n < 1) {
    if (n >= 1) {
      return(x)
    }
    else{
      return(NULL)
    }
  }
  else{
    return(newton.raphson.rec(x, tol, n - 1, f, df))
  }
}

newton.raphson(0.5, 1e-8, 100, f, df)
newton.raphson.rec(0.5, 1e-8, 100, f, df)

#########
#Secante#
#########

secante <- function(x0, x1, tol, n, f) {
  i <- 2
  x0.temp <- x0
  x1.temp <- x1
  while (i <= n) {
    x <- x1.temp - ((x0.temp - x1.temp) * f(x1.temp)) / (f(x0.temp) - f(x1.temp))
    if (abs(x - x0.temp) < tol) {
      return(x)
    }
    x0.temp <- x1.temp
    x1.temp <- x
    i <- i + 1
  }
  return(NULL)
}

secante.rec <- function(x0, x1, tol, n, f) {
  x <- x1 - ((x0 - x1) * f(x1)) / (f(x0) - f(x1))
  if (abs(x - x0) < tol || n < 1) {
    if (n >= 1) {
      return(x)
    }
    else{
      return(NULL)
    }
  }
  else{
    return(secante.rec(x1, x, tol, n - 1, f))
  }
}

g <- function(x) cos(x) - x^2

secante(0, 1, 1e-8, 100, f)
secante.rec(0, 1, 1e-8, 100, f)

############
#steffensen#
############

steffensen <- function(x0, tol, n, f) {
  i <- 2
  x0.temp <- x0
  while (i <= n) {
    x1 <- f(x0.temp)
    x2 <- f(x1)
    x <- x0.temp - (x1 - x0.temp)^2 / (x2 - 2 * x1 + x0.temp)
    if (abs(x - x0.temp) < tol) {
      return(x)
    }
    i <- i + 1
    x0.temp <- x
  }
  return(NULL)
}

steffensen.rec <- function(x0, tol, n, f) {
  x1 <- f(x0)
  x2 <- f(x1)
  x <- x0 - (x1-x0)^2/(x2-2*x1+x0)
  if(abs(x - x0) < tol || n < 1){
    if(n>=1){
      return(x)
    }
    else{
      return(NULL)
    }
  }
  else{
    return(steffensen.rec(x, tol, n-1, f))
  }
}

steffensen(1, 1e-8, 100, pba)
steffensen.rec(1, 1e-8, 100, pba)

##############
#Regula Falsi#
##############

regula.falsi <- function(a,b, tol, n, f){
  i <- 1
  a.temp <- a
  b.temp <- b
  while (i <= n) {
    x <- (a.temp*f(b.temp) - b.temp*f(a.temp))/(f(b.temp) - f(a.temp))
    if(f(a.temp)*f(x) > 0){
      a.temp <- x
    }
    else{
      b.temp <- x
    }
    if(abs(b.temp-a.temp) < tol) {
      return(x)
    }
    i <- i + 1
  }
  return(NULL)
}

regula.falsi.rec(0, 1, 1e-8, 100, f)

regula.falsi.rec <- function(a, b, tol, n, f) {
  a.temp <- a
  b.temp <- b
  x <- (a.temp*f(b.temp) - b.temp*f(a.temp))/(f(b.temp) - f(a.temp))
  if(abs(b.temp-a.temp)<tol || n<1){
    if(n>=1){
      return(x)
    }
    else{
      return(NULL)
    }
  }
  else{
    if(f(a)*f(x) > 0){
      return(regula.falsi.rec(x, b.temp, tol, n-1, f))
    }
    else{
      return(regula.falsi.rec(a.temp, x, tol, n-1, f))
    }
  }
}

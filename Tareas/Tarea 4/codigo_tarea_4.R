library(tidyverse)

# Neville: interpolacion en un punto x a partir de un vector de nodos (xi) y otro de valores (f(xi))

# Retorna una lista con el valor interpolado y la matriz completa Q (la que posee todos los polinomios utilizados para la construccion del polinomio completo)

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

equis <- c(2, 2.2, 2.4, 2.6, 2.8)
ye <- c(0.5103757, 0.5207843, 0.5104147, 0.4813306, 0.4359160)

neville(equis, ye, 2.5)

graficar.polinomio <- function(nodos,
                               a,
                               b,
                               metodo,
                               f = NULL,
                               valores = NULL,
                               df = NULL,
                               derivadas.clamped = NULL) {
  stopifnot(is.numeric(nodos), length(nodos) >= 2)
  # Validación: exactamente uno de f o valores
  if (is.null(f) == is.null(valores)) {
    stop("Debe proveer exactamente uno: 'f' (función) o 'valores' (numérico).")
  }
  # Obtener valores en nodos según el caso
  if (!is.null(f)) {
    stopifnot(is.function(f))
    valores_nodos <- f(nodos)
  } else {
    stopifnot(is.numeric(valores), length(valores) == length(nodos))
    valores_nodos <- valores
  }
  
  # Derivadas (opcional). Acepta función o vector numérico.
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
  
  # Wrapper vectorizado para el método (con o sin derivadas)
  if (!is.null(derivadas_nodos)) {
    H <- function(x)
      vapply(x, function(xx)
        metodo(nodos, valores_nodos, derivadas_nodos, xx)$valor, numeric(1))
  } else {
    if(is.null(derivadas.clamped)){
      H <- function(x)
        vapply(x, function(xx)
          metodo(nodos, valores_nodos, xx)$valor, numeric(1))
    } else{
      H <- function(x)
        vapply(x, function(xx)
          spline.sujeto(nodos, valores_nodos, derivadas.clamped, xx)$valor, numeric(1))
    }
  }
  
  # Malla y data frames
  xi <- seq(a, b, length.out = 400)
  df_plot <- data.frame(x  = xi,
                        Hx = H(xi),
                        fx = if (!is.null(f))
                          f(xi)
                        else
                          NA_real_)
  df_nodos <- data.frame(x = nodos, y = valores_nodos)
  
  # Gráfico: con f (dos curvas) o solo interpolación
  p <- ggplot(df_plot, aes(x = x))
  if (!is.null(f)) {
    p <- p +
      geom_line(aes(y = fx, color = "Original"), linewidth = 1) +
      geom_line(aes(y = Hx, color = "Interpolación"),
                linewidth = 1,
                linetype = "dashed") +
      scale_color_manual(values = c(
        "Original" = "blue",
        "Interpolación" = "red"
      ))
  } else {
    p <- p +
      geom_line(aes(y = Hx, color = "Interpolación"), linewidth = 1) +
      scale_color_manual(values = c("Interpolación" = "red"))
  }
  p +
    geom_point(
      data = df_nodos,
      aes(x = x, y = y),
      shape = 21,
      size = 3,
      fill = "white"
    ) +
    labs(title = paste0("Interpolación por ", deparse(substitute(metodo))), y = "Valor", color = "Serie") +
    theme_minimal(base_size = 14)
}

graficar.polinomio(equis, 2, 3, neville, valores = ye)

nodos.x <- c(-2, -1, 0, 1, 2)
func <- function(x)
  3^x
der <- function(x)
  3^x * log(3)
nodos.y <- func(nodos.x)
deriv.y <- der(nodos.x)

neville(nodos.x, nodos.y, 0.3)

graficar.polinomio(nodos.x, -3, 3, neville, f = func)


# Lagrange.Newton: interpolacion en un punto x function# Lagrange.Newton: interpolacion en un punto x a partir de un vector de nodos y otro de valores

# Retorna una lista con el valor interpolado, la matriz para la construccion de los coeficientes, y un vector de los coeficientes del polinomio de lagrange

lagrange.newton <- function(nodos, valores, x) {
  stopifnot(is.numeric(nodos),
            is.numeric(valores),
            length(nodos) == length(valores))
  n <- length(nodos)
  Q <- matrix(NA_real_, nrow = n, ncol = n)
  
  Q[, 1] <- valores
  
  # Construccion de la tabla de diferencias de Newton
  for (i in 2:n) {
    for (j in 2:i) {
      numerador <- Q[i, j - 1] - Q[i - 1, j - 1]
      denominador <- nodos[i] - nodos[i - j + 1]
      Q[i, j] <-  numerador / denominador
    }
  }
  coeficientes <- diag(Q)
  
  valor = Q[1, 1]
  producto = 1
  
  for (i in 2:n) {
    producto <- producto * (x - nodos[i - 1])
    valor <- valor + coeficientes[i] * producto
  }
  return(list(
    valor = valor,
    tabla = Q,
    coeficientes = coeficientes
  ))
}

lagrange.newton(equis, ye, 2.5)
graficar.polinomio(nodos.x, func, -3, 3, lagrange.newton)

# Hermite.newton: interpolacion en un punto x a partir de un vector de nodos, otro de valores y otro de derivadas

# Retorna una lista con el valor interpolado, la matriz para la construccion de los  y un vector de los coeficientes

hermite.newton <- function(nodos, valores, derivadas, x) {
  stopifnot(
    is.numeric(nodos),
    is.numeric(valores),
    is.numeric(derivadas),
    length(nodos) == length(valores),
    length(nodos) == length(derivadas)
  )
  n <- length(nodos)
  
  Z <- numeric(2 * n)
  Q <- matrix(0, nrow = 2 * n, ncol = 2 * n)
  
  
  # Set-up inicial de la matriz
  for (i in 1:n) {
    z0 <- 2 * i - 1
    z1 <- 2 * i
    Z[z0] <- nodos[i]
    Z[z1] <- nodos[i]
    Q[z0, 1] <- valores[i]
    Q[z1, 1] <- valores[i]
    Q[z1, 2] <- derivadas[i]
    if (i != 1) {
      Q[z0, 2] <- (Q[z0, 1] - Q[z0 - 1, 1]) / (Z[z0] - Z[z0 - 1])
    }
  }
  # Rellenar el resto de la matriz a partir de estos valores
  for (i in 3:(2 * n)) {
    for (j in 3:i) {
      Q[i, j] <- (Q[i, j - 1] - Q[i - 1, j - 1]) / (Z[i] - Z[i - j + 1])
    }
  }
  
  coeficientes <- diag(Q)
  
  valor <- Q[1, 1]
  producto <- 1
  
  for (i in 2:(2 * n)) {
    producto <- producto * (x - Z[i - 1])
    valor <- valor + coeficientes[i] * producto
  }
  return(list(
    valor = valor,
    tabla = Q,
    coeficientes = coeficientes
  ))
}

hermite.newton(nodos.x, nodos.y, deriv.y, 0.5)

graficar.polinomio(nodos.x, func, -5, 5, hermite.newton, der)
graficar.polinomio(c(0, pi / 2, pi, 3 * pi / 2, 2 * pi), function(x)
  sin(x), -0, 8, hermite.newton, function(x)
    cos(x))

graficar.polinomio(c(1, 5, 10, 20), function(x)
  x^10 - x^5 + 1000, 0, 30, hermite.newton, function(x)
    10 * x^9 - 5 * x^4)

# spline.natural: funcion que calcula la interpolacion por splines a partir de unos nodos y sus valores

# Retorna los valores de los coeficientes a, b, c, d de cada una de las n-1 ecuaciones generadas
spline.natural <- function(nodos, valores, x) {
  # Note que los valores de a corresponden a los valores de los nodos en la funcion, por lo que cuando aparezca el vector "valores" se debe entender que equivale al vector "a"
  n <- length(nodos)
  
  h <- numeric(n - 1)
  alfa <- numeric(n - 1)
  alfa[1] <- 0
  
  # Paso 1 y 2: definir los h's y alfas
  for (i in 1:(n - 1)) {
    h[i] <- nodos[i + 1] - nodos[i]
    if (i != 1) {
      alfa[i] <- (3 / h[i]) * (valores[i + 1] - valores[i]) - (3 / h[i - 1]) *
        (valores[i] - valores[i - 1])
    }
  }
  
  #Paso 3: definir valores iniciales de l, m, y z
  l <- numeric(n) # creo que el tamano de esto puede ser n-1
  l[1] <- 1
  m <- numeric(n - 1)
  m[1] <- 0
  z <- numeric(n) # creo que el tamano de esto puede ser n-1
  z[1] <- 0
  
  # Paso 4: rellenar vectores l, m, z
  for (i in 2:(n - 1)) {
    l[i] <- 2 * (nodos[i + 1] - nodos[i - 1]) - h[i - 1] * m[i - 1]
    m[i] <- h[i] / l[i]
    z[i] <- (alfa[i] - h[i - 1] * z[i - 1]) / l[i]
  }
  
  # Paso 5: definir valores finales
  l[n] <- 1 #creo que esto no hace falta
  z[n] <- 0 #esto tampoco
  c <- numeric(n)
  c[n] <- 0
  b <- numeric(n)
  d <- numeric(n)
  
  # Paso 6: sustitucion hacia atras
  for (j in (n - 1):1) {
    c[j] <- z[j] - m[j] * c[j + 1]
    b[j] <- (valores[j + 1] - valores[j]) / h[j] - (h[j] / 3) * (c[j + 1] +
                                                                   2 * c[j])
    d[j] <- (c[j + 1] - c[j]) / (3 * h[j])
  }
  
  # Paso extra: evaluar la interpolacion en el punto x especificado
  
  ## Encontramos los dos nodos que estan prensando al intervalo
  indice <- NULL
  for (i in 1:(n - 1)) {
    if (x >= nodos[i] && x < nodos[i + 1]) {
      indice <- i
    }
  }
  if(x == nodos[n]){
    indice <- n
  }
  if (is.null(indice)) {
    return("El valor de interpolacion debe estar entre dos nodos")
  }
  ## evaluamos en la funcion asociada
  
  valor <- valores[indice] + b[indice] * (x - nodos[indice]) + c[indice] * (x - nodos[indice])^2 + d[indice] * (x - nodos[indice])^3
  
  return(list(
    a = valores,
    b = b,
    c = c,
    d = d,
    valor = valor
  ))
}
spline.natural(c(-2, -1, 0, 1, 2), c(exp(-2), exp(-1), 1, exp(1), exp(2)), 1.1)



# spline.sujeto: funcion que calcula la interpolacion por splines a partir de unos nodos y sus valores

# Retorna los valores de los coeficientes a, b, c, d de cada una de las n-1 ecuaciones generadas
spline.sujeto <- function(nodos, valores, derivadas, x) {
  # Note que los valores de a corresponden a los valores de los nodos en la funcion, por lo que cuando aparezca el vector "valores" se debe entender que equivale al vector "a"
  stopifnot(length(derivadas) == 2)
  n <- length(nodos)
  
  h <- numeric(n - 1)
  alfa <- numeric(n)
  
  # Paso 1 y 2: definir los h's y alfas
  for (i in 1:(n - 1)) {
    h[i] <- nodos[i + 1] - nodos[i]
    if (i != 1) {
      alfa[i] <- (3 / h[i]) * (valores[i + 1] - valores[i]) - (3 / h[i - 1]) *
        (valores[i] - valores[i - 1])
    }
  }
  
  alfa[1] <- 3 * ((valores[2] - valores[1]) / h[1] - derivadas[1])
  alfa[n] <- 3 * (derivadas[2] - (valores[n] - valores[n - 1]) / h[n - 1])
  
  #Paso 3: definir valores iniciales de l, m, y z
  l <- numeric(n) # creo que el tamano de esto puede ser n-1
  l[1] <- 2 * h[1]
  m <- numeric(n - 1)
  m[1] <- 1 / 2
  z <- numeric(n) # creo que el tamano de esto puede ser n-1
  z[1] <- alfa[1] / l[1]
  
  # Paso 4: rellenar vectores l, m, z
  for (i in 2:(n - 1)) {
    l[i] <- 2 * (nodos[i + 1] - nodos[i - 1]) - h[i - 1] * m[i - 1]
    m[i] <- h[i] / l[i]
    z[i] <- (alfa[i] - h[i - 1] * z[i - 1]) / l[i]
  }
  
  # Paso 5: definir valores finales
  l[n] <- h[n - 1] * (2 - m[n - 1])
  z[n] <- (alfa[n] - h[n - 1] * z[n - 1]) / l[n]
  c <- numeric(n)
  c[n] <- z[n]
  b <- numeric(n)
  d <- numeric(n)
  
  # Paso 6: sustitucion hacia atras
  for (j in (n - 1):1) {
    c[j] <- z[j] - m[j] * c[j + 1]
    b[j] <- (valores[j + 1] - valores[j]) / h[j] - (h[j] / 3) * (c[j + 1] + 2 * c[j])
    d[j] <- (c[j + 1] - c[j]) / (3 * h[j])
  }
  
  # Paso extra: evaluar la interpolacion en el punto x especificado
  
  ## Encontramos los dos nodos que estan prensando al intervalo
  indice <- NULL
  for (i in 1:(n - 1)) {
    if (x >= nodos[i] && x < nodos[i + 1]) {
      indice <- i
    }
  }
  if (x == nodos[n]) {
    indice <- n
  }
  if (is.null(indice)) {
    return("El valor de interpolacion debe estar entre dos nodos")
  }
  ## evaluamos en la funcion asociada
  
  valor <- valores[indice] + b[indice] * (x - nodos[indice]) + c[indice] * (x - nodos[indice])^2 + d[indice] * (x - nodos[indice])^3
  
  return(list(
    a = valores,
    b = b,
    c = c,
    d = d,
    valor = valor
  ))
}

spline.sujeto(c(1, 2, 3), c(2, 3, 5), c(2, 1), 2.5)
graficar.polinomio(nodos = c(1, 2, 3), 1, 3, spline.sujeto, valores = c(2, 3, 5), derivadas.clamped = c(2,1))


graficar.polinomio(c(-2, -1, 0, 1, 2), -2, 2, neville, function(x) tanh(x))
graficar.polinomio(c(-2, -1, 0, 1, 2), -2, 2, lagrange.newton, function(x) tanh(x))
graficar.polinomio(c(-2, -1, 0, 1, 2), -2, 2, hermite.newton, f = function(x) tanh(x), df = function(x) 1/cosh(x))
graficar.polinomio(c(-2, -1, 0, 1, 2), -2, 2, spline.natural, function(x) tanh(x))
graficar.polinomio(c(-2, -1, 0, 1, 2), -2, 2, spline.sujeto, f = function(x) tanh(x), derivadas.clamped = c(1/cosh(-2), 1/cosh(2)))


sistema.minimos <- function(nodos, valores, n){
  stopifnot(is.numeric(nodos), is.numeric(valores),
            length(nodos) == length(valores),
            is.numeric(n), n >= 0, n < length(nodos))
  m <- length(nodos)
  A <- matrix(NA_real_, nrow = n + 1, ncol = n + 1)
  B <- numeric(n + 1)
  for (j in 0:n){
    for (k in 0:n) {
      A[j + 1, k + 1] <- sum(nodos^(j + k)) 
    }
    B[j + 1] <- sum(nodos^(j) * valores)
  }
  return(list(A = A, B = B))
}

ajuste.minimos.cuadrados <- function(nodos, valores, n, a, b){
  stopifnot(is.numeric(nodos), is.numeric(valores),
            length(nodos) == length(valores),
            is.numeric(a), is.numeric(b), a < b,
            is.numeric(n), n >= 0)
  sistema <- sistema.minimos(nodos, valores, n)
  
  coef <- as.numeric(solve(sistema$A, sistema$B))
  
  # Polinomio ajustado
  f_hat <- function(z) vapply(z, function(zz) sum(coef * zz^(0:n)), numeric(1))
  
  # Datos para el gráfico
  xi <- seq(a, b, length.out = 400)
  df_plot <- data.frame(x = xi, Hx = f_hat(xi))
  df_pts  <- data.frame(x = nodos, y = valores)
  
  # Gráfico 
  p <- ggplot(df_plot, aes(x = x)) +
    geom_line(aes(y = Hx, color = "Ajuste (MC)"),
              linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("Ajuste (MC)" = "red")) +
    geom_point(data = df_pts, aes(x = x, y = y),
               shape = 21, size = 3, fill = "white") +
    labs(title = paste0("Interpolación por ajuste de\nmínimos cuadrados (grado ", n, ")"),
         y = "Valor", color = "Serie") +
    theme_minimal(base_size = 14)
  
  print(p)
  invisible(list(coeficientes = coef, f_ajuste = f_hat, grafico = p))
}


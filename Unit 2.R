# simple table

P <- c(1, 2, 3, 4, 5)   # basınç (Kg/cm²)
C <- c(1, 1, 2, 2, 4)   # sıkışma (mm)

simple_table <- rbind(P, C)
colnames(simple_table) <- paste0("obs", 1:5)   # isteğe bağlı
simple_table

# frequency table

t <- table(P, C) # absolute frequency

addmargins(t) # with the sums

rel_tab <- prop.table(t) # relative frequency


# marginal distribution
abs_tab <- table(P, C)

rowSums(abs_tab)    # P'nin marjinal (mutlak) dağılımı
colSums(abs_tab)    # C'nin marjinal (mutlak) dağılımı


rel_tab <- prop.table(abs_tab) # relative

rowSums(rel_tab)    # P: 0.2 0.2 0.2 0.2 0.2
colSums(rel_tab)    # C: 0.4 0.4 0.2


# conditional dimentions

tab <- matrix(c( 2,  7, 12, 10,  4,
                 5, 14, 23, 15,  7,
                 12, 31, 23,  8,  3,
                 20, 18,  8,  2,  1),
              nrow = 4, byrow = TRUE,
              dimnames = list(X = 0:3, Y = -2:2))

# Y | X = 2
table(tab["2", ])

# X | Y = -1
table(tab[, "-1"])

# moment of order (r,s) considering the point (a,b):

moment_of_order <- function(tab, r, s, a, b){
  rows <- nrow(tab)
  cols <- ncol(tab)
  
  x <- as.numeric(rownames(tab))
  y <- as.numeric(colnames(tab))
  
  total <- 0
  
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      total <- total + (x[i] - a)^r*(y[j]-b)^s*tab[i,j]
    }
  }
  
  n <- sum(tab)
  return(total/n)
  
}

moment_of_order(t, 1,2,4,5)


# ordinary moment when (a,b) = (0,0)
# central moment when (a,b) = (x_u, y_u) mean of x and mean of y


mean_x <- function(tab){
  x <- as.numeric(row.names(tab))
  x_frequencies <- rowSums(tab)
  return(sum(x*x_frequencies)/sum(x_frequencies))
  
}

mean_y <- function(tab){
  y <- as.numeric(colnames(tab))
  y_frequencies <- colSums(tab)
  return(sum(y*y_frequencies)/sum(y_frequencies))
  
}


# m00 = 1, m10 = mean of x, m01 = mean of y -> ordinary moments

# center of gravity: (m10, m01) (ordinary moments)

# covariance from a frequency table: sum of x_i*y_j*f_ij - x_u*y_u

moment_of_order(t, 1, 1, mean_x(t), mean_y(t))


# s^2_n = rs = 20 central moment

moment_of_order(t, 2, 0, mean_x(t), mean_y(t))

# M12(3,2)

moment_of_order(t, 1, 2, 3, 2)

# independency

isIndependent <- function(tab){
  n   <- sum(tab)
  f_i <- rowSums(tab) / n
  f_j <- colSums(tab) / n
  
  for (i in 1:nrow(tab)){
    for (j in 1:ncol(tab)){
      if (abs(tab[i, j] / n - f_i[i] * f_j[j]) > 1e-9){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# correlation

correlation <- function(tab){
  sxy <- moment_of_order(tab, 1, 1, mean_x(tab), mean_y(tab))
  sx  <- sqrt(moment_of_order(tab, 2, 0, mean_x(tab), mean_y(tab)))
  sy  <- sqrt(moment_of_order(tab, 0, 2, mean_x(tab), mean_y(tab)))
  sxy / (sx * sy)
}


#If X and Y are independent then s2_n (X + Y ) = s2_n (X−Y ) = s2_n (X ) + s2_n (Y ).

tab <- matrix(c(10, 6, 0, 0,
                0, 0, 0, 6,
                0, 0, 3, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(X = 1:3, Y = 1:4))

# X, Y'ye fonksiyonel bağımlı mı? Her sütunda tek sıfırdan farklı hücre
depends_X_on_Y <- function(tab){
  all(apply(tab, 2, function(col) sum(col != 0) == 1))
}

# Y, X'e bağımlı mı? Her satırda tek sıfırdan farklı hücre
depends_Y_on_X <- function(tab){
  all(apply(tab, 1, function(row) sum(row != 0) == 1))
}

depends_X_on_Y(tab)   # TRUE
depends_Y_on_X(tab)   # FALSE



#############################################
# Modelling Regression lines
#############################################

# least squares objective function
# y = f(x) for some xi
# general regression curve of Y on X is defined as the function that assigns to each value xi
# of the variable X, the mean of the variable Y /xi.











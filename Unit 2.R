options(digits = 15)  
l4 <- 10000
# simple table

P <- c(1, 2, 3, 4, 4, 5)   # basınç (Kg/cm²)
C <- c(1, 1, 2, 2, 2, 4)   # sıkışma (mm)


# frequency table

t <- table(P, C) # absolute frequency
t
addmargins(t) # with the sums

rel_tab <- prop.table(t) # relative frequency

###########################
# marginal distributions
###########################

rowSums(t)    # P'nin marjinal (mutlak) dağılımı
colSums(t)    # C'nin marjinal (mutlak) dağılımı


rel_tab <- prop.table(t) # relative

rowSums(rel_tab)    # P: 0.2 0.2 0.2 0.2 0.2
colSums(rel_tab)    # C: 0.4 0.4 0.2

##########################
# conditional dimentions
##########################

tab <- matrix(c( 2,  7, 12, 10,  4,
                 5, 14, 23, 15,  7,
                 12, 31, 23,  8,  3,
                 20, 18,  8,  2,  1),
              nrow = 4, byrow = TRUE,
              dimnames = list(X = 0:3, Y = -2:2))

addmargins(tab)

tab
# Y | X = 2
tab["2", ] # conditional dimention of Y conditioned at x = 2

# X | Y = -1
tab[, "-1"]

#######################################################
# moment of order (r,s) considering the point (a,b):
#######################################################

# tab is a matrix or table(data)
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

# ordinary moment m_rs when (a,b) = (0,0): moment_of_order(t, r, s, 0, 0)

# central moment µ_rs when (a,b) = (x_u, y_u) mean of x and mean of y: moment_of_order(t, r, s, mean_x(t), mean_y(t))


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

# center of gravity: (m10, m01) (ordinary moments) or (mean of x, mean of y)

#####################################################################
# covariance from a frequency table: sum of x_i*y_j*f_ij - x_u*y_u
#####################################################################

covariance <- function(tab) {
  moment_of_order(tab, 1, 1, mean_x(tab), mean_y(tab)) 
}


# s^2_n = rs = 20 of central moment (variance of x)

variance_of_x <- function(tab){
  moment_of_order(tab, 2, 0, mean_x(tab), mean_y(tab))
}


variance_of_y <- function(tab){
  moment_of_order(tab, 0, 2, mean_x(tab), mean_y(tab))
}

# M12(3,2)
# moment_of_order(t, 1, 2, 3, 2)

##################################
# independency
##################################

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

#####################
# correlation
#####################


correlation <- function(tab){
  sxy <- moment_of_order(tab, 1, 1, mean_x(tab), mean_y(tab)) #covariance
  sx  <- sqrt(moment_of_order(tab, 2, 0, mean_x(tab), mean_y(tab))) # sd of x
  sy  <- sqrt(moment_of_order(tab, 0, 2, mean_x(tab), mean_y(tab))) # sd of y
  sxy / (sx * sy)
}


#If X and Y are independent then s2_n (X + Y ) = s2_n (X−Y ) = s2_n (X ) + s2_n (Y ).

tab <- matrix(c(10, 6, 0, 0,
                0, 0, 0, 6,
                0, 8, 3, 0),
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


# ============================================================
# Modelling and Regression
# ============================================================


# ===============================================================
# variables and explanations
# correlation: a measure of the degree of dependence between variables
# Regression aims to find an approximate model of the dependence between variables



#===============================================================
# functions depending on x and y vectors
#===============================================================

###################################################################################################
# residual vector ei = yi - yi_est where yi_est = f(xi) given by y = f(x) is the fitting function
###################################################################################################

residual_vector<- function(y, yest){
  return(y - yest)
}

####################################################
# mean squared error mse: average of errors made
####################################################

sse <- function(z, zest, f = rep(1, length(z))){
  sum((z-zest)^2)
}

mse <- function(y, yest, f = rep(1, length(y))){
  return(sum((y-yest)^2 * f) / sum(f))
}

##########################################################
# coefficient of determination (R^2): measure of fit
# R^2 = 1 - (MSE / s^2n(Y))
# R^2 = Vex / variance of y
# R^2 = 1 means perfect fit ( variance explained by the fit is equal to the variance of y)
####################################################

r_squared <- function(y, yest, f = rep(1, length(y))){
  e <- y - yest
  mse_cal <- mse(y,yest, f)
  mean_y <- sum(f * y) / sum(f)
  var_y <- sum((y - mean_y)^2 * f) / sum(f)
  
  return(1- mse_cal / var_y)
}

##########################################
# variance explained by the fit:
# Vex = R^2 * S^2n(Y)
##########################################

variance_explained <- function( y, yest, f){
  mean_y <- sum(f * y) / sum(f)
  var_y <- sum((y - mean_y)^2 * f) / sum(f)
  r_2 <- r_squared(y, yest, f)
  return(r_2 * var_y)
}


#===========Simple Linear Regression============
# Ry/x: regression line of Y on X:
# Ry/x = mean(y) + Cov(x,y) / variance(X)* (x - mean(x))
# Rx/y: regression line of X on Y:
# Rx/y = mean(x) + Cov(x,y) / variance(Y) * (y - mean(y))
#===============================================

# Weighted mean of x with frequencies f.
mean_f <- function(x, f = rep(1, length(x))) {
  sum(x * f) / sum(f)
}

# Variance s_n^2(X) = (sum f*x^2)/n - mean^2 (divides by n, not n-1).
var_n <- function(x, f = rep(1, length(x))) {
  sum(f * (x - mean_f(x, f))^2) / sum(f)
}

# Covariance Cov(X,Y) = m11 - mean(x)*mean(y) (divides by n).
cov_n <- function(x, y, f = rep(1, length(x))) {
  sum(f * x * y) / sum(f) - mean_f(x, f) * mean_f(y, f)
}
#==========================================

# Slope b of the line Y on X: b = Cov(X,Y) / s_n^2(X).
slope_y_on_x <- function(x, y, f = rep(1, length(x))) {
  cov_n(x, y, f) / var_n(x, f)
}

# Slope b' of the line X on Y: b' = Cov(X,Y) / s_n^2(Y).
slope_x_on_y <- function(x, y, f = rep(1, length(x))) {
  cov_n(x, y, f) / var_n(y, f)
}

# Regression line Y on X: y = a + b*x with a = mean(y) - b*mean(x).
# Returns the coefficients a, b.
reg_y_on_x <- function(x, y, f = rep(1, length(x))) {
  b <- slope_y_on_x(x, y, f)
  a <- mean_f(y, f) - b * mean_f(x, f)
  c(a = a, b = b)
}

# Regression line X on Y: x = a' + b'*y with a' = mean(x) - b'*mean(y).
# Returns the coefficients a', b'.
reg_x_on_y <- function(x, y, f = rep(1, length(x))) {
  b <- slope_x_on_y(x, y, f)
  a <- mean_f(x, f) - b * mean_f(y, f)
  c(a = a, b = b)
}


#===========Simple Linear Regression============
# Linear regression functions calculated on a table
#===============================================

x <- c(1,2,3,4,5)
y <- c(3,7,2,5,7)

xytab <- table(x, y)

# slope b of the line y on x
slope_y_on_x_table <- function(tab){
  cov <- covariance(tab)
  varx <- variance_of_x(tab)
  cov/varx
}

# slope b' of the line x on y
slope_x_on_y_table <- function(tab){
  cov <- covariance(tab)
  vary <- variance_of_y(tab)
  cov/vary
}

# regression line of y on x calculated on a table
reg_y_on_x_table <- function(tab) {
  b <- slope_y_on_x_table(tab)
  a <- mean_y(tab) - b * mean_x(tab)
  c(a = a, b = b)
}


# regression line X on Y calculated on a table
# returns the coefficients a', b'.
reg_x_on_y_table <- function(tab) {
  b <- slope_x_on_y_table(tab)
  a <- mean_x(tab) - b * mean_y(tab)
  c(a = a, b = b)
}

#===============================================
# pearson's line of coefficient (r): cov(x,y) / sn(x)*sn(y)
# 0<r≤1 →direct linear correlation.
# -1≤r<0 →inverse linear correlation.
# r=0 →There is no linear relationship between the variables
# (uncorrelated variables).
# r=1 or r=-1 perfect linear correlation (direct or inverse)
#===============================================

# calculated on a table
pearson_linear_coefficient_table<- function(tab){
  c <- covariance(tab)
  sd_x <- sqrt(variance_of_x(tab))
  sd_y <- sqrt(variance_of_y(tab))
  c/(sd_x * sd_y)
}

# for x and y lines
pearson_linear_coefficient <- function(x, y, f = rep(1, length(x))){
  c <- cov_n(x, y, f)
  sd_x <- sqrt(var_n(x, f))
  sd_y <- sqrt(var_n(y, f))
  c / (sd_x * sd_y)
}

#===============================================
# MSE = variance(Y) - b^2 * variance(x)
# Since MSE = (1−R2)s2n(Y ) for the line y on x
#===============================================

mse_y_on_x_shortcut <- function(tab){
  variance_of_y(tab) * (1 - pearson_linear_coefficient_table(tab)^2)
}

#===============================================
# r² = b · b′ ≤ 1.
# Cov = b · variance(x)
# σ²ᵧ (variance(y)) = Cov / b′ 


# We say that the regression line RY /X is ”better”than the regression line RX /Y if
#          MSE_Y/X ≤ MSE_X/Y ⇐⇒s 2n (Y) ≤s2n(X )


# regression lines pass through the center of gravity (mean x, meany)


intersection_of_lines <- function(a, b, a_prime, b_prime) {
  # y = a + b*x ,  x = a' + b'*y  ->  x = a' + b'*(a + b*x)
  x <- (a_prime + b_prime * a) / (1 - b * b_prime)
  y <- a + b * x
  c(x = x, y = y)
}

# ---------- Multiple linear regression ----------

# Design matrix M with a column of ones followed by the predictors (each predictor = one column).
design_matrix <- function(...) {
  cbind(1, ...)
}

# Least squares solution of the normal equations: (M'M) A = M'Y.
# Returns the coefficients (a0, a1, ..., ak).
multi_fit <- function(M, y) {
  as.vector(solve(t(M) %*% M, t(M) %*% y))
}

# Fitted values: Yest = M * A.
multi_predict <- function(M, A) {
  as.vector(M %*% A)
}

# Residual variance, i.e. the MSE of the multiple fit: MSE = SSE / n.
multi_mse <- function(M, y, A) {
  mse(y, multi_predict(M, A))
}

# Coefficient of determination of the multiple fit: R^2 = 1 - MSE / s_n^2(Y).
multi_r2 <- function(M, y, A) {
  r_squared(y, multi_predict(M, A))
}

# Fit of a plane z = a + b*x + c*y; returns coefficients and SSE, MSE, R^2.
fit_plane <- function(x, y, z) {
  M <- design_matrix(x, y)
  A <- multi_fit(M, z)
  zest <- multi_predict(M, A)
  list(coef = setNames(A, c("a", "b", "c")),
       SSE = sse(z, zest), MSE = mse(z, zest), R2 = r_squared(z, zest))
}

# =======================ex 1 from the slides=======================

# z = a + bx + cy
x <- c(0, 1, 0, -1, 0, 1, 1, 2)
y <- c(0, 1, 1,  1, 3, 2, 2, 0)
z <- c(2, 3, 4,  5, 7, 4, 5, 0)


M <- design_matrix(x, y)
M

A <- multi_fit(M,z)
A # a, b, c

zest <- multi_predict(M,A)
trunc(zest * 10000) / 10000 # predicted values for z

sse <- function(z, zest, f = rep(1, length(z))){
  sum((z-zest)^2)
}
sse(z, zest)

resz <- z- zest

var_n(resz) # residual variance = mse 

multi_mse(M,z, A)

r_squared(z, zest)

# ======================= ex 2 =======================
x <- c(5, 5, 5, 5, 10, 10, 10, 10, 15, 15, 15, 15, 20, 20, 20, 20)
y <- c(1, 2, 3, 4,  1,  2,  3,  4,  1,  2,  3,  4,  1,  2,  3,  4)
z <- c(28, 30, 48, 74, 29, 50, 57, 42, 20, 24, 31, 47, 9, 18, 22, 31)

M <- design_matrix(x, y)
M

A <- multi_fit(M,z)
trunc(A* 10000) / 10000


zest <- multi_predict(M,A)
zest
var_n(zest) # variance of estimated z values
sse(z, zest) # sum of squared errors

r_squared(z, zest)

#======================= ex 3 =======================

x <- c(0, 1, 0, -1, 0, 1, 1, 2)
x_sin <- sin(x)
x_cos <- cos(x)

y <- c(0, 1, 1, 1, 3, 2, 2, 0)

M <- cbind(x_sin, x_cos)
M

t(M) %*% M 
t(M) %*% y

A <- multi_fit(M,y)
A

yest <- multi_predict(M,A)
yest

r_squared(y,yest) # returns how good the fit is : 0.17 very poor


# part 2 

M <- design_matrix(x, x^2)
M

A <- multi_fit(M,y)
A # c, b, a

yest <- multi_predict(M,A)
yest

r_squared(y, yest)


# ---------- Variable Transformation ----------

# Exponential Fit y = ae^(bx), ln(y) = lna + bx
# Y = ln(y)
# A = ln(a)
# Y = A + bx

x <- c(0, 1, 2, 3, 6)
y <- c(7, 5, 4, 3.5, 3)

M <- cbind(1, x ) # ln(a) b
M

A <- multi_fit(M, log(y))
A

lna <- A[1]
lna
a <- exp(lna)
a
b <- A[2]
b

yest <- a * exp(b * x)
yest
mse(y, yest)
r_squared(y, yest) # 0.8122

# hyperbolic fit: y = 1 / (a+bx)
# 1/y = a + bx
# Y = 1/y
# solve Y = a + bx

# ex with the same x and y vectors


M <-cbind(1, x)
M
y_swapped <- 1 / y

A <- multi_fit(M,y_swapped)
A

Yest <- multi_predict(M,A)
Yest

yest <- 1 / Yest
yest
mse(y, yest)
r_squared(y,yest) # 0.8572


# linear fit y = a + bx

M <- cbind(1, x)
M

A <- multi_fit(M, y)
A

yest <- multi_predict(M,A)
yest
r_squared(y, yest) # 0.7 


# ex: hypercolic fit

x <- c(1, 2, 3, 4, 5)
y <- c(1, 0.5, 0.33, 0.25, 0.2)
y_swapped <- 1 / y
M <-cbind(1, x)
M
A <- multi_fit(M, y_swapped)
A
Yest <- multi_predict(M,A)
Yest
yest <- 1 / Yest
yest
r_squared(y, yest)

# ex: the fit y = ax^b
# lny = lna + b * lnx

x <- c(1, 2, 3, 4, 5)
y <- c(0.5, 2, 4.5, 8, 12.5)

lnx <- log(x)
lny <- log(y)

M <- cbind(1, lnx) # lna , b
M
A <- multi_fit(M, lny)
A

a <- exp(A[1])
a
b <- A[2]
b

yest <- a * x^b
yest
r_squared(y , yest)


# ex slide 110

x <- c(1,2,3,4,5)
y <- c(1,1,2,4,8)

# linear fit y = a + bx

M <- cbind(1,x)

A <- multi_fit(M, y)
A

yest <- multi_predict(M,A)
r_squared(y, yest) # 0.8304: a convinient fit

# the fit y= a*b^x
# lny = lna + x * lnb

M <- cbind(1, x) 

lny <- log(y)

A <- multi_fit(M, lny) # returns the values lna and lnb

a <- exp(A[1])
b <- exp(A[2])

yest <- a * (b^x)
r_squared(y, yest) # 0.9620, a better approximation


yest_on_x_6_linear <- A[1] + A[2] * 6
yest_on_x_6_linear # 8.3

yest_on_x_6 <- a* b^6
yest_on_x_6 # 12.1257


# ex slide 111

# linear regression y = a + bx

tab111 <- matrix(c( 2, 0, 0,
                    1, 3, 2,
                    1, 3, 2,
                    2, 0, 0),
              nrow = 4, byrow = TRUE,
              dimnames = list(X = c(20, 30, 40, 50), Y = 0:2))
tab111

# returning to avectors 
xs <- as.numeric(rownames(tab111))   # 20 30 40 50
ys <- as.numeric(colnames(tab111))   # 0 1 2

g <- expand.grid(x = xs, y = ys)     # every (x, y) combination
f <- as.vector(tab111)               # its frequency, in the same order

x <- rep(g$x, f)
y <- rep(g$y, f)

# y = a + bx


isIndependent(tab111)



















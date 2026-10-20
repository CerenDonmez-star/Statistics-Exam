# ============================================================
# Unit 2 - Statistical Modelling: from "Coefficient of Determination" on
# All variances / covariances use the population formula (divide by n),
# as in the slides (s_n^2, m11, m20, m02).
# Every function accepts optional frequencies f (default: all 1).
# ============================================================


# ---------- Basic building blocks ----------

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


# ---------- MSE, R^2 and explained variance ----------

# Residual vector e_i = y_i - yest_i.
residuals_fit <- function(y, yest) {
  y - yest
}

# Sum of squared errors: SSE = sum(e_i^2).
sse <- function(y, yest, f = rep(1, length(y))) {
  sum(f * (y - yest)^2)
}

# Mean squared error: MSE = sum(e_i^2) / n.
mse <- function(y, yest, f = rep(1, length(y))) {
  sse(y, yest, f) / sum(f)
}

# Coefficient of determination: R^2 = 1 - MSE / s_n^2(Y).
# Can be negative if the model fits worse than the horizontal line y = mean(y).
r_squared <- function(y, yest, f = rep(1, length(y))) {
  1 - mse(y, yest, f) / var_n(y, f)
}

# Variance explained by the fit: Vex = R^2 * s_n^2(Y)  (so s_n^2(Y) = MSE + Vex).
variance_explained <- function(y, yest, f = rep(1, length(y))) {
  r_squared(y, yest, f) * var_n(y, f)
}

# Check of the decomposition s_n^2(Y) = MSE + Vex (should be TRUE).
check_variance_decomposition <- function(y, yest, f = rep(1, length(y))) {
  isTRUE(all.equal(var_n(y, f), mse(y, yest, f) + variance_explained(y, yest, f)))
}


# ---------- Simple linear regression ----------

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

# Normal equations of Y on X in matrix form: [n sum(x); sum(x) sum(x^2)] (a,b)' = (sum(y), sum(xy))'
# Solves them numerically; gives the same result as reg_y_on_x.
normal_eq_y_on_x <- function(x, y, f = rep(1, length(x))) {
  A <- matrix(c(sum(f), sum(f * x), sum(f * x), sum(f * x^2)), nrow = 2, byrow = TRUE)
  rhs <- c(sum(f * y), sum(f * x * y))
  sol <- solve(A, rhs)
  names(sol) <- c("a", "b")
  sol
}

# Normal equations of X on Y in matrix form (swap the roles of x and y).
normal_eq_x_on_y <- function(x, y, f = rep(1, length(x))) {
  normal_eq_y_on_x(y, x, f)
}

# Predicted values of the line Y on X.
predict_line <- function(coef, x) {
  coef["a"] + coef["b"] * x
}


# ---------- Linear correlation coefficient ----------

# Pearson's linear correlation coefficient: r = Cov(X,Y) / (s_n(X) * s_n(Y)).
pearson_r <- function(x, y, f = rep(1, length(x))) {
  cov_n(x, y, f) / (sqrt(var_n(x, f)) * sqrt(var_n(y, f)))
}

# r from the two slopes: r^2 = b*b', r has the sign of the covariance.
r_from_slopes <- function(b, b_prime) {
  sign(b) * sqrt(b * b_prime)
}

# Slope of Y on X from r: b = r * s_n(Y) / s_n(X).
slope_from_r <- function(r, sx, sy) {
  r * sy / sx
}

# MSE of the linear fit Y on X from r: MSE = s_n^2(Y) * (1 - r^2).
mse_linear_from_r <- function(y, r, f = rep(1, length(y))) {
  var_n(y, f) * (1 - r^2)
}

# MSE of the linear fit from the slope: MSE = s_n^2(Y) - b^2 * s_n^2(X).
mse_linear_from_b <- function(x, y, f = rep(1, length(x))) {
  b <- slope_y_on_x(x, y, f)
  var_n(y, f) - b^2 * var_n(x, f)
}

# In the linear case R^2 = r^2 (TRUE if it holds for the data).
check_r2_equals_r_squared <- function(x, y, f = rep(1, length(x))) {
  co <- reg_y_on_x(x, y, f)
  isTRUE(all.equal(r_squared(y, co["a"] + co["b"] * x, f), pearson_r(x, y, f)^2))
}

# Is the line Y on X better than X on Y?  (MSE_Y/X <= MSE_X/Y  <=>  s_n^2(Y) <= s_n^2(X))
y_on_x_is_better <- function(x, y, f = rep(1, length(x))) {
  var_n(y, f) <= var_n(x, f)
}

# Both regression lines meet at the centre of gravity (mean(x), mean(y)).
centre_of_gravity <- function(x, y, f = rep(1, length(x))) {
  c(x = mean_f(x, f), y = mean_f(y, f))
}


# ---------- Regression from the two lines (slide tasks) ----------

# Given the lines y = a + b*x and x = a' + b'*y, find the centre of gravity (intersection).
intersection_of_lines <- function(a, b, a_prime, b_prime) {
  # y = a + b*x ,  x = a' + b'*y  ->  x = a' + b'*(a + b*x)
  x <- (a_prime + b_prime * a) / (1 - b * b_prime)
  y <- a + b * x
  c(x = x, y = y)
}

# r, s_n(Y) and Cov from both regression lines and s_n(X):
# r = sign(b)*sqrt(b*b'), s_n(Y) = b*s_n(X)/r (from b = r*sy/sx), Cov = b*s_n^2(X).
stats_from_lines <- function(b, b_prime, sx) {
  r <- sign(b) * sqrt(b * b_prime)
  sy <- abs(b) * sx / abs(r)
  list(r = r, sy = sy, var_y = sy^2, cov = b * sx^2)
}

# Non-central moments m20, m02, m11 from the centre of gravity and s_n, Cov:
# m20 = s_n^2(X) + mean(x)^2, m02 = s_n^2(Y) + mean(y)^2, m11 = Cov + mean(x)*mean(y).
raw_moments_from_stats <- function(mx, my, sx, sy, cov) {
  c(m20 = sx^2 + mx^2, m02 = sy^2 + my^2, m11 = cov + mx * my)
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


# ---------- Fits with arbitrary basis functions ----------

# Fit y = sum_k c_k * g_k(x) with a list of basis functions g_k.
# Works for y = a*sin(x) + b*cos(x), the parabola y = a*x^2 + b*x + c, etc.
fit_basis <- function(x, y, basis) {
  M <- sapply(basis, function(g) rep_len(g(x), length(x)))   # constants are recycled to length(x)
  A <- multi_fit(M, y)
  yest <- as.vector(M %*% A)
  list(coef = A, yest = yest, SSE = sse(y, yest), MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Fit of the parabola y = a + b*x + c*x^2 (returns a, b, c, SSE, MSE, R^2).
fit_parabola <- function(x, y) {
  fit_basis(x, y, list(function(x) rep(1, length(x)), function(x) x, function(x) x^2))
}


# ---------- Variable transformations (linearisation) ----------
# In all of them the line is fitted in the transformed variables and the
# MSE and R^2 are computed with the ORIGINAL y and the fitted y (as in the slides).

# Exponential fit y = a*exp(b*x): fit ln(y) = A + b*x, a = exp(A).
fit_exponential <- function(x, y) {
  co <- reg_y_on_x(x, log(y))
  a <- exp(co["a"]); b <- co["b"]
  yest <- a * exp(b * x)
  list(a = unname(a), b = unname(b), yest = yest, SSE = sse(y, yest),
       MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Hyperbolic fit y = 1/(a + b*x): fit Y = 1/y = a + b*x.
fit_hyperbolic <- function(x, y) {
  co <- reg_y_on_x(x, 1 / y)
  yest <- 1 / (co["a"] + co["b"] * x)
  list(a = unname(co["a"]), b = unname(co["b"]), yest = yest, SSE = sse(y, yest),
       MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Power fit y = a*x^b: fit ln(y) = ln(a) + b*ln(x), a = exp(A).
fit_power <- function(x, y) {
  co <- reg_y_on_x(log(x), log(y))
  a <- exp(co["a"]); b <- co["b"]
  yest <- a * x^b
  list(a = unname(a), b = unname(b), yest = yest, SSE = sse(y, yest),
       MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Exponential fit in the form y = a*b^x: fit ln(y) = ln(a) + x*ln(b).
fit_exp_base <- function(x, y) {
  co <- reg_y_on_x(x, log(y))
  a <- exp(co["a"]); b <- exp(co["b"])
  yest <- a * b^x
  list(a = unname(a), b = unname(b), yest = yest, SSE = sse(y, yest),
       MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Linear fit packed in the same output format (to compare models easily).
fit_linear <- function(x, y) {
  co <- reg_y_on_x(x, y)
  yest <- co["a"] + co["b"] * x
  list(a = unname(co["a"]), b = unname(co["b"]), yest = yest, SSE = sse(y, yest),
       MSE = mse(y, yest), R2 = r_squared(y, yest))
}

# Compare models: table with SSE, MSE and R^2 of each fit; the best has the largest R^2.
compare_fits <- function(fits) {
  data.frame(model = names(fits),
             SSE = sapply(fits, function(m) m$SSE),
             MSE = sapply(fits, function(m) m$MSE),
             R2  = sapply(fits, function(m) m$R2),
             row.names = NULL)
}


# ============================================================
# Examples with the data of the slides
# ============================================================

# Exponential / hyperbolic / linear fit (pages 103-107).
x <- c(0, 1, 2, 3, 6)
y <- c(7, 5, 4, 3.5, 3)

fits <- list(exponential = fit_exponential(x, y),   # a = 5.8846, b = -0.1307, MSE = 0.3755, R2 = 0.8123
             hyperbolic  = fit_hyperbolic(x, y),    # a = 0.1690, b = +0.0306 (slide typo: sign), MSE = 0.2855, R2 = 0.8572
             linear      = fit_linear(x, y))        # SSE = 2.63
compare_fits(fits)                                  # hyperbolic is the best fit

pearson_r(x, y)                                     # r = -0.8585 (r^2 = 0.7370)

# Plane fit (page 99): z = 2.2045 - 1.0568x + 1.6591y, SSE = 0.6477, R2 = 0.9794
xp <- c(0, 1, 0, -1, 0, 1, 1, 2)
yp <- c(0, 1, 1, 1, 3, 2, 2, 0)
zp <- c(2, 3, 4, 5, 7, 4, 5, 0)
fit_plane(xp, yp, zp)

# Steel bars task (page 100): z = a + b*x + c*y
xs <- rep(c(5, 10, 15, 20), each = 4)
ys <- rep(1:4, times = 4)
zs <- c(28, 30, 48, 74, 29, 50, 57, 42, 20, 24, 31, 47, 9, 18, 22, 31)
fit_plane(xs, ys, zs)

# Sin/cos and parabola fits (page 101).
xq <- c(0, 1, 0, -1, 0, 1, 1, 2)
yq <- c(0, 1, 1, 1, 3, 2, 2, 0)
fit_basis(xq, yq, list(sin, cos))
fit_parabola(xq, yq)

# Traffic tunnel task (page 87): scatter plot and r.
D <- c(43, 55, 40, 52, 39, 33, 50, 33, 44, 21)
V <- c(27.0, 23.8, 30.7, 24.0, 34.8, 41.4, 27.0, 40.4, 31.7, 51.2)
plot(D, V, main = "Traffic flow", xlab = "D (veh/km)", ylab = "V (km/h)")
pearson_r(D, V)                                     # close to -0.968

# Lines x - 2y = 4 (X on Y: x = 4 + 2y) and 2x - 9y = 8 (Y on X: y = (2x - 8)/9).
stats_from_lines(b = 2 / 9, b_prime = 2, sx = 3)

# Lines y = x/2 + 2 and x = 2y - 4 (page 90): centre of gravity and r.
intersection_of_lines(a = 2, b = 1 / 2, a_prime = -4, b_prime = 2)


#ex

D <- c(43, 55, 40, 52, 39, 33, 50, 33, 44, 21)
V <- c(27.0, 23.8, 30.7, 24.0, 34.8, 41.4, 27.0, 40.4, 31.7, 51.2)

# 1) Scatter diagram (with the regression line of V on D)
plot(D, V, xlab = "D (veh/km)", ylab = "V (km/h)")
co <- reg_y_on_x(D, V)
abline(a = co["a"], b = co["b"], col = "red")

# 2) Choose r from {0.968, -0.968, -0.198}:
#    the points go down as D increases (inverse) and lie close to a line,
#    so the correct value is r = -0.968.

# 3) Verify by calculating r
r <- pearson_r(D, V)
r                                   # -0.9683


# Extra information from the other functions
co                                  # a = 67.485, b = -0.8362
cov_n(D, V)                         # covariance (negative -> inverse relation)
r^2                                 # 0.9376
check_r2_equals_r_squared(D, V)     # TRUE: R^2 = r^2 in the linear case
mse(V, co["a"] + co["b"] * D)       # 4.390



# ex 2 


N   <- 10
vx  <- 9          # sigma_x^2 = 9

# --- Which line is which? ---
# x - 2y = 4   ->  x = 4 + 2y        : X on Y, slope b' = 2
# 2x - 9y = 8  ->  y = (2x - 8)/9    : Y on X, slope b  = 2/9
# (the other assignment gives r^2 = 9/4 > 1, impossible)
b       <- 2 / 9
b_prime <- 2

# --- 1) r, variance of Y and covariance ---
res <- stats_from_lines(b, b_prime, sx = sqrt(vx))
res$r        # 0.6667
res$var_y    # 1
res$cov      # 2

# Centre of gravity (intersection of both lines)
G <- intersection_of_lines(a = -8 / 9, b = 2 / 9, a_prime = 4, b_prime = 2)
G            # x = 4, y = 0

# --- 2) Remove the point (2, -1) and recompute the lines ---
# Removes one point from a data set described only by its n, means, variances and covariance.
remove_point <- function(n, mx, my, vx, vy, cov, px, py) {
  m <- raw_moments_from_stats(mx, my, sqrt(vx), sqrt(vy), cov)   # m20, m02, m11
  k <- n - 1
  mx2 <- (n * mx - px) / k
  my2 <- (n * my - py) / k
  list(n = k, mx = mx2, my = my2,
       vx  = (n * m["m20"] - px^2) / k - mx2^2,
       vy  = (n * m["m02"] - py^2) / k - my2^2,
       cov = (n * m["m11"] - px * py) / k - mx2 * my2)
}

new <- remove_point(N, G["x"], G["y"], vx, res$var_y, res$cov, px = 2, py = -1)

b_new       <- new$cov / new$vx            # slope of Y on X
a_new       <- new$my - b_new * new$mx
b_prime_new <- new$cov / new$vy            # slope of X on Y
a_prime_new <- new$mx - b_prime_new * new$my
r_new       <- new$cov / sqrt(new$vx * new$vy)

c(a = a_new, b = b_new)                    # y = -0.7662 + 0.2078 x
c(a = a_prime_new, b = b_prime_new)        # x = 4 + 2 y
r_new                                      # 0.6447

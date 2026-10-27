# =====================================================================
#  Topic 3: TIME SERIES  -  R functions for every formula in the slides
#  (Univ. de Malaga, Dept. of Applied Mathematics)
#
#  Layout of the file
#    1. Moving averages (order k, odd / even)
#    2. Decomposition: trend-cycle, E*A, seasonal indices, random part
#    3. Deseasonalisation
#    4. Trend by least squares (regression line on X/E)
#    5. Prediction (trend * seasonal index)
#    6. Autocorrelation (cycle detection)
#    7. Index numbers + deflation
#    8. Whole pipeline in one call + car-sales exercise
#
#  Every place marked   >>> SLIDE CODING TASK <<<   is code that appears
#  literally in the slides ("Code in R" slides) - run those first.
#  Expected values from the slides are written in comments (# ==> ...).
#  Only base R / stats is needed (no tidyverse): stats::filter is used
#  explicitly so it does not clash with dplyr::filter.
# =====================================================================


# =====================================================================
# 1. MOVING AVERAGES
# =====================================================================

# Moving average of order k:  Y_t = (X_t + X_{t+1} + ... + X_{t+k-1}) / k
#  * k odd : plain filter  rep(1,k)/k, result sits on the middle point
#            (order 3 -> first value is assigned to t = 2).
#  * k even: not centred, so correct it by weighting the two extremes
#            with 0.5:  c(.5, 1, ..., 1, .5)/k  (same as a 2nd MA of order 2).
#  The first/last floor(k/2) values are lost -> NA.

# p lot ( uk )
# lines ( desestM , col =4)

# best way is to use stats::filter(x, c(...))
xx <- c (2.6767 , 3.2050 , 4.8714 , 4.6275 , 5.9008 ,
         5.5200 , 5.1186 , 4.9507 , 4.6612 , 6.2572)

stats::filter(xx, c(1,1,1) / 3) # moving average of order 3
stats:: filter(xx, c(0.5, 1, 1, 1, 0.5)/4) # moving average of order 4 centered

# slide 18
x59 <- c(6, 6.2, 5.1, 4.9, 5.2, 5.8, 7)
names(x59) <- 1959:1965

# If we know the period of seasonal factors, it is advisable to use
# that value as the order for the moving averages, e.g. for monthly
# data use k = 12, for quarterly data k = 4 (four measurements in a year), etc.

# take moving average of k where k is the period for the data (mothly weekly etc.)

xx <- c (178.2 , 153.2 , 185.9 , 163.6 ,
         196.3 , 156.9 , 197.9 , 166.4 ,
         197.3 , 159.7 , 202.6 , 175.6 ,
         209.5 , 169.5 , 202.4 , 179.8 ,
         200.0 , 168.6 , 216.1 , 178.5
)
# indicating the frequency
uk <- ts ( xx , start = c (1980 , 1) , frequency = 4)
uk # X

# trend component, using moving average of order 4
tend <- stats::filter(uk, c(0.5, 1, 1, 1, 0.5)/ 4)
tend # estimates trend-cycle T*C

# eliminating Trend-cycle divide x by Tren-component
# only remains E * A
# MULTİPLİCATİVE CASE
# if additive subtract it x = T + E + C + A

est_ea <- uk / tend # E*A
trunc(est_ea* 10000) / 10000 # gives E * A (the other factors are eliminated)


# =====================================================================
# SEASONAL VARIATION INDEX (E)
# =====================================================================
# 1. unnormalized seasonal component
seasonal_component_unnormalized <- colMeans(matrix(est_ea, ncol = 4, byrow = T), na.rm = T)
seasonal_component_unnormalized

# 2. normalized seasonal component (corrected or adjusted)
# the average of seasonal averages should be equal to 1 or 100 in percentage
# in order to correct it: multiply 1 / current average or their sum / current sum

normalized_seasonal_component <- seasonal_component_unnormalized / mean(seasonal_component_unnormalized)
normalized_seasonal_component

# =====================================================================
# Random Component (A)
# =====================================================================
# E*A / E (normalized seasonal component)

random_component <- est_ea / normalized_seasonal_component
random_component


# =====================================================================
# Deseasonalisation
# =====================================================================
# X / E (after calculating seasonal component (E))
# removing seasonal component
# remain is called true or corrected trend
# try to fit a regression line on X/E

deseas <- uk / normalized_seasonal_component
trunc(deseas * 10000) / 10000


# =====================================================================
# Regression Line (linear trend of deseasonalized series)
# =====================================================================
# Solution
# To find the equation of the regression line in the form:
# Y= a + b·t, where:
# Y is the deseasonalized value of the series in a specific period.
# t is the assigned time number (1, 2, . . . , 20).
# a is the intercept of the trend.
# b is the slope that shows the change per time period.
xx <- seq_along(uk) # puts a vector 1: n where n is the length of time series vector
model <- lm(deseas ~ xx )
model
a <- coef(model)[["(Intercept)"]]   # the constant
b <- coef(model)[["xx"]]  

regL <- ts( predict.lm ( model ) , start = c (1969 , 1) ,
             frequency =4)
regL # holds the predicted values of the regression line

plot ( deseas )
lines ( regL , col =2)



# ==================== ex slide 56 ====================

ex56 <- c(33742, 35401, 36147, 38067, 34149, 36176, 37485,
       39047, 34783, 37042, 38008, 40132, 36466, 38680, 39976, 42273,
       39131, 40780, 41852, 43684, 38729, 40427, 41576, 43886, 39131,
       40394, 40956, 42959, 38714, 40062, 41152, 43460, 38695, 39780,
       40923, 44093, 40777, 41778, 43160, 45897, 41947, 44061, 44378,
       47237, 43315, 43396, 44843, 46835, 42833, 43548, 44637, 47107,
       42552, 43526, 45039, 47940, 43740, 45007, 46667, 49325, 44878,
       46234, 47055, 50318, 46354, 47260, 48883, 52605, 48527, 50237,
       51592, 55152, 50451, 52294, 54633, 58802, 53990, 55477, 57850,
       61978)

uk56 <- ts (ex56 , start = c(1969 , 1) , frequency = 4)
uk56

tend56 <- stats::filter(uk56, c(0.5, 1,1,1,0.5) / 4)
tend56 # trend cycle

ea56 <- uk56 / tend56
ea56 # E*A 

seasonal_component_unnormalized56 <- colMeans(matrix(ea56, ncol = 4, byrow = T), na.rm = T)
seasonal_component_unnormalized56

normalized_seasonal_component56 <- seasonal_component_unnormalized56 / mean(seasonal_component_unnormalized56)
normalized_seasonal_component56 # E (corrected)

deseasonalized_trend <- ex56 / normalized_seasonal_component56
deseasonalized_trend

random_component56 <- ea56 / normalized_seasonal_component56
random_component56

xx56 <- seq_along(uk56) # puts a vector 1: n where n is the length of time series vector
model <- lm(deseasonalized_trend ~ xx56 )
model

# the trend line: a + bt
a <- coef(model)[["(Intercept)"]]   # the constant
b <- coef(model)[["xx56"]]  

regL <- ts( predict.lm( model ) , start = c (1969 ,1),
            frequency =4)
regL # holds the predicted values of the regression line


# ======================= Prediction =======================

# Autocorrelation

data <- c (1.2 , 1.4 , 2.0 , 1.23 , 1.37 , 2.07 ,
           1.25 , 1.47 , 2.09)
# Function that calculates autocorrelation for a frequency n

autoCorrel <- function ( data , n ) {
  # Remove the last n data points
  xx1 <- head ( data , -n )
  # Remove the first n data points
  xx2 <- tail ( data , -n )

  sqrt( summary(lm( xx1 ~ xx2 )) $ r.squared )
}

autoCorrel(data, 3) # autocorrelation of data with freq 3
# sliding the data by 3 cycles


# ============== Index Numbers ================
# elementary or simple index number the percentage
# of variation (as a percentage or per unit) of the studied variable
# with respect to the value obtained for the base reference period
# X_a/b


# ex: If in 1970 a statistical variable X is 123 and in 1980 it is 148,
# calculate the percentage of variation with 1970 as the base period

# X_1980/1970: 148 / 123 (value at time a / value at time b)

# =====================================================================
# 5. PREDICTION  (trend extrapolated x seasonal index)
# =====================================================================

# Season (1..p) of period number t:  ((t-1) mod p) + 1
season_of <- function(t, period) ((t - 1) %% period) + 1

# Forecast for period t:
#   multiplicative:  X^_t = (a + b*t) * E_season(t)
#   additive:        X^_t = (a + b*t) + E_season(t)
# coef: a and b values of regressin line, i_idx: corrected seasonal indices, t: period of forecast
forecast_ts <- function(coef, e_idx, t, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  tr <- unname(coef["a"] + coef["b"] * t)
  e  <- e_idx[season_of(t, length(e_idx))]
  if (model == "multiplicative") tr * e else tr + e
}

forecast_ts(c(a = a,b  = b), normalized_seasonal_component56, 81:85)




# 2. DECOMPOSITION  X = T*E*C*A   (multiplicative)  /  X = T+E+C+A (additive)
# =====================================================================

# Trend-cycle T*C: moving average whose order = seasonal period
# (quarterly p = 4, monthly p = 12, ...).
trend_cycle <- function(x, period) moving_average(x, period)

# Remove trend-cycle -> seasonal * random part.
#   multiplicative:  E*A = X / (T*C)        additive:  E+A = X - (T*C)
seasonal_random <- function(x, tc, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  if (model == "multiplicative") x / tc else x - tc
}

# Arrange a series into a table with one column per season (rows = years).
# Pads with NA at the end if the last year is incomplete.
season_matrix <- function(x, period) {
  x <- as.numeric(x)
  pad <- (-length(x)) %% period
  matrix(c(x, rep(NA, pad)), ncol = period, byrow = TRUE)
}

# UNCORRECTED seasonal variation indices: average of every column of the
# E*A table (removes the random part A).
#   I_season = mean over years of (X / T*C),   NA values ignored.
#   They are usually written as percentages (x100).
seasonal_index_uncorrected <- function(ea, period) {
  colMeans(season_matrix(ea, period), na.rm = TRUE)
}

# CORRECTED (normalised) indices: multiplicative -> mean must be 1 (sum = p,
# i.e. 400 % for quarters):   I' = I * p / sum(I) = I / mean(I)
#   e.g. 400/399.53 * 109.03 % = 109.1606 %
# additive -> they must sum to 0:  I' = I - mean(I)
seasonal_index_corrected <- function(i_unc, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  if (model == "multiplicative") i_unc / mean(i_unc) else i_unc - mean(i_unc)
}

# Same correction written with the slides' formula  400/399.53 * I  (generic p).
seasonal_index_corrected_slide <- function(i_unc) {
  length(i_unc) / sum(i_unc) * i_unc
}

# Random component:  A = (E*A) / E   (additive: (E+A) - E).
# The seasonal indices (length p) are recycled along the series.
random_component <- function(ea, e_idx, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  if (model == "multiplicative") ea / e_idx else ea - e_idx
}

# Read an index as text:  1.0903 -> "+9.03 % w.r.t. the annual average"
interpret_index <- function(i) {
  d <- (i - 1) * 100
  sprintf("%.2f %% %s the annual average", abs(d), ifelse(d >= 0, "above", "below"))
}

# ---- >>> SLIDE CODING TASK (p.~23-42): the company "Muja" <<<
xx <- c(178.2, 153.2, 185.9, 163.6,     # 1980 Spring..Winter
        196.3, 156.9, 197.9, 166.4,     # 1981
        197.3, 159.7, 202.6, 175.6,     # 1982
        209.5, 169.5, 202.4, 179.8,     # 1983
        200.0, 168.6, 216.1, 178.5)     # 1984
# Reading the data, frequency = 4 (quarters)
uk <- ts(xx, start = c(1980, 1), frequency = 4)
uk
plot(uk)

# Trend component, using a moving average of order 4  (slide code)
tend <- stats::filter(uk, c(.5, 1, 1, 1, .5) / 4)
tend
# ==> Autumn 1980 = 172.4875, Winter 1980 = 175.2125, ... Summer 1984 = 190.9625
plot(uk)
lines(tend, col = 2)                 # trend in red

# Seasonal (E) * Random (A) = X / TC  (slide code)
est_aleM <- uk / tend
est_aleM
# ==> 1981: 1.1079 0.8764 1.1025 0.9246 ...

# Unnormalised seasonal component (slide code)
est1M <- colMeans(matrix(est_aleM, ncol = 4, byrow = T), na.rm = T)
est1M
# ==> 1.0903 0.8817 1.0859 0.9373   (sum 3.9953 = 399.53 %, should be 400 %)

# Normalise the seasonal component (slide code)
estM <- est1M / mean(est1M)
estM
# ==> 1.091606 0.882746 1.087217 0.938432  (= 109.16 %, 88.27 %, 108.72 %, 93.84 %)

# Random component (slide code)
aleM <- est_aleM / estM
aleM
# ==> e.g. Autumn 1983 = 0.9844  (accidental factors lowered X by 1.56 %)

# --- the same with my functions ---
tc    <- trend_cycle(uk, 4)
ea    <- seasonal_random(uk, tc)                    # E*A
i_unc <- seasonal_index_uncorrected(ea, 4)          # 1.0903 0.8817 1.0859 0.9373
sum(i_unc) * 100                                    # 399.53 -> must be 400
i_cor <- seasonal_index_corrected(i_unc)            # 1.0916 0.8827 1.0872 0.9384
seasonal_index_corrected_slide(i_unc)               # identical
interpret_index(i_cor[1])                           # "9.16 % above the annual average"
interpret_index(i_cor[2])                           # "11.73 % below the annual average"
a_rand <- random_component(ea, i_cor)               # random part A (mean ~ 1)

# Additive version of the same data (X = T+E+C+A): subtract, centre at 0
ea_add  <- seasonal_random(uk, tc, "additive")
e_add   <- seasonal_index_corrected(seasonal_index_uncorrected(ea_add, 4), "additive")
sum(e_add)                                          # 0


# =====================================================================
# 3. DESEASONALISATION    X / E = T*C*A   (additive: X - E)
# =====================================================================

# Divide each value by its corrected seasonal index (indices recycled).
deseasonalise <- function(x, e_idx, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  if (model == "multiplicative") x / e_idx else x - e_idx
}

# ---- >>> SLIDE CODING TASK (p.~46): deseasonalised series <<<
desestM <- uk / estM
desestM
# ==> Spring 1980: 178.2 / 1.091606 = 163.2458 ; Summer 1980: 173.5494 ...
plot(uk)
lines(desestM, col = 4)              # deseasonalised series in blue
deseasonalise(uk, i_cor)             # same, with my function


# =====================================================================
# 4. TREND BY LEAST SQUARES (regression line on the deseasonalised series)
# =====================================================================
#   Y = a + b*t ,  t = 1,2,...,n (1 = first period)
#   b = cov(t,Y)/var(t)  ,  a = mean(Y) - b*mean(t)     (cov/var with the same n)
trend_line <- function(y) {
  y <- as.numeric(y); t <- seq_along(y)
  b <- (mean(t * y) - mean(t) * mean(y)) / (mean(t^2) - mean(t)^2)
  a <- mean(y) - b * mean(t)
  c(a = a, b = b)
}

# Value of the trend line at times t:  T_t = a + b*t
trend_values <- function(coef, t) coef["a"] + coef["b"] * t

# ---- >>> SLIDE CODING TASK (p.~53): regression line of X/E <<<
# This puts a vector 1:n in xx, where n is the length of uk
xx <- seq_along(uk)
# Linear regression
model <- lm(desestM ~ xx)
model                                                  # ==> 169.2985 + 1.30027 * t
# Predicted values of the line.  (The slide writes start = c(1969,1); for the
# Muja data the series starts in 1980, so use c(1980,1).)
regL <- ts(predict.lm(model), start = c(1980, 1), frequency = 4)
regL
# ==> 170.599 171.899 173.199 174.5 | 175.8 177.1 178.4 179.701 | ... | 195.304
# Plot the deseasonalised series with the line
plot(desestM)
lines(regL, col = 2)

# --- same with my functions ---
co <- trend_line(desestM)                              # a = 169.2985, b = 1.30027
trend_values(co, 1)                                    # 170.599  (t = 1: 1.30027*1 + 169.2985)


# =====================================================================
# 5. PREDICTION  (trend extrapolated x seasonal index)
# =====================================================================

# Season (1..p) of period number t:  ((t-1) mod p) + 1
season_of <- function(t, period) ((t - 1) %% period) + 1

# Forecast for period t:
#   multiplicative:  X^_t = (a + b*t) * E_season(t)
#   additive:        X^_t = (a + b*t) + E_season(t)
forecast_ts <- function(coef, e_idx, t, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  tr <- unname(coef["a"] + coef["b"] * t)
  e  <- e_idx[season_of(t, length(e_idx))]
  if (model == "multiplicative") tr * e else tr + e
}

# Next 1980-84 data end at t = 20 -> forecast the four quarters of 1985:
forecast_ts(co, i_cor, 21:24)
# ==> 214.614 (Spring)  174.699 (Summer)  216.579 (Autumn)  188.160 (Winter)


# =====================================================================
# 6. AUTOCORRELATION  (detecting a cycle of period T)
# =====================================================================
# Linear correlation between X_t and X_{t+T}; r ~ 1  =>  cycle of period T.
# Slide version: sqrt(R^2) of lm(xx1 ~ xx2)  (always >= 0).
auto_correl <- function(data, n) {
  xx1 <- head(data, -n)            # remove the last n data points
  xx2 <- tail(data, -n)            # remove the first n data points
  sqrt(summary(lm(xx1 ~ xx2))$r.squared)
}

# Signed Pearson version:  r = cov(X_t, X_{t+T}) / (sd * sd)
auto_correl_signed <- function(data, n) cor(head(data, -n), tail(data, -n))

# Correlation for lags 1..max_lag (to find the best period quickly)
auto_correl_table <- function(data, max_lag) {
  sapply(1:max_lag, function(l) auto_correl_signed(data, l))
}

# ---- >>> SLIDE CODING TASK (p.~61-62): cycle of period 3 <<<
data <- c(1.2, 1.4, 2.0, 1.23, 1.37, 2.07, 1.25, 1.47, 2.09)
auto_correl(data, 3)               # ==> 0.9936708   (~1  => cycle of period 3)
auto_correl_signed(data, 3)        # 0.9936708
auto_correl_table(data, 5)         # lag 3 is the largest


# =====================================================================
# 7. INDEX NUMBERS  and  DEFLATION
# =====================================================================

# Simple (elementary) index of period a with base b:   X_{a/b} = X_a / X_b
# (1.2033 = +20.33 %). Vectorised: pass a whole series, base = position/name.
simple_index <- function(x_a, x_b) x_a / x_b

# All indices of a series with respect to the base period (position `base`).
index_series <- function(x, base = 1) x / x[base]

# Properties of simple indices
index_inverse   <- function(i_ab) 1 / i_ab              # X_{a/b} = 1 / X_{b/a}
index_chain     <- function(...) prod(...)              # X_{a/b}*X_{b/c}*...*X_{y/z} = X_{a/z}
index_cyclic_ok <- function(i_ab, i_bc, i_ca) isTRUE(all.equal(i_ab * i_bc * i_ca, 1))

# Change the base of an index series: I_{t/t1} = I_{t/t0} / I_{t1/t0}
rebase_index <- function(i, new_base = 1) i / i[new_base]

# Percentage of variation:  (X_{a/b} - 1) * 100
pct_change <- function(i) (i - 1) * 100

# Deflation (nominal -> real):  real value_{t0} = nominal value_t / CPI_{t/t0}
# (CPI as "so many per one": 135 % -> 1.35)
deflate <- function(nominal, cpi) nominal / cpi

# Inverse operation (real -> nominal, e.g. 1990 pesetas in 1998 pesetas)
inflate <- function(real, cpi) real * cpi

# ---- Example: X1970 = 123, X1980 = 148
simple_index(148, 123)                     # ==> 1.2033  (+20.33 %)
pct_change(simple_index(148, 123))         # 20.33

# ---- Example: eyeglass frames, base 2000
frames <- c(`2000` = 1254, `2001` = 1345, `2002` = 1408, `2003` = 1451)
index_series(frames, 1)                    # ==> 1, 1.0726, 1.1228, 1.1571
index_inverse(1.0726)                      # 0.9323 = index of 2000 with base 2001
index_chain(1.0726, 1.0468)                # example of the "modified cyclic" property

# ---- Example: salary +50 %, CPI doubled (1970-1980)
real_1970 <- deflate(150, 2)               # ==> 75  -> real salary fell by 25 %
pct_change(real_1970 / 100)                # -25

# ---- Exercise: CPI_{1998/1990} = 135 %
deflate(1000, 1.35)                        # 1) real value in 1990 of 1000 pts of 1998 = 740.74
inflate(1000, 1.35)                        # 2) real value in 1998 of 1000 pts of 1990 = 1350

# ---- Exercise: car bought in 2006 for 24,000 EUR; I90/69 = 1000 %, I90/06 = 75 %
#   I_{06/69} = I_{06/90} * I_{90/69} = (1/0.75) * 10 = 13.33
i_06_69 <- index_chain(index_inverse(0.75), 10)
deflate(24000, i_06_69)                    # ==> 1800 euros (in 1969 euros)


# =====================================================================
# 8. WHOLE PIPELINE IN ONE CALL
# =====================================================================
# Runs: trend-cycle (MA of order p) -> E*A -> uncorrected & corrected seasonal
# indices -> random part -> deseasonalised series -> trend line -> fitted trend.
decompose_ts <- function(x, period, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  x   <- as.numeric(x)
  tc  <- as.numeric(trend_cycle(x, period))
  ea  <- seasonal_random(x, tc, model)
  iu  <- seasonal_index_uncorrected(ea, period)
  ic  <- seasonal_index_corrected(iu, model)
  a   <- random_component(ea, ic, model)
  des <- deseasonalise(x, ic, model)
  co  <- trend_line(des)
  list(trend_cycle = tc, seasonal_random = ea,
       index_uncorrected = iu, index_corrected = ic,
       random = a, deseasonalised = des,
       trend_coef = co, trend = unname(co["a"] + co["b"] * seq_along(x)))
}

res <- decompose_ts(xx <- c(178.2,153.2,185.9,163.6,196.3,156.9,197.9,166.4,197.3,159.7,
                            202.6,175.6,209.5,169.5,202.4,179.8,200.0,168.6,216.1,178.5), 4)
res$index_corrected      # 1.091606 0.882746 1.087217 0.938432
res$trend_coef           # a = 169.2985, b = 1.30027


# ---- >>> SLIDE EXERCISE (p.~56): quarterly car sales from Q1 1969 <<<
# Do the whole decomposition (trend-cycle, seasonal indices, deseasonalised
# series, regression line, forecast) for these 80 quarters.
xx <- c(33742, 35401, 36147, 38067, 34149, 36176, 37485,
        39047, 34783, 37042, 38008, 40132, 36466, 38680, 39976, 42273,
        39131, 40780, 41852, 43684, 38729, 40427, 41576, 43886, 39131,
        40394, 40956, 42959, 38714, 40062, 41152, 43460, 38695, 39780,
        40923, 44093, 40777, 41778, 43160, 45897, 41947, 44061, 44378,
        47237, 43315, 43396, 44843, 46835, 42833, 43548, 44637, 47107,
        42552, 43526, 45039, 47940, 43740, 45007, 46667, 49325, 44878,
        46234, 47055, 50318, 46354, 47260, 48883, 52605, 48527, 50237,
        51592, 55152, 50451, 52294, 54633, 58802, 53990, 55477, 57850,
        61978)
cars <- ts(xx, start = c(1969, 1), frequency = 4)
plot(cars)

tend_c <- stats::filter(cars, c(.5, 1, 1, 1, .5) / 4)       # trend-cycle
lines(tend_c, col = 2)
est_ale_c <- cars / tend_c                                  # E*A
est1_c <- colMeans(matrix(est_ale_c, ncol = 4, byrow = T), na.rm = T)
est_c  <- est1_c / mean(est1_c)                             # corrected indices
est_c                       # ==> 0.9555 0.9821 1.0045 1.0579  (Q1 is the weakest quarter)
ale_c   <- est_ale_c / est_c                                # random part
desest_c <- cars / est_c                                    # deseasonalised
xt <- seq_along(cars)
model_c <- lm(desest_c ~ xt)
reg_c <- ts(predict.lm(model_c), start = c(1969, 1), frequency = 4)
plot(desest_c); lines(reg_c, col = 2)
# Forecast of the next 4 quarters (t = 81..84)
forecast_ts(trend_line(desest_c), est_c, 81:84)
# Is there a yearly cycle in the raw data?  (lag 4)
auto_correl_signed(as.numeric(cars), 4)

# Same exercise through the pipeline function:
res_c <- decompose_ts(xx, 4)
res_c$index_corrected
# Unit 1 functions

array <- c(1,2,3,3)

# creating a frequency table:

table <- table(array)

names <- names(table) # returns the elements in the table

freq <- as.numeric(table) # returns the frequencies in the same order with the elements

# creating a relative frequency table

f_table <- prop.table(table)
print(f_table)

# bar chart and pie chart

d <- c( a = 2, b=3)
print(d)
barplot(d)

pie(d)

# creating a histogram:
# hist(x,
#    main = "Distribution of Values",
#    xlab = "Values",
#    ylab = "Frequency",
#    col = "lightblue",
#    border = "black")

hist(array)



# simple arithmetic mean

my_mean <- function(v){
  sum <- sum(v)
  N <- length(v)
  m <- sum/N
  return(trunc(m * 10000) / 10000) # for 4 decimals 
}

my_mean_f <- function(x, f){
  return(sum(x*f)/sum(f))
}


# quadratic mean or root mean:

quadratic_mean <- function(v){
  s <- 0
  for (i in v) {
    s <- s + i**2
  }
  return(trunc(sqrt(s/length(v)) * 100) / 100) # two digits
}

quadratic_mean_f <- function(x, f) {
  return(sqrt(sum((x**2) * f) / sum(f)))
}


# mode:
a <- c(7, 11, 11, 8, 12, 9, 6, 20)

get_mode <- function(x) {
  frequencies <- table(x)
  return(as.numeric(names(frequencies)[frequencies == max(frequencies)]))
}

# median:

m <- median(a)

# quartiles: 0.25, 0.5, 0.75

q1 <- quantile(a, 0.25)
q2 <- quantile(a, 0.5)

# deciles: 0.1, 0.2, ... , 0.9

# percentiles: 0.01, 0.02, .. , 0.99

# Interquartile range: Q3 - Q1

R_Q <- function(x){
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  return(q3-q1)
}

# Intercentile range: P99 - P1

R_P <- function(x){
  p1 <- quantile(x, 0.01)
  q99 <- quantile(x, 0.99)
  return(p99-p1)
}

################################
# deviatian from p:
################################

deviation <- function(x, p){
  s <- 0
  for (i in x){
    s <- s + abs(i - p)
  }
  return(s/length(x))
}

deviation_f <- function(x, f, p){
  d <- sum(abs(x-p)*f)/sum(f)
  return(d)
}

################################
# mean deviation:
################################

deviation(a, mean(a))

###################################
# mean squared error or deviation:
###################################

mean_squared_dev <- function(x, p){
  s <- 0
  for (i in x){
    s <- s + (i-p)**2
  }
  
  return(s/length(x))
}

mean_squared_dev_f <- function(x, f, p){
  s <- sum(((x-p)**2) * f)
  return(s/sum(f))
}

#################################################################
# VARIANCE: 
# mean of squares of deviation concerning the mean
# population variance: sum((x - m)^2 * f) / sum(f)

# sample variance (sn**2): sum((x - m)^2 * f) / sum(f) 

# sample quasi-variance (s2): sum((x - m)^2 * f) / (sum(f)-1)
#################################################################

variance_f <- function(x, f){
   m <- my_mean_f(x, f)
   return(sum(((x-m)**2)*f)/sum(f))
}



variance <- function(x){
  m <- mean(x)
  return(sum((x-m)**2)/length(x))
}

###############################################
# standard deviation: squareroot of variance
# tells how spread data is
###############################################

standard_dev_f <- function(x, f){
  m <- my_mean_f(x, f)
  s <- sum((x-m)^2 * f) / sum(f)
  return(sqrt(s))
}

standard_dev <- function(x){
  return(sqrt(variance(x)))
}


#############################################################
# STANDARDISED VARIABLE
# compare data with different means and standard deviations
#  zᵢ = (xᵢ − x̄) / sₙ
# if |z| < 1: close to the mean
# if |z| near 2: not unusually far
# if |z| > 3: possibly and outlier
#############################################################

z_var <- function(x){
  m <- mean(x)
  s <- standard_dev(x)
  return((x-m)/s)
}

###########################################################
# pearson's coefficient of variace: CV 
# CV = sd/abs(mean)
# measures the spreadness of data based on the mean
# low CV: data is homogenous and mean is dependable
# high CV: data is heterogenous and mean is a weak measure
###########################################################

CV <- function(x){
  s <- standard_dev(x)
  return(s/abs(mean(x)))
}

#################################################################
# moment of order r with respect to the point c:
# measures the shape they measure deviation from the center
# high r value means measurement is very sensitive to far points
#################################################################

moment_of_order <- function(x, r, c){
  return(sum((x-c)**r)/length(x))
}

moment_of_order_freq <- function(x, r, c, f){
  return(sum((x-c)^r*f) / sum(f))
}

####################################
# ordinary moment of order
# c = 0
####################################


ordinary_moment <- function(x, r){
  return(sum(x**r)/length(x))
}

ordinary_moment_freq <- function(x, r, f){
  return(sum(x^r * f) / sum(f))
}

#########################################################
# central moment of order r : 
# arithmetic mean of rth powers of deviatian from mean
#########################################################

central_moment <- function(x, r){
  return(sum((x-mean(x))**r)/length(x))
}

central_moment_freq <- function(x, r, f){
  m <- my_mean_f(x, f)
  return(sum((x - m)^r * f) / sum(f))
}

##############################################
# pearson's skewness coefficient
##############################################

pearson_skewness <- function(x){
  mo <- get_mode(x)
  s <- standard_dev(x)
  me <- mean(x)
  return((me - mo) / s)
  
}

###############################################
# fisher skewness coefficient
# g₁ = μ₃ / sₙ³ (central moment of order 3 / cube of sd)
# g₁ = 0: symmetric. The two sides balance.
# g₁ > 0: positive (right) skew. A long tail to the right, i.e. a few unusually large values. Usually mean > median.
# g₁ < 0: negative (left) skew. A long tail to the left, i.e. a few unusually small values. Usually mean < median.
###############################################

fisher_g1 <- function(x){
  return(central_moment(x,3)/(standard_dev(x)**3))
}

fisher_g1_freq <- function(x, f){
  return(central_moment_freq(x, 3, f) / (standard_dev_f(x, f)^3) )
}

#########################################################################
# fisher's curtoises coefficient
# g₂ = μ₄ / sₙ⁴ − 3( (central moment of order 4 / fourth power of sd) - 3
# g2 = 0	Mesokurtic: normal distribution. extreme values about as often as a normal distribution would give
# g2 > 0 Leptokurtic: sharp peak, heavy tails.	most values packed near the mean, but more outliers than normal
# g2 < 0 Platykurtic: Flat, light tails. values spread evenly across the range, with few or no extremes
#########################################################################

fisher_g2 <- function(x){
  return(central_moment(x,4)/(standard_dev(x)**4) -3)
}

fisher_g2_freq <- function(x, f){
  return(central_moment_freq(x,4,f)/(standard_dev_f(x, f)**4) -3)
}


########################################
# outliers
# to detect a specific outlier calculate standardized variable z for whole vector
# if z>3 definetely an outlier
# interquartile range rule
########################################

iqr_outliers <- function(x) {
  q1  <- quantile(x, 0.25)
  q3  <- quantile(x, 0.75)
  iqr <- q3 - q1
  list(
    outliers = x[x < q1 - 1.5*iqr | x > q3 + 1.5*iqr],
    extreme  = x[x < q1 - 3*iqr   | x > q3 + 3*iqr]
  )
}




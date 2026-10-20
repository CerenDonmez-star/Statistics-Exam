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
a <- c(7, 11, 11, 8, 12, 9, 6)

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
  q1 <- quantile(x, 0.01)
  q3 <- quantile(x, 0.99)
  return(q3-q1)
}


# deviatian from p:

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

# mean deviation:

deviation(a, mean(a))


# mean squared error or deviation:

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



# variance: mean of squares of deviation concerning the mean

variance <- function(x, f){
   m <- my_mean_f(x, f)
   return(sum(((x-m)**2)*f)/sum(f))
}

variance_vector <- function(x){
  m <- mean(x)
  return(sum((x-m)**2)/length(x))
}

# standard deviation: squareroot of variance

standard_dev <- function(x, f){
  return(sqrt(variance(x,f)))
}

standard_dev_f <- function(x){
  return(sqrt(variance_vector(x)))
}


# population variance: sum((x - m)^2 * f) / sum(f)

# sample variance (sn2): sum((x - m)^2 * f) / (sum(f) - 1)

# sample quasi-variance (s2): sum((x - m)^2 * f) / (sum(f) - 1)


# standardised variable: compare datas with different means and standard deviations

z <- function(x){
  m <- mean(x)
  t <- table(x)
  s <- standard_dev(as.numeric(names(t)), as.numeric(t))
  return((x-m)/s)
}

# pearson's coefficient of variace: CV as percentage sd/abs(mean)

CV <- function(x){
  t <- table(x)
  s <- standard_dev(as.numeric(names(t)), as.numeric(t))
  return(s/abs(mean(x))*100)
}


# moment of order r with respect to the point c:

moment_of_order <- function(x, r, c){
  return(sum((x-c)**r)/length(x))
}



# ordinary moment of order

ordinary_moment <- function(x, r){
  return(sum(x**r)/length(x))
}

# central moment of order r : arithmetic mean of rth powers of deviatian from mean

central_moment <- function(x, r){
  return(sum((x-mean(x))**r)/length(x))
}

# pearson's skewness coefficient

A_p <- function(x){
  mo <- get_mode(x)
  s <- standard_dev_f(x)
  me <- mean(x)
  return((me - mo) / s)
  
}


# fisher skewness coefficient

fisher_g1 <- function(x){
  return(central_moment(x,3)/(standard_dev_f(x)**3))
}

# fisher's curtoises coefficient

fisher_g2 <- function(x){
  return(central_moment(x,4)/(standard_dev_f(x)**4) -3)
}
















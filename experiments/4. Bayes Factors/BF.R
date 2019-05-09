f1 <- "b0 > b1 > b2"
f2 <- "b0 < b1"
f3 <- "b0 >= b1"
f4 <- "b0 <= b1"
f5 <- "b0 == b1"
f6 <- "b0 != b1"
f6 <- "b0 - b1 < 0"

# Keep it simple; only one operator allowed!

# Split at operator
operators <- c(">", "<", ">=", "<=", "!=", "==")

# Check if illegal operators found
# (basically if any of the allowed operators occur)
library(stringr)
detect_operator <- str_detect(f1, operators)

# Check
if(all(detect_operator == FALSE)) {
  stop("Illegal operator passed")
}

# If multiple found
if(length(which(detect_operator)) > 1) {
  stop("Can only pass one operator in each hypothesis")
}

# Split at operator
operator <- operators[which(detect_operator)]
op_split <- str_split(f1, operator)[[1]] %>%
  # Trim whitespace
  trimws()

# Look at parameters in left side

# Look at parameters in right side

# Check if all variables in data


## EXAMPLE

mscale <- function(b, y) {
  (x - mu)/sd
}

# Simulate data
pr1 <- rnorm(1000, 10, 3)
pr2 <- rnorm(1000, 9, 3)
pr3 <- rnorm(1000, 5, 15)

r <- (pr2 > pr3)
l <- (pr1 > pr2)

t <- (l + r) == 2

c <- mean(t)

for(i in 1:5) print(paste0(pr1[i], ", ", pr2[i], ", ", pr3[i], ", ", t[i]))
# Lesson: this does not stack! > > >

## Scale
mp1 <- mean(pr1)
sd1 <- sd(pr1)

pr1 <- mscale(pr1, mp1, sd1)
pr2 <- mscale(pr2, mp1, sd1)

# Hypothesis:
#  u1 > u2
#  (u1 - u2) < 1
c1 <- mean(pr1 > pr2)
c2 <- mean(abs(pr1 - pr2) < .1)

# Posteriors
pos1 <- rnorm(1000, mean=10, sd=1)
pos2 <- rnorm(1000, mean=9.9, sd=1)

# Scale
mpos1 <- mean(pos1)
sdpos1 <- sd(pos1)

pos1 <- mscale(pos1, mpos1, sdpos1)
pos2 <- mscale(pos2, mpos1, sdpos1)

# Fit?
f1 <- mean(pos1 > pos2)
f2 <- mean((abs(pos1 - pos2) < 1.5))

# BF = fit / complexity
f1 / c1
f2 / c2

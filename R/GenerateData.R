
# Influential Outliers

set.seed(123)

### -----------------------------
### Main cluster
### -----------------------------
n1 <- 200
x1_main <- rnorm(n1, mean = 0, sd = 1)
x2_main <- rnorm(n1, mean = 0, sd = 1)

# True logistic model
eta1 <- -1 + 2*x1_main + 1.5*x2_main
p1   <- plogis(eta1)
y1   <- rbinom(n1, size = 1, prob = p1)

main <- data.frame(
  x1 = x1_main,
  x2 = x2_main,
  y  = y1,
  cluster = "main"
)

### -----------------------------
### Influential cluster
### -----------------------------
n2 <- 40
x1_inf <- rnorm(n2, mean = 5, sd = 0.6)   # far away in predictor space
x2_inf <- rnorm(n2, mean = 5, sd = 0.6)

# Contradict the pattern: make all y = 0
y2 <- rep(0, n2)

inf <- data.frame(
  x1 = x1_inf,
  x2 = x2_inf,
  y  = y2,
  cluster = "influential"
)

### Combine
df <- rbind(main, inf)


cols <- c("0" = "orange", "1" = "blue")

plot(df$x1, df$x2,
     col = cols[as.character(df$y)],
     pch = ifelse(df$cluster == "main", 19, 17),
     xlab = "x1",
     ylab = "x2",
     main = "Two-Predictor Binomial Data (Colored by y)")

legend("topleft",
       legend = c("y = 0", "y = 1", "main cluster", "influential cluster"),
       col    = c("orange", "blue", "black", "black"),
       pch    = c(19, 19, 19, 17),
       bty = "n")


fit_with  <- glm(y ~ x1 + x2, data = df, family = binomial)
fit_main  <- glm(y ~ x1 + x2, data = subset(df, cluster=="main"), family = binomial)

summary(fit_with)
summary(fit_main)

library(car)
influencePlot(fit_with)

## Fit logistic model
fit <- glm(y ~ x1 + x2, data = df, family = binomial)

## Create grid to evaluate predicted probabilities
x1_seq <- seq(min(df$x1), max(df$x1), length = 100)
x2_seq <- seq(min(df$x2), max(df$x2), length = 100)
grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)

## Predict probabilities on grid
grid$pred <- predict(fit, newdata = grid, type = "response")

## ---- PLOT: Data + fitted probability contours ----

cols <- c("0" = "orange", "1" = "blue")

plot(df$x1, df$x2,
     col = cols[as.character(df$y)],
     pch = ifelse(df$cluster == "main", 19, 17),
     xlab = "x1", ylab = "x2",
     main = "Data with Fitted Logistic Model Contours")

legend("topleft",
       legend = c("y = 0", "y = 1", "main cluster", "influential cluster"),
       col    = c("orange", "blue", "black", "black"),
       pch    = c(19, 19, 19, 17),
       bty = "n")

## Add probability contour lines
contour(x1_seq, x2_seq, matrix(grid$pred, 100, 100),
        levels = c(0.1, 0.25, 0.5, 0.75, 0.9),
        add = TRUE, drawlabels = TRUE, col = "darkgreen", lwd = 2)

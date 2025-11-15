### ---------------------------
### 1. IRLS logistic regression
### ---------------------------

irls_logistic <- function(X, y, max_iter = 25, tol = 1e-6) {
  p <- ncol(X)
  beta <- rep(0, p)

  paras <- matrix(NA_real_, nrow = max_iter + 1, ncol = p)
  colnames(paras) <- colnames(X)
  paras[1, ] <- beta

  converged <- FALSE

  for (k in 1:max_iter) {
    eta <- as.vector(X %*% beta)
    mu  <- 1 / (1 + exp(-eta))
    W   <- mu * (1 - mu)
    W[W == 0] <- 1e-8
    z <- eta + (y - mu) / W

    WX   <- X * W
    XtWX <- t(X) %*% WX
    XtWz <- t(X) %*% (W * z)

    beta_new <- solve(XtWX, XtWz)
    paras[k + 1, ] <- beta_new

    if (max(abs(beta_new - beta)) < tol) {
      paras <- paras[1:(k + 1), , drop = FALSE]
      converged <- TRUE
      break
    }
    beta <- beta_new
  }

  list(beta = beta_new, paras = paras, iter = k, converged = converged)
}


### ---------------------------
### 2. Run IRLS
### ---------------------------

# (Assume X and y already defined earlier)
result <- irls_logistic(X, y)
paras  <- result$paras

iters  <- 0:(nrow(paras) - 1)
xrange <- range(iters)
yrange <- range(paras)


### ---------------------------
### 3. Generate PNG frames
### ---------------------------

for (k in 1:nrow(paras)) {
  png(sprintf("irls_%02d.png", k), width = 700, height = 450)

  matplot(
    x = iters[1:k],
    y = paras[1:k, , drop = FALSE],
    type = "l",
    lty  = 1,
    lwd  = 2,
    xlim = xrange,
    ylim = yrange,
    xlab = "Iteration",
    ylab = "Parameter value",
    main = paste("IRLS Parameter Updates: iteration", iters[k])
  )

  legend("topright",
         legend = colnames(paras),
         col    = 1:ncol(paras),
         lty    = 1,
         lwd    = 2)

  dev.off()
}


### ---------------------------
### 4. Generate GIF
### ---------------------------

# detect magick command
if (system("magick -version", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0) {
  # new Windows ImageMagick syntax
  system("magick convert -delay 20 -loop 0 irls_*.png irls.gif")
} else {
  # Linux/macOS / old ImageMagick
  system("convert -delay 20 -loop 0 irls_*.png irls.gif")
}


### ---------------------------
### 5. (Optional) remove PNG frames
### ---------------------------

# unlink("irls_*.png")

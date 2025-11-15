#' IRLS for logistic regression with GIF output
#'
#' @param X Design matrix (including intercept column if desired).
#' @param y Response vector (0/1).
#' @param gif_file Output GIF file name.
#' @param max_iter Maximum number of IRLS iterations.
#' @param tol Convergence tolerance on max |beta_new - beta|.
#' @param interval Time interval between frames in seconds.
#' @param width GIF width in pixels.
#' @param height GIF height in pixels.
#'
#' @return A list with components:
#'   \item{beta}{Final coefficient estimates.}
#'   \item{paras}{Matrix of coefficients per iteration.}
#'   \item{iter}{Number of iterations run.}
#'   \item{converged}{Logical, whether convergence criterion met.}
#' @export
irls_logistic_gif <- function(X, y,
                              gif_file  = "irls.gif",
                              max_iter  = 25,
                              tol       = 1e-6,
                              interval  = 0.2,
                              width     = 700,
                              height    = 450) {
  # ensure matrix
  X <- as.matrix(X)
  y <- as.numeric(y)

  # ---- IRLS core ----
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
      converged <- TRUE
      paras <- paras[1:(k + 1), , drop = FALSE]
      break
    }

    beta <- beta_new
  }

  if (!converged) {
    # if never broke, truncate to full matrix
    paras <- paras[1:(max_iter + 1), , drop = FALSE]
    k <- max_iter
  }

  iters  <- 0:(nrow(paras) - 1)
  xrange <- range(iters)
  yrange <- range(paras)

  # ---- GIF generation using animation::saveGIF ----
  animation::saveGIF(
    expr = {
      for (i in 1:nrow(paras)) {
        graphics::matplot(
          x = iters[1:i],
          y = paras[1:i, , drop = FALSE],
          type = "l",
          lty  = 1,
          lwd  = 2,
          xlim = xrange,
          ylim = yrange,
          xlab = "Iteration",
          ylab = "Parameter value",
          main = paste("IRLS Parameter Updates – up to iteration", iters[i])
        )
        graphics::legend(
          "topright",
          legend = colnames(paras),
          col    = seq_len(ncol(paras)),
          lty    = 1,
          lwd    = 2
        )
      }
    },
    movie.name = gif_file,
    interval   = interval,
    ani.width  = width,
    ani.height = height
  )

  # return results
  list(
    beta      = beta,
    paras     = paras,
    iter      = k,
    converged = converged,
    gif_file  = gif_file
  )
}

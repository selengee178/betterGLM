#' IRLS for logistic regression and GIF generation (base R)
#'
#' Runs IRLS for logistic regression, saves a sequence of PNG frames
#' showing parameter updates, and then uses an external ImageMagick
#' command (`magick convert` or `convert`) to create a GIF.
#'
#' @param X Design matrix (include intercept column if desired).
#' @param y Binary response (0/1).
#' @param gif_file Name of the output GIF file (relative to working dir).
#' @param max_iter Maximum number of IRLS iterations.
#' @param tol Convergence tolerance on max |beta_new - beta|.
#' @param delay Delay between frames in the GIF (ImageMagick `-delay` units).
#' @param width,height Width/height of PNG frames in pixels.
#' @param cleanup Logical, whether to delete the intermediate PNG frames.
#'
#' @return A list with elements:
#'   \item{beta}{Final coefficient estimates.}
#'   \item{paras}{Matrix of coefficients at each iteration (rows = iterations).}
#'   \item{iter}{Number of iterations run.}
#'   \item{converged}{Logical, whether it converged.}
#'   \item{gif_file}{Path to the generated GIF file.}
#' @export
irls_logistic_gif <- function(X, y,
                              gif_file = "irls.gif",
                              max_iter = 25,
                              tol      = 1e-6,
                              delay    = 20,
                              width    = 700,
                              height   = 450,
                              cleanup  = TRUE) {
  # ensure basic types
  X <- as.matrix(X)
  y <- as.numeric(y)

  n <- nrow(X)
  p <- ncol(X)

  # ---- IRLS core ----
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
    # if we used all iterations without breaking, keep full matrix
    paras <- paras[1:(max_iter + 1), , drop = FALSE]
    k <- max_iter
  }

  iters  <- 0:(nrow(paras) - 1)
  xrange <- range(iters)
  yrange <- range(paras)

  # ---- 1) Create PNG frames with base R (grDevices + graphics) ----
  frame_files <- character(nrow(paras))
  for (i in seq_len(nrow(paras))) {
    frame_files[i] <- sprintf("irls_%02d.png", i)

    grDevices::png(frame_files[i], width = width, height = height)

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

    grDevices::dev.off()
  }

  # ---- 2) Call ImageMagick via system() (base R) ----
  # Try "magick convert" first (common on Windows),
  # then plain "convert" (common on Linux/macOS).
  # Use shQuote for safety.
  png_args <- paste(shQuote(frame_files), collapse = " ")
  gif_arg  <- shQuote(gif_file)

  # check which command is available
  has_magick <- (system("magick -version",
                        ignore.stdout = TRUE,
                        ignore.stderr = TRUE) == 0)
  cmd <- NULL
  if (has_magick) {
    cmd <- paste("magick convert -delay", delay, "-loop 0",
                 png_args, gif_arg)
  } else {
    has_convert <- (system("convert -version",
                           ignore.stdout = TRUE,
                           ignore.stderr = TRUE) == 0)
  }

  if (has_magick) {
    cmd <- paste("magick convert -delay", delay, "-loop 0",
                 png_args, gif_arg)
  } else if (has_convert) {
    cmd <- paste("convert -delay", delay, "-loop 0",
                 png_args, gif_arg)
  } else {
    cmd <- NULL
  }

  if (is.null(cmd)) {
    warning("Neither 'magick' nor 'convert' found. GIF not created.")
  } else {
    system(cmd)
  }

  if (cleanup) {
    unlink(frame_files)
  }

  list(
    beta      = beta,
    paras     = paras,
    iter      = k,
    converged = converged,
    gif_file  = gif_file
  )
}






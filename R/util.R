
spectral.colors <- function (n) 
{
  col <- c("#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C")
  # colors obtained as rev(brewer.pal(5, "Spectral"))
  palette <- grDevices::colorRampPalette(col)
  palette(n)
}

bl2gr.colors <- function (n) 
{
  palette <- grDevices::colorRampPalette(c("#084081", "#0868AC", "#2B8CBE", 
                                           "#4EB3D3", "#7BCCC4", "#A8DDB5", 
                                           "#CCEBC5", "#E0F3DB"), 
                                         space = "Lab")
  palette(n)
}

blue2grey.colors <- function(n) 
{
  basecol <- c("#E6E6E6", "#bcc9d1", "#6c7f97", "#3e5264")
  palette <- grDevices::colorRampPalette(basecol, space = "Lab")
  palette(n)
}

persp3D <- function(x, y, z, theta = 30, phi = 20, d = 5, expand = 2/3, 
                    xlim = range(x, finite = TRUE), 
                    ylim = range(y, finite = TRUE), 
                    zlim = range(z, finite = TRUE), 
                    levels = pretty(zlim, nlevels), nlevels = 20, 
                    color.palette = spectral.colors, border = NA, 
                    ticktype = "detailed", 
                    xlab = NULL, ylab = NULL, zlab = NULL, ...)
{
#----------------------------------------------------------------------------#  
# 3D plot, i.e. perspective plot, with different levels in different colors
#
# Example
# y <- x <- seq(-10, 10, length=60)
# f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
# z <- outer(x, y, f)
# persp3D(x, y, z, theta = 30, phi = 30, expand = 0.5)
# persp3D(x, y, z, color.palette = heat.colors, phi = 30, theta = 225, box = TRUE, border = NA, shade = .4)
# persp3D(x, y, z, color.palette = terrain.colors, phi = 30, theta = 225, box = FALSE, border = NA, shade = .4)
#
# x1 = seq(-3,3,length=50)
# x2 = seq(-3,3,length=50)
# y = function(x1, x2) sin(x1)+cos(x2)
# persp3D(x1, x2, outer(x1,x2,y), zlab="y", theta = 150, phi = 20, expand = 0.6)
#
#----------------------------------------------------------------------------#
  
  if(is.null(xlab)) 
    xlab <- if(!missing(x)) 
      deparse(substitute(x))
  else "X"
  if(is.null(ylab)) 
    ylab <- if(!missing(y)) 
      deparse(substitute(y))
  else "Y"
  if(is.null(zlab)) 
    zlab <- if(!missing(z)) 
      deparse(substitute(z))
  else "Z"
  if(missing(z))
  { if(!missing(x)) 
  { if(is.list(x)) 
  { z <- x$z
  y <- x$y
  x <- x$x }
    else 
    { z <- x
    x <- seq.int(0, 1, length.out = nrow(z)) }
  }
    else stop("no 'z' matrix specified")
  }
  else if(is.list(x))
  { y <- x$y
  x <- x$x }
  if(any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")

  # browser()

  # getting the value of the midpoint
  zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  # set colors for levels
  cols <- color.palette(length(levels)-1)
  zzz <- cut(zz, breaks = levels, include.lowest = TRUE, labels = cols)
  # plot
  out <- persp(x, y, z, theta = theta, phi = phi, d = d, expand = expand,
               col = as.character(zzz),
               xlim = xlim, ylim = ylim, zlim = zlim,
               border = border, ticktype = ticktype, 
               xlab = xlab, ylab = ylab, zlab = zlab, ...)
  # add breaks and colors for a legend
  out <- list(persp = out, levels = levels, colors = cols)
  invisible(out)
}

# This function is superseded by logsumexp() in mclust >= 6.1
# logsumexp <- function(x, v = NULL)
# { 
#   x <- as.matrix(x)
#   v <- if(is.null(v)) double(ncol(x)) else as.vector(v)
#   if(length(v) != ncol(x))
#     stop("Non-conforming arguments in logsumexp")
#   # as.vector(logsumexp_Rcpp(x,v))
#   logsumexp_Rcpp(x,v)
# }

# This function is superseded by logsumexp() in mclust >= 6.1
# softmax <- function(x, v = NULL)
# { 
#   x <- as.matrix(x)
#   v <- if(is.null(v)) double(ncol(x)) else as.vector(v)
#   if(length(v) != ncol(x))
#     stop("Non-conforming arguments in logsumexp")
#   softmax_Rcpp(x,v)
# }

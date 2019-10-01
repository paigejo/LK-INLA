overMap = function (points, database = "world", regions = ".", exact = FALSE, boundary = TRUE, 
          interior = TRUE, projection = "", parameters = NULL, orientation = NULL, 
          fill = FALSE, col = 1, plot = TRUE, add = FALSE, namesonly = FALSE, 
          xlim = NULL, ylim = NULL, wrap = FALSE, resolution = if (plot) 1 else 0, 
          type = "l", bg = par("bg"), mar = c(4.1, 4.1, par("mar")[3], 
                                              0.1), myborder = 0.01, namefield = "name", lforce = "n", 
          ...) 
{
  if (resolution > 0 && !plot) 
    stop("must have plot=TRUE if resolution is given")
  if (!fill && !boundary && !interior) 
    stop("one of boundary and interior must be TRUE")
  doproj <- !missing(projection) || !missing(parameters) || 
    !missing(orientation)
  if (doproj && !requireNamespace("mapproj", quietly = TRUE)) {
    stop("Please install the package 'mapproj' for projections.")
  }
  coordtype <- maps:::maptype(database)
  if (coordtype == "unknown") 
    stop("missing database or unknown coordinate type")
  if (doproj && coordtype != "spherical") 
    stop(paste(database, "database is not spherical; projections not allowed"))
  if (is.character(database)) 
    as.polygon = fill
  else as.polygon = TRUE
  coord <- maps:::map.poly(database, regions, exact, xlim, ylim, boundary, 
                    interior, fill, as.polygon, namefield = namefield)
  if (is.na(coord$x[1])) 
    stop("first coordinate is NA.  bad map data?")
  if (length(wrap) >= 2) {
    antarctica <- if (length(wrap) == 2) 
      -89.9
    else wrap[3]
    coord <- map.wrap.poly(coord, xlim = wrap[1:2], poly = fill, 
                           antarctica = antarctica)
  }
  if (lforce == "e") {
    coord <- map.clip.poly(coord, xlim, ylim, poly = fill)
  }
  if (plot) {
    .map.range(coord$range)
  }
  if (doproj) {
    nam <- coord$names
    coord <- mapproj::mapproject(coord, projection = projection, 
                                 parameters = parameters, orientation = orientation)
    coord$projection = projection
    coord$parameters = parameters
    coord$orientation = orientation
    coord$names <- nam
    if (!is.null(xlim) && !is.null(ylim) && lforce %in% c("s", 
                                                          "l")) {
      prange <- mapproj::mapproject(x = rep(xlim, 2), y = rep(ylim, 
                                                              each = 2))
      if (lforce == "s") {
        xlim <- c(max(prange$x[c(1, 3)]), min(prange$x[c(2, 
                                                         4)]))
        ylim <- c(max(prange$y[c(1, 2)]), min(prange$y[c(3, 
                                                         4)]))
      }
      else {
        xlim <- c(min(prange$x[c(1, 3)]), max(prange$x[c(2, 
                                                         4)]))
        ylim <- c(min(prange$y[c(1, 2)]), max(prange$y[c(3, 
                                                         4)]))
      }
    }
    if (plot && coord$error) 
      if (all(is.na(coord$x))) 
        stop("projection failed for all data")
    else warning("projection failed for some data")
  }
  if (length(wrap) == 1 && wrap) 
    coord <- map.wrap(coord)
  if (plot) {
    if (!add) {
      opar = par(bg = bg)
      if (!par("new")) 
        plot.new()
      if (is.null(xlim) || (doproj && !(lforce %in% c("s", 
                                                      "l")))) 
        xrange <- range(coord$x, na.rm = TRUE)
      else xrange <- xlim
      if (is.null(ylim) || (doproj && !(lforce %in% c("s", 
                                                      "l")))) 
        yrange <- range(coord$y, na.rm = TRUE)
      else yrange <- ylim
      if (coordtype != "spherical" || doproj) {
        aspect <- c(1, 1)
      }
      else aspect <- c(cos((mean(yrange) * pi)/180), 1)
      d <- c(diff(xrange), diff(yrange)) * (1 + 2 * myborder) * 
        aspect
      if (coordtype != "spherical" || doproj) {
        plot.window(xrange, yrange, asp = 1/aspect[1])
      }
      else {
        par(mar = mar)
        p <- par("fin") - as.vector(matrix(c(0, 1, 1, 
                                             0, 0, 1, 1, 0), nrow = 2) %*% par("mai"))
        par(pin = p)
        p <- par("pin")
        p <- d * min(p/d)
        par(pin = p)
        d <- d * myborder + ((p/min(p/d) - d)/2)/aspect
        usr <- c(xrange, yrange) + rep(c(-1, 1), 2) * 
          rep(d, c(2, 2))
        par(usr = usr)
      }
      on.exit(par(opar))
    }
    if (type != "n") {
      if (!as.polygon && resolution != 0) {
        pin <- par("pin")
        usr <- par("usr")
        resolution <- resolution * min(diff(usr)[-2]/pin/100)
        coord[c("x", "y")] <- mapthin(coord, resolution)
      }
      if (fill) 
        polygon(coord, col = col, ...)
      else lines(coord, col = col, type = type, ...)
      ##### check to see if the points are over the polygon
      return(1)
      test = in.poly(points, cbind(coord$x, coord$y))
    }
  }
  class(coord) = "map"
  value <- if (namesonly) 
    coord$names
  else coord
  if (plot) 
    invisible(value)
  else value
}

myMakePoly = function (xy, gonsize, keep) 
{
  x <- xy$x
  y <- xy$y
  n <- length(x)
  gonsize <- gonsize[-length(gonsize)]
  discard <- seq(length(x))[is.na(x)]
  if (length(discard) > 0) {
    dups = c(discard - 1, n)
    i = which(diff(c(0, dups)) > 2)
    discard <- c(discard, dups[i])
  }
  if (length(gonsize) > 0) 
    discard <- discard[-cumsum(gonsize)]
  if (length(discard) > 0) {
    x <- x[-discard]
    y <- y[-discard]
  }
  keep <- rep(keep, diff(c(0, seq(length(x))[is.na(x)], length(x))))
  closed.polygon(list(x = x[keep], y = y[keep]))
}
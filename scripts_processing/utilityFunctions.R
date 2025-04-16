
no.legend <- theme(legend.position="none")

geom_diagline <- function(linetype='solid',size=0.1,colour="grey20",...) {
    geom_abline(slope=1,intercept=0,linetype=linetype,colour=colour)
}

scientific_10 <- function(x) {
    xout <- gsub("1e", "10^{", format(x,scientific=T),fixed=TRUE)
    xout <- gsub("{-0", "{-", xout,fixed=TRUE)
    xout <- gsub("{+", "{", xout,fixed=TRUE)
    xout <- gsub("{0", "{", xout,fixed=TRUE)
    xout <- paste(xout,"}",sep="")
    return(parse(text=xout))
}

scale_x_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_x_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_y_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_y_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_loglog <- function(xname=waiver(), yname=waiver(), ...) {
    list(scale_x_log10nice(name=xname, ...),scale_y_log10nice(name=yname,...))
}

scale_x_log2nice <- function(name=waiver(),omag=seq(-6,6),...) {
    breaks2 <- 2^omag
    scale_x_log10(name,breaks=breaks2,
                  labels=parse(text=paste("2^{",omag,"}")),...)
}

scale_y_log2nice <- function(name=waiver(),omag=seq(-6,6),...) {
    breaks2 <- 2^omag
    scale_y_log10(name,breaks=breaks2,
                  labels=parse(text=paste("2^{",omag,"}")),...)
}

scale_loglog2 <- function(...) {
    list(scale_x_log2nice(...),scale_y_log2nice(...))
}

theme_density <- function(...) {
    list(scale_y_continuous(expand=c(0.01,0.01)),
         theme(panel.border=element_blank(),
               panel.background=element_blank(),
               axis.line=element_blank(),
               axis.text.y=element_blank(),
               axis.title.y=element_blank(),
               axis.ticks.y=element_blank()),...)
}

# No legend for plots
# Usage: ggplot(d, aes(...)) + no.legend
no.legend <- theme(legend.position="none")


logfinite <- function(x) {
  xl <- log(x)
  xl[!is.finite(xl)] <- NA
  xl
}

corfinite <- function(x,y=NULL,use='pairwise.complete.obs',method="pearson") {
    # wrapper to calculate correlation coefficients for only finite values
    # useful for correlations of log-transformed values
  if (!is.null(y)) {
    niceinds <- which(is.finite(x) & is.finite(y))
    res <- cor(x[niceinds],y[niceinds],method=method)
  } else {
    x[!is.finite(x)] <- NA
    res <- cor(x,use=use,method=method)

  }
  res
}

logcor <- function(x,y=NULL,use='pairwise.complete.obs',method="pearson") {
    # wrapper for correlations of log-transformed values
  x <- as.matrix(x)
  ly <- y
  if (!is.null(y)) {ly <- log(as.matrix(y))}
  corfinite(log(x),ly,use=use,method=method)
}

odds <- function(p) {
  p/(1-p)
}

odds2p <- function(o) {
  o/(1+o)
}

p2odds <- function(p) {
  res <- odds(p)
  res[!is.finite(res)] <- NA
  res
}

p2logodds <- function(p) {
  res = log(p2odds(p))
  res
}

logodds <- function(x) {
  # log odds, a shortcut
  res <- log(odds(x))
  res[!is.finite(res)] <- NA
  res
}

invlogodds <- function(x) {
  # inverse log odds: a = invlogodds(logodds(a))
  y <- exp(x)
  y/(1+y)
}

# Use a fixed window size for x. 
# n is the number of windows, or (if a list is given)
# the center of each window
# Window size is 2*range(x)/n
windowfunction_stratum = function(x, y, n, fn=mean, logx=F, minn=20, minwidth_fraction=NULL, na.rm=TRUE, ...) {
  
  if (length(n)>1) {
    # centers are given
    range_x = n
    n_windows = length(range_x)
    range_xin = range_x
  } else {
    n_windows = n
    range_xin = NULL
    range_x = NULL
  }
  
  xin = x
  if (logx) {
    xin = log(xin)
  }
  minx = min(xin, na.rm=na.rm)
  maxx = max(xin, na.rm=na.rm)
  
  if (is_null(range_x)) {
    range_x = seq(minx, maxx, length.out=n_windows)
  } else {
    if (logx) {
      range_x = log(range_x)
    }
  }
  # Calculate window size.
  if (is_null(minwidth_fraction)) {
    windowsize = 2*(maxx-minx)/n_windows
  } else {
    # use minwidth_fraction, which is a fraction of the range
    windowsize = minwidth_fraction*(maxx-minx)
  }
  
  d = tibble(x=xin, y=y)
  yvals = bind_rows(lapply(range_x, function(center_x){
    # Lower bound of window
    lower_x = max(minx, center_x-windowsize, na.rm=na.rm)
    # Upper bound of window
    upper_x = min(maxx, center_x+windowsize, na.rm=na.rm)
    # Filter to isolate values in window
    wind = d |> filter(x>=lower_x & x<=upper_x) |> select(y)
    ny = nrow(wind |> filter(!is.na(y)))
    # Compute value of function in window, or NA if fewer than minn values
    if (nrow(wind)<minn) {res = c(y=NA, n=nrow(wind))} else {res = c(y=fn(unlist(wind), na.rm=na.rm, ...), nx=nrow(wind), ny=ny)}
    res
  }))
  xout = range_x
  if (logx) {
    xout = exp(range_x)
  }
  if (!is_null(range_x)) {
    xout = range_xin
  }
  bind_cols(tibble(x=xout), yvals)
}

# Rolling function (e.g. mean)
windowfunction = function(d, var, window_var, by, fn=mean, ...) {
  # Isolate levels of "by" variable
  strata = unique(d[[by]])
  # Apply the strata and bind resulting rows
  res = bind_rows(lapply(strata, function(stratum){
    # Filter by stratum
    x = d |> filter(!!as.name(by)==stratum)
    # Calculate windowed function for this stratum
    wm = windowfunction_stratum(x[[window_var]], x[[var]], fn=fn, ...) |> 
      mutate(!!by:=stratum) |> rename(!!window_var:=x, !!var:=y, !!paste(window_var,"n",sep="_"):=nx, !!paste(var,"n",sep="_"):=ny)
    wm
  }))
  as_tibble(res)
}




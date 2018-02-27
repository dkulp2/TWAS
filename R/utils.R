# returns a string containing the first n elements of v
# vec2str(letters) => "a,b,c,d,e,..."
vec2str <- function(v, n=5, sep=',') {
  x <- v[1:min(length(v),n)]
  if (n < length(v)) {
    paste(c(x, "..."), collapse=sep)
  } else {
    paste(x, collapse=',')
  }
}

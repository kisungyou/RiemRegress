############# Auxiliary R Functions #########################
# (01) rmreg_stack3d     : stack data into 3d cube
# (02) rmreg_equivariant : make an equivariant embedding stacked as rows
# (03) rmreg_10foldidx   : create 10-fold labels
# (04) rmreg_cubesubset  : remove anomaly in subsetting cube by third dimension
# (05) rmreg_cube2list   : as it says

# (01) rmreg_stack3d     : stack data into 3d cube
#' @keywords internal
#' @noRd
rmreg_stack3d <- function(riemdata){
  msize = riemdata$size
  ndata = length(riemdata$data)
  
  matdata = array(0,c(msize[1], msize[2], ndata))
  for (i in 1:ndata){
    matdata[,,i] = (riemdata$data[[i]])
  }
  return(matdata)
}

# (02) rmreg_equivariant : make an equivariant embedding as rows
#' @keywords internal
#' @noRd
rmreg_equivariant <- function(input){
  mfddata = rmreg_stack3d(input)
  mfdname = input$name
  output  = cpp_equivariant(mfddata, mfdname)
  return(output)
}

# (03) rmreg_10foldidx   : create 10-fold labels
#' @keywords internal
#' @noRd
rmreg_10foldidx <- function(ndata){
  tmpfd <- cut(seq(1:ndata),breaks=10,labels=FALSE)
  folds <- tmpfd[sample(1:ndata,ndata,replace=FALSE)]
  return(folds)
}

# (04) rmreg_cubesubset  : remove anomaly in subsetting cube by third dimension
#' @keywords internal
#' @noRd
rmreg_cubesubset <- function(datacube, subsetid){
  nrow = dim(datacube)[1]
  ncol = dim(datacube)[2]
  nslice = length(subsetid)
  
  output = array(0,c(nrow,ncol,nslice))
  for (i in 1:nslice){
    tmp = datacube[,,subsetid[i]]
    if (is.vector(tmp)){
      tmp = as.matrix(tmp)
    }
    output[,,i]=tmp
  }
  return(output)
}

# (05) rmreg_cube2list   : as it says
#' @keywords internal
#' @noRd
rmreg_cube2list <- function(datacube){
  n = dim(datacube)[3]
  output = list()
  for (i in 1:n){
    tmp = datacube[,,i]
    if (is.vector(tmp)){
      tmp = as.matrix(tmp)
    }
    output[[i]] = tmp
  }
  return(output)
}
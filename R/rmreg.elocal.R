#' Extrinsic Local Regression on Manifold-Valued Data
#' 
#' 
#' @export
rmreg.elocal <- function(y, x=NULL, h=seq(from=0.1,to=2,by=0.1), parallel=TRUE,
                         predict=FALSE, x.predict=NULL){
  #-------------------------------------------------------
  # PREPROCESSING
  # 1. response must be of 'riemdata' class with number of data and mfd name and size
  if ((class(y))!="riemdata"){
    stop("* rmreg.elocal : the response 'y' must be of 'riemdata' class. Use 'RiemBase::riemfactory' first to manage your data.")
  }
  ndata   = length(y$data)
  mfdsize = y$size; dnrow = mfdsize[1]; dncol = mfdsize[2];
  mfdname = y$name
  # 2. x : (ndata-by-p)
  if ((length(x)==0)&&(is.null(x))){
    x = matrix((1:ndata))
  } else if (is.vector(x)){
    x = matrix(x)
  }
  if ((nrow(x)!=ndata)||(!is.matrix(x))){
    stop("* rmreg.elocal : 'x' must be either a vector or a matrix with corresponding number of observations.")
  }
  # 3. h
  if (any(h<=0)){
    stop("* rmreg.elocal : kernel bandwidths must be nonnegative.")
  }
  if ((ndata<=10)&&(length(h)>1)){
    stop("* rmreg.elocal : with the small number of data, please use a single plug-in for a bandwidth parameter.")
  }
  # 4. Equivariant Embedding of the Data stacked as ROWS
  y.cube  <- rmreg_stack3d(y)
  y.equiv <- rmreg_equivariant(y)
  
  # 5. parallel setup
  if (parallel){
    nCore = parallel::detectCores()
  } else {
    nCore = 1
  }

  #-------------------------------------------------------
  # MAIN 1 : CV or Single 'h'
  if (length(h)>1){ ## case of CV to get optimal value
    if (nCore==1){ 
      CVscore = rep(0,rep(length(h)))
      for (i in 1:length(h)){
        CVscore[i] = rmreg.elocal.cvsingle(y.cube, x, h[i], mfdname)
        print(h[i])
      }
    } else {
      cl = makeCluster(nCore)
      registerDoParallel(cl)
      # let's care about parallel implementation issue.
      itforeach=NULL
      CVscore  = foreach (itforeach=1:length(h), .combine=cbind) %dopar% {
        rmreg.elocal.cvsingle(y.cube, x, h[itforeach], mfdname)
      }
      stopCluster(cl)
    }
    CVscore = as.vector(CVscore)
    idxmin  = which.min(CVscore)
    h.optimal  = h[idxmin]
  } else {
    h.optimal  = h # simply, that one.
  }

  # MAIN 2 : Compute Fitted Values into a Cube
  yres.fit = elocal_fit(h.optimal,y.cube,x,mfdname)
  yfit.list = rmreg_cube2list(yres.fit$fitted) # as a list
  yfit.mse  = yres.fit$mse
  
  # MAIN 3 : Prediction
  if (predict==TRUE){
    if ((length(x.predict)==0)&&(is.null(x.predict))){
      stop("* rmreg.elocal : given predict option 'TRUE' with no 'x.predict', it's sufficient to check fitted values only.")
    }
    if (is.vector(x.predict)){
      x.predict = as.matrix(x.predict)
    }
    if (ncol(x.predict)!=ncol(x)){
      stop("* rmreg.elocal : the number of independent variables for 'x.predict' does not match to our original setting.")
    }
    ypredict = elocal_predict(h.optimal,y.cube,x,x.predict,mfdname) # fitted values only
  }

  #-------------------------------------------------------
  # RETURN RESULTS
  output = list()
  output$fit = RiemBase::riemfactory(yfit.list,name=mfdname)
  output$mse = yfit.mse
  if (length(h)>1){
    output$h.optimal = h.optimal
  }
  if (predict==TRUE){
    output$predict = RiemBase::riemfactory(rmreg_cube2list(ypredict),name=mfdname)
  }
  return(output)
}


#' @keywords internal
#' @noRd
rmreg.elocal.cvsingle <- function(y, x, h.val, mfdname){ # y now accepts cube only
  ndata = nrow(x)
  cvidx = rmreg_10foldidx(ndata)
  score = rep(0,10)
  for (i in 1:10){
    testidx  = which(cvidx==i)
    trainidx = setdiff(1:ndata,testidx)
    
    ytrain = rmreg_cubesubset(y,trainidx)
    ytest  = rmreg_cubesubset(y,testidx)
    
    xtrain = x[trainidx,]
    xtest  = x[testidx,]
    if (is.vector(xtrain)){
      xtrain = matrix(xtrain)
    }
    if (is.vector(xtest)){
      xtest = matrix(xtest)
    }
    
    score[i] = tryCatch({elocal_fit_cv(h.val, ytrain, ytest, xtrain, xtest, mfdname)},
                        warning = function(w){return(NA)},
                        error = function(e){return(NA)})
  }
  return(as.matrix(mean(score)))
}

# 
# 
# pole.exp <- function(t, Ft){
#   datmat = array(0,c(3,length(t)))
#   datmat[1,] = t
#   datmat[2,] = Ft
#   datmat[3,] = 1
# 
#   outmat = array(0,c(3,length(t)))
# 
#   x = c(0,0,1)
#   for (i in 1:length(t)){
#     u = as.vector(datmat[,i])
#     d = u-(x*(sum(x*u)))
# 
#     nrm_td = norm(matrix(d),type="F")
#     if (nrm_td < 1e-10){
#       y = x
#     } else {
#       y = cos(nrm_td)*x + ((sin(nrm_td))/nrm_td)*d
#     }
#     outmat[,i] = y
#   }
#   return(outmat)
# }
# ntime=100
# a = 1
# b = 1
# error = 0.1
# t      = seq(from=-0.75,to=0.75,length.out = ntime)
# Ft.pop = sin(2*pi*t)+rnorm(ntime,sd=error) # true curve
# 
# curve.error = RiemBase::riemfactory(pole.exp(t,Ft.pop), name="sphere")
# curve.fit   = rmreg.elocal(curve.error)
# 
# 
# 
#   rmreg.elocal <- function(y, x=NULL, h=seq(from=0.1,to=2,by=0.1), parallel=TRUE,
#                            predict=FALSE, x.predict=NULL){
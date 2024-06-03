#
#       PanDA: Pan-omics Discriminant Analysis;
#
#
#       [eigvector, eigvalue] = PanDA(fea, labels, Dim, meu)
#
#              Input:
#                  fea     - cell array containing different omics data matrix (X_v).
#                            Each row vector of X_v is a sample vector.
#                  labels  - Colunm vector contaiing label information for each
#                            data point.
#                  Dim     - no. of extracted eigenvectors.
#                  meu     - regularization parameter
#
#              output:
#                  eigvector - Each column is an embedding function, for a new
#                              sample vector (row vector) x,  y = x*eigvector
#                              will be the embedding result of x.
#                  eigvalue  - The sorted eigvalue of the eigen-problem.
#
#    Reference:
#
#
#    version 0.1.0 --Oct/2022
#
#    Written by M. Aminu (muhammadaminu47 AT gmail.com)
#

#' @param data
#'
#' @param labels
#' @param Dim
#' @param meu
#'
#' @return
#' @export
#'
#' @importFrom geigen geigen

PanDA <- function (data,labels,Dim, meu){
  if (missing(data)){
    stop("data not provided")
  }

  if (missing(labels)){
    stop("class labels needed but not provided")
  }

  if (missing(Dim)){
    Dim <- 2
    message("number of extracted components not provided, default = 2 components will be extracted")
  }

  # user specified parameters
  lambda = 1e-1 # regularizer to handle singularity problem of the total scatter matrix

  if (missing(meu)){
    meu <- 0.2
    message("regularization parameter not provided, default = 0.2 will be used")
  }

  # Calculate the number of samples in each class from all views
  labels <- as.numeric(labels)

  C = length(unique(labels))
  V = length(data)
  n <- matrix(0, nrow = C, ncol = V)

  for (i in 1:C) {
    for (j in 1:V) {
      n[i,j] <- sum(labels==i)
    }

  }

  # Compute the mean of samples in each class for all views
  u <- matrix(list(), 1, V)
  for (j in 1:V) {
    m = 0
    view = data[[j]]
    ii = nrow(view)
    uu <- matrix(0, nrow = ii, ncol = C)
    for (i in 1:C) {
      mOld <- m
      m = m + n[i,j]
      mOld = mOld + 1
      uu[,i] = rowMeans(view[,c(mOld:m)])
    }
    u[[j]] <- uu
  }

  # Compute the between class scatter matrix
  eta = 0
  for (i in 1:C) {
    eta = eta + (sum(n[i,]))^2
  }
  delta = eta/((sum(n[,1]))-1)^2

  BB = vector()
  for (j in 1:V) {
    bigB = vector()
    view1 = data[[j]]
    jj = nrow(view1)
    for (r in 1:V) {
      view2 <- data[[r]]
      rr = nrow(view2)
      bb = matrix(0, nrow = jj, ncol = rr)
      mm = matrix(0, nrow = jj, ncol = 1)
      nn = matrix(0, nrow = rr, ncol = 1)
      uij = u[[j]]
      uir = u[[r]]
      for (i in 1:C) {
        bb = bb + (n[i,j]*n[i,r]/(sum(n[,1])-1)^2*uij[,i] %*% t(uir[,i]))
        mm = mm + n[i,j]*uij[,i]
        nn = nn + n[i,r]*uir[,i]
      }
      tempB = bb - delta*mm %*% t(nn)
      bigB = cbind(bigB,tempB)
    }
    BB = rbind(BB,bigB)
  }

  # Compute the correspondence matrix M
  M <- matrix(list(), V, V)
  P <- matrix(list(), 1, V)
  for (i in 1:V) {
    view <- data[[i]]
    P[[i]] <- ginv(view)
  }

  MM = vector()
  for (i in 1:V) {
    tempM = vector()
    for (j in 1:V) {
      M = vector()
      if (i == j) {
        M <- 2*(V-1)*t(P[[i]]) %*% P[[i]]
      }
      else {
        M <- -2*t(P[[i]]) %*% P[[j]]
      }
      tempM = cbind(tempM,M)
    }
    MM = rbind(MM,tempM)
  }

  # Compute the multi-view total scatter matrix C
  CC = vector()
  for (i in 1:V) {
    tempC = vector()
    view1 <- data[[i]]
    for (j in 1:V) {
      C = vector()
      view2 <- data[[j]]
      if (i == j) {
        C <- 1/sum(n[,1])*view1 %*% t(view2)
      }
      else {
        C <- matrix(0, nrow = nrow(view1 %*% t(view2)), ncol = ncol(view1 %*% t(view2)))
      }
      tempC = cbind(tempC,C)
    }
    CC = rbind(CC,tempC)
  }

  # Compute the projection matrix W
  # BB = pmax(BB, t(BB))
  # MM = pmax(MM, t(MM))
  # CC = pmax(CC, t(CC))

  DD = BB - meu*MM
  CC = CC + lambda*diag(nrow(CC));
  jeig <- geigen(DD, CC, symmetric = TRUE)
  ind <- order(-jeig[["values"]])
  eigvector <- jeig[["vectors"]][,ind[c(1:Dim)]]

  # EXtract PanDA components
  dataType <- names(data)
  compNames <- paste(rep("DC",Dim),1:Dim,sep = " ")
  PanDAComponents = list()
  projMatrices = list()
  d = 0
  for (k in 1:V) {
    projMatrices$comp <- eigvector[c(d+1:nrow(data[[dataType[k]]])),]
    names(projMatrices)[k] <- paste("W",dataType[k],sep = "")
    PanDAComponents$comp <- t(data[[dataType[k]]]) %*% projMatrices[[k]]
    colnames(PanDAComponents[[k]]) <- compNames
    names(PanDAComponents)[k] <- paste(dataType[k],"Components",sep = "")
    d = d + nrow(data[[dataType[k]]])
  }

  PanDAModel <- list(projMatrices,PanDAComponents)
  names(PanDAModel) <- c("projMatrices","PanDAComponents")

  return(PanDAModel)

}


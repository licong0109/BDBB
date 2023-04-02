######------------------------------mulitiple bicluster model-------------------------------########
####MCE蒙特卡洛误差,链内方差#####
#' Title
#'
#' @param a input data matrix
#' @param iteration iteration Number of initial iterations
#' @param burn_in Number of burnings
#' @param bic_row_labels Labels for rows in the patterns found, all with initial values of 0
#' @param number_bicluster Number of preset patterns
#'
#' @return
#' @export
#'
#' @examples
#'
#'

multiple_GSB <- function(a,iteration,burn_in,bic_row_labels,number_bicluster){
  out = list()
  for (c in 1:number_bicluster) {
    GSB_single_object <- GSB_SINGLE(a = a,iteration = iteration,burn_in = burn_in, bic_row_labels = bic_row_labels)

    assign(paste( c,"bicluster", sep="_"),GSB_single_object)
    out[[paste( c,"bicluster", sep="_")]] <- get(paste( c,"bicluster", sep="_"))
    save(out,file =  "./out")
    if( sum(GSB_single_object$row_labels)==0|sum(GSB_single_object$col_labels)==0 ){
      break;
    }else{
      bic_row_labels <- bic_row_labels + GSB_single_object$row_labels # 更新bicluster 标签

    }



  }
  return(out)

}
MCE <- function(samp,v){
  b <- 0
  mean_b <- 0
  k <- length(samp)/v
  for (i in 1:k) {
    # b[i] <- samp[((i-1)*v+1):(i*v)]

    mean_b[i] <- mean(samp[((i-1)*v+1):(i*v)]) # 组内均值


  }
  mce <- sqrt(var(mean_b)/(k-1))
  return(mce)

}
# debug(MCE)
# tt <- MCE(multiple_GSB_object$`1_bicluster`$bic_pro[501:1000],v=50)



GSB_SINGLE <- function(a,iteration,burn_in,bic_row_labels){
  # initialize row and column lables

  n <- dim(a)[1]
  m<- dim(a)[2]
  matrix_pro_row <- matrix(0,n,1)# construct a matrix to save the value of problity of each row of iteration
  matrix_pro_col <- matrix(0,m,1)
  matrix_labels_row <- matrix(0,n,1)# construct a matrix to save the value of labels of each row of iteration
  matrix_labels_col <- matrix(0,m,1)
  k <- 1
  row_theta <- c(2,3)
  col_theta <- c(4,10)
  pro_row <- rbeta(1,row_theta[1],row_theta[2])
  pro_col <- rbeta(1,col_theta[1],col_theta[2])
  row_labels <- c(rbinom(n = dim(a)[1], size = 1, prob = pro_row))
  col_labels <- c(rbinom(n = dim(a)[2], size = 1, prob = pro_col))
  # initialize bicluster and background parameters
  bic_theta <- c(4,10)
  bic_shap_para.1_chain <- bic_theta[1]
  bic_shap_para.2_chain <- bic_theta[2]
  aa <- 5
  bb <- 10
  cc <- 9
  dd <- 3
  bag_theta <- c(9,3)
  index_row <- c(which(bic_row_labels==0))
  # 去除位点中已分配给其它bicluster的位点
  row_labels[c(which(bic_row_labels!=0))] <- 0
  likelihood <-0

  # ini_pro_bic <- rbeta(1,bic_theta[1],bic_theta[2])
  # ini_pro_bag <- rbeta(1,bag_theta[1],bag_theta[2])
  arrary_pro_row <- 0
  arrary_pro_col <- 0
  repeat{
    for (i in index_row) {

      # label of i equal 1
      bic_row1 <- i # row index of bicluster
      bic_col1 <- c(which(col_labels==1))    # column index of bicluster


      f.1 <- sum(log(dbeta(a[i,bic_col1],bic_theta[1],bic_theta[2])/dbeta(a[i,bic_col1],bag_theta[1],bag_theta[2])))



      # v <- table(row_labels[-i]) # statistic number of each comphonent of all the row labeles exclude i
      # v_ibar <- v[[2]]# the number of 1 in the row lables exclude
      v_ibar <- length(which(row_labels[-i]==1))
      f.2 <- log((v_ibar + row_theta[1])/(n-v_ibar+row_theta[2]-1))

      odd_i <- f.1+f.2
      update_pro_rowi <- 1/(1+exp(-odd_i))
      if(is.na(update_pro_rowi)){
        update_pro_rowi <- 0.5
      }
      arrary_pro_row[i] <- update_pro_rowi
      # matrix_pro_row[i,k] <- update_pro_rowi# save the problity of label of i-th row
      #mask
      # if(is.na(rbinom(1,1,update_pro_rowi))){browser()}
      # if(bic_row_labels[i]== 0){
      #   row_labels[i] <- rbinom(1,1,update_pro_rowi)
      #
      # }else{
      #   row_labels[i] <- 0
      # }

      i_label <- rbinom(1,1,update_pro_rowi)
      row_labels[i] <- i_label

    }
    matrix_pro_row <- cbind(matrix_pro_row,arrary_pro_row)
    matrix_labels_row <- cbind(matrix_labels_row, row_labels)
    if(length(which(row_labels==1))==0){
      iteration <- k
      burn_in <- 1
      break;
    }
    update_bic_row <- c(which(row_labels==1))
    bic_col <- c(which(col_labels==1))
    update_bag_row <- c(which(row_labels==0))
    bag_col <- c(which(col_labels==0))


    means_bic <- mean(c(a[update_bic_row,bic_col]))
    vars_bic <- var(c(a[update_bic_row,bic_col]))

    com_bic <- means_bic * (1 - means_bic) / vars_bic - 1 # 贝塔分布的参数距估计
    update_a <- means_bic * com_bic
    update_b <- (1 - means_bic) * com_bic
    bic_theta[1] <- update_a
    bic_theta[2] <-update_b
    # aa[k+1] <- update_a
    # bb[k+1] <- update_b
    means_bag <- mean(c(a[update_bag_row,],a[update_bic_row,bag_col]))
    vars_bag <- var(c(a[update_bag_row,],a[update_bic_row,bag_col]))
    com_bag <- means_bag * (1 - means_bag) / vars_bag - 1
    update_c <- means_bag * com_bag
    update_d <- (1 - means_bag) * com_bag
    bag_theta[1] <- update_c
    bag_theta[2] <- update_d
    # cc[k+1] <- update_c
    # dd[k+1] <- update_d



    for (j in 1:m) {

      bic_col2 <- j
      bic_row2 <- c(which(row_labels==1))


      g.1 <- sum(log(dbeta(a[bic_row2,j],bic_theta[1],bic_theta[2])/dbeta(a[bic_row2,j],bag_theta[1],bag_theta[2])))




      # w <- table(col_labels[-j])
      # w_jbar <- w[[2]]
      w_jbar <- length(which(col_labels[-j]==1))
      g.2 <- log((w_jbar+col_theta[1])/(m-w_jbar+col_theta[2]-1))
      odd_j <- g.1 + g.2
      update_pro_colj <- 1/(1+exp(-odd_j))
      if(is.na(update_pro_colj)){
        update_pro_colj <- 0.5
      }
      arrary_pro_col[j] <- update_pro_colj
      # matrix_pro_col[j,k] <- update_pro_colj
      j_label <- rbinom(1,1,update_pro_colj)
      col_labels[j] <- j_label


    }
    matrix_pro_col <- cbind(matrix_pro_col,arrary_pro_col)
    matrix_labels_col <- cbind(matrix_labels_col, col_labels)
    # if(length(which(col_labels==1))==0){
    #   break;
    # }

    bic_row <- update_bic_row
    update_bic_col <- c(which(col_labels==1))
    bag_row <- update_bag_row
    update_bag_col <- c(which(col_labels==0))


    means_bic <- mean(c(a[bic_row,update_bic_col]))
    vars_bic <- var(c(a[bic_row,update_bic_col]))

    com_bic <- means_bic * (1 - means_bic) / vars_bic - 1 # 贝塔分布的参数距估计
    update_a <- means_bic * com_bic
    update_b <- (1 - means_bic) * com_bic
    bic_theta[1] <- update_a
    bic_theta[2] <-update_b
    aa[k+1] <- update_a
    bb[k+1] <- update_b
    means_bag <- mean(c(a[bag_row,],a[bic_row,update_bag_col]))
    vars_bag <- var(c(a[bag_row,],a[bic_row,update_bag_col]))
    com_bag <- means_bag * (1 - means_bag) / vars_bag - 1
    update_c <- means_bag * com_bag
    update_d <- (1 - means_bag) * com_bag
    bag_theta[1] <- update_c
    bag_theta[2] <- update_d
    cc[k+1] <- update_c
    dd[k+1] <- update_d

    # 存储bicluster的两个形状参数，用于收敛判断
    bic_shap_para.1_chain[k+1] <- update_a
    bic_shap_para.2_chain[k+1] <- update_b
    # if(is.na(bag_theta[1])){browser()}
    # 计算每次迭代的loglikelihood
    # bicluster的 likelihood
    index.bic.row <- c(which(row_labels==1))
    index.bic.col <- c(which(col_labels==1))
    l_bi <- log(dbeta(a[index.bic.row,index.bic.col],update_a,update_b))

    l_bi[which(l_bi == -Inf)] <- -.Machine$double.xmin

    if(length(l_bi)==0){
      l_bi <- 0
    }
    likelihood.1 <- sum(l_bi)
    #background 的likelihood
    index.all.row <- c(1:n)
    index.all.col <- c(1:m)
    index.bag.row <- setdiff(index.all.row,index.bic.row)
    index.bag.col <- setdiff(index.all.col,index.bic.col)
    l_ba <- log(dbeta(a[index.bag.row,index.bag.col],update_c,update_d))


    l_ba[which(l_ba == -Inf)] <- -.Machine$double.xmin
    if(length(l_ba)==0){
      l_ba <- 0
    }
    likelihood.2 <- sum(l_ba)
    likelihood[k] <- sum(likelihood.1,likelihood.2)

    #收敛判断
    # if(k==iteration){
    #
    #   if(MCE(samp = bic_shap_para.1_chain[(iteration-500+1):iteration],v = 50)< 0.01|iteration==3000)  #为了得到高表达的bicluster放大v减小阈值
    #   {
    #     break
    #   }else{
    #     iteration <- iteration+500
    #   }
    #
    #
    # }
    if(k==iteration){
      if(iteration<500){
        break
      }else{
        # browser()

        #likelihood的方差在一定范围内
        # if(var(bic_pro_chain[(iteration-500+1):iteration])< 0.01|iteration==2000)
        if(MCE(samp = bic_shap_para.1_chain[(iteration-500+1):iteration],v = 50)< 1e-5|iteration==2000)  #为了得到高表达的bicluster放大v减小阈值
        {
          break
        }else{
          iteration <- iteration+500
        }

      }



    }
    k <- k+1




  }
  # 如果运算过程中，得到的行的标签row_labels全部为0，说明此时没有比较尖锐的贝塔分布，对于剩下的数据算法认为其是相对较平坦的贝塔分布，令其行标签全部为0，作为背景。
  if(length(which(row_labels==1))==0){
    update_row_labels <-  rep(0,n)
    update_col_labels <-  rep(0,m)
    # update_col_labels <-  apply(matrix_labels_col[,(burn_in+1):iteration], 1, sum)/(iteration-burn_in)

  }else{
    update_row_labels <-  apply(matrix_labels_row[,(burn_in+1):iteration], 1, sum)/(iteration-burn_in)
    update_col_labels <-  apply(matrix_labels_col[,(burn_in+1):iteration], 1, sum)/(iteration-burn_in)

  }


  # update_row_labels <-  apply(matrix_labels_row[,(burn_in+1):iteration], 1, sum)/(iteration-burn_in)
  # update_col_labels <-  apply(matrix_labels_col[,(burn_in+1):iteration], 1, sum)/(iteration-burn_in)
  # update_row_labels <- apply(matrix(mapply(function(x,y) x*y, matrix_labels_row[,(burn_in+1):iteration],matrix_pro_row[,(burn_in+1):iteration]),nr=n),1,sum)/(iteration-burn_in)
  # update_col_labels <- apply(matrix(mapply(function(x,y) x*y, matrix_labels_col[,(burn_in+1):iteration],matrix_pro_col[,(burn_in+1):iteration]),nr=n),1,sum)/(iteration-burn_in)
  update_row_labels[which(update_row_labels >= 0.75)] <- 1
  update_row_labels[which(update_row_labels < 0.75)] <- 0
  update_col_labels[which(update_col_labels >= 0.75)] <- 1
  update_col_labels[which(update_col_labels < 0.75)] <- 0


  out <- list()
  out[['row_labels']] <- update_row_labels
  out[['col_labels']] <- update_col_labels
  out[['matrix_row_labels']] <- matrix_labels_row
  out[['matrix_col_labels']] <- matrix_labels_col
  out[["bic_a"]] <- aa
  out[["bic_b"]] <- bb
  out[["bag_a"]] <- cc
  out[["bag_b"]] <- dd
  out[["log_likelihood"]] <- likelihood
  out[["bic.shapparameter.1"]] <- bic_shap_para.1_chain
  out[["bic.shapparameter.2"]] <- bic_shap_para.2_chain

  return(out)





}







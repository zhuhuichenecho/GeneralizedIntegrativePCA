
EPCAJIVEMissbio = function(data, rankj, ranka, D, family, tol,max.iter = 500,n.size = NULL,lambda = 0){
  
  
  K = rankj
  dimention = c(0,cumsum(D))
  
  X_joint = do.call(cbind,data)
  N = nrow(X_joint)
  P = ncol(X_joint)
  
  ######## Initialize starting value ########
  
  Familyx = rep(family,times = D)
  
  dg = function(theta,n.size) {
    
    
    dgfun = NULL
    for(d in 1:length(D)){
      if(family[d]=="binomial")
        dgfun[d] = paste("dg",family[d],"(theta[",dimention[d]+1,":",dimention[d+1],"],n.size[",
                         dimention[d]+1,":",dimention[d+1],"])",sep = "")
      else
      dgfun[d] = paste("dg",family[d],"(theta[",dimention[d]+1,":",dimention[d+1],"])",sep = "")
      
    }
    
    eval(parse(text = paste("c(",paste(dgfun,collapse = ","),")",sep = "")))
    
  }
  
  MU = function(theta,n.size){
    
    MUfun = NULL
    for(d in 1:length(D)){
      
      if(family[d]=="binomial")
        MUfun[d] = paste("MU",family[d],"(theta[",dimention[d]+1,":",dimention[d+1],"],n.size[",
                         dimention[d]+1,":",dimention[d+1],"])",sep = "")
      else
      MUfun[d] = paste("MU",family[d],"(theta[",dimention[d]+1,":",dimention[d+1],"])",sep = "")
      
    }
    
    eval(parse(text = paste("c(",paste(MUfun,collapse = ","),")",sep = "")))
    
  }
  
  loglikfun = NULL
  
  for(d in 1:length(D)){
    if(family[d]=="binomial"){
      
      loglikfun = c(loglikfun,paste("loglik",family[d],"(data[[",d,"]],Theta[,",dimention[d]+1,":",
                                    dimention[d+1],"],n.size[,",dimention[d]+1,":",dimention[d+1],"])",sep = ""))
      
    }else{
      
      loglikfun = c(loglikfun,paste("loglik",family[d],"(data[[",d,"]],Theta[,",dimention[d]+1,":",
                                    dimention[d+1],"])",sep = ""))
      
      
    }
    
    
  }
  
  
  loglik_pre = -Inf
  loglik_cur = -Inf
  l = 1
  # set.seed(1234)
  U.joint = cbind(rep(1,N),matrix(runif(N*K,-1,1),N))
  # set.seed(1234)
  V.joint = matrix(runif((K+1)*P,-1,1),(K+1))
  
  U.ind = list()
  V.ind = list()
  
  for(d in 1:length(D)){
    # set.seed(1234)
    U.ind[[d]] = matrix(runif(sum(rowSums(is.na(data[[d]]))<D[d])*ranka[d],-1,1),ncol = ranka[d])
    # set.seed(1234)
    V.ind[[d]] = matrix(runif(ranka[d]*D[d],-1,1),nrow = ranka[d])
    
  }
  n.size = n.size
  loglik.sep = NULL
  
  Theta.joint = U.joint%*%V.joint
  Theta.inds = matrix(0,nrow(X_joint),ncol(X_joint))
  Theta = Theta.joint+Theta.inds
  eval(parse(text = paste("loglik_cur = ",paste(loglikfun,collapse = "+"),sep = "")))
  
  loglik = NULL
  
  loglik = c(loglik,loglik_cur)
  
  
  while ((abs(loglik_cur-loglik_pre)>tol||l==1)&l<=max.iter){
    
    # print(l)
    V_pre.joint = V.joint
    loglik_pre = loglik_cur
    
    if(rankj == 0){
      
      V_new.joint = do.call(cbind,sapply(1:ncol(X_joint),function(i){
        # print(i)
        if (Familyx[i]=="binomial1"){
          
          V_update = V.joint[,i]
          V_pre = V.joint[,i]+1
          l1=1
          while(max(abs(V_pre-V_update))>1e-2|l1<=2){
            
            # print(l1)
            
            mu = MUbinomial1(Theta.joint[,i]+Theta.inds[,i])
            DG = dgbinomial1(mu)
            V_pre = matrix(V_update,nrow = K+1)
            index = is.na(X_joint[,i])==0
            
            U_trans = t(U.joint[index,]*sqrt(1/DG)[index])
            X_trans = sqrt(1/DG)[index]*(matrix(U.joint[index,],ncol = K+1)%*%
                                           matrix(V_pre,nrow = K+1)+matrix((X_joint[index,i]-mu[index])*DG[index],ncol = 1))
            
            
            V_update = solve(U_trans%*%t(U_trans))%*%U_trans%*%X_trans
            
            # print(V_update)
            
            Theta.joint[,i] = U.joint%*%matrix(V_update,nrow = K+1)
            
            l1 = l1+1
          }
          #  print(j)
          matrix(V_update,nrow = K+1)
          
          
        }else if (Familyx[i]=="binomial"){
          
          glm(cbind(X_joint[,i],n.size[,i]-X_joint[,i])~-1+U.joint+offset(Theta.inds[,i]),family=Familyx[i])$coeff
          
        }else if(Familyx[i]=="binomial2"){
          
          glm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]),family="binomial")$coeff
          
        }else if(Familyx[i]=="gaussian"){
          
          lm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]))$coef
          
        }else{
          
          glm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]),family=Familyx[i])$coeff
          
        }
        
        
      },simplify = F))
      V.joint = V_new.joint
      
      Theta.joint = matrix(V_new.joint,nrow = nrow(X_joint),ncol = ncol(X_joint),byrow = T)
      
      project.matrix = diag(nrow(X_joint))
      
      Theta = Theta.joint + Theta.inds
      
      # print(angle(t(svd(scale(Theta.joint,center = T,scale = F),nu = K,nv = K)$v),t(V)))
      
      eval(parse(text = paste("loglik_cur = ",paste(loglikfun,collapse = "+"),sep = "")))
      # loglik = c(loglik,loglik_cur)
      
    }else{
      
      V_new.joint = do.call(cbind,sapply(1:ncol(X_joint),function(i){
        # print(i)
        if (Familyx[i]=="binomial1"){
          
          V_update = V.joint[,i]
          V_pre = V.joint[,i]+1
          l1=1
          while(max(abs(V_pre-V_update))>1e-2|l1<=2){
            
            # print(l1)
            
            mu = MUbinomial1(Theta.joint[,i]+Theta.inds[,i])
            DG = dgbinomial1(mu)
            V_pre = matrix(V_update,nrow = K+1)
            index = is.na(X_joint[,i])==0
            
            U_trans = t(U.joint[index,]*sqrt(1/DG)[index])
            X_trans = sqrt(1/DG)[index]*(matrix(U.joint[index,],ncol = K+1)%*%
                                           matrix(V_pre,nrow = K+1)+matrix((X_joint[index,i]-mu[index])*DG[index],ncol = 1))
            
            
            V_update = solve(U_trans%*%t(U_trans))%*%U_trans%*%X_trans
            
            # print(V_update)
            
            Theta.joint[,i] = U.joint%*%matrix(V_update,nrow = K+1)
            
            l1 = l1+1
          }
          #  print(j)
          matrix(V_update,nrow = K+1)
          
          
        }else if (Familyx[i]=="binomial"){
          
          glm(cbind(X_joint[,i],n.size[,i]-X_joint[,i])~-1+U.joint+offset(Theta.inds[,i]),family=Familyx[i])$coeff
          
        }else if(Familyx[i]=="binomial2"){
          
          glm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]),family="binomial")$coeff
          
        }else if(Familyx[i]=="gaussian"){
          
          lm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]))$coef
          
        }else{
          
          glm(X_joint[,i]~-1+U.joint+offset(Theta.inds[,i]),family=Familyx[i])$coeff
          
        }
        
        
      },simplify = F))
      Theta.joint = U.joint%*%V_new.joint
      Theta = Theta.joint+Theta.inds
      eval(parse(text = paste("loglik_cur = ",paste(loglikfun,collapse = "+"),sep = "")))
      
      # loglik = c(loglik,loglik_cur)
      
      U_new.joint = do.call(rbind,sapply(1:nrow(X_joint),function(j){
        # print(j)
        
        if(sum(family=="gaussian")==length(D)){
          
          
          c(1,lm(X_joint[j,]~-1+t(matrix(V_new.joint[-1,],nrow = K))+offset(Theta.inds[j,]+V_new.joint[1,]))$coef)
          
        }else if(sum(family=="binomial")==length(D)){
          
          c(1,glm(cbind(X_joint[j,],n.size[j,]-X_joint[j,])~-1+t(matrix(V_new.joint[-1,],
                                                                         nrow = K))+
                    offset(Theta.inds[j,]+V_new.joint[1,]),family = "binomial")$coef)
          
        }else{
    
          
          
          if(eval(parse(text = paste(sapply(1:length(D),function(d)  
            paste("sum(is.na(X_joint[j,(dimention[",d,"]+1):dimention[",d+1,"]]))>=D[",d,"]",sep = "")),collapse = "|")))){
            
            
            index = which(sapply(1:length(D),function(d)  
              eval(parse(text = paste("sum(is.na(X_joint[j,(dimention[",d,"]+1):dimention[",d+1,"]]))>=D[",d,"]",sep = ""))))!=TRUE)
            
            
            if(family[index]=="binomial")
              c(1,glm(cbind(X_joint[j,(dimention[index]+1):dimention[index+1]],n.size[j,(dimention[index]+1):dimention[index+1]]
                        -X_joint[j,(dimention[index]+1):dimention[index+1]])~-1+t(matrix(V_new.joint[-1,(dimention[index]+1):dimention[index+1]],
                                                                                         nrow = K))+
                    offset(Theta.inds[j,(dimention[index]+1):dimention[index+1]]+
                             V_new.joint[1,(dimention[index]+1):dimention[index+1]]),family = "binomial")$coef)
            else
              c(1,glm(X_joint[j,(dimention[index]+1):dimention[index+1]]~-1+t(matrix(V_new.joint[-1,(dimention[index]+1):dimention[index+1]],
                                                                                         nrow = K))+
                            offset(Theta.inds[j,(dimention[index]+1):dimention[index+1]]+
                                     V_new.joint[1,(dimention[index]+1):dimention[index+1]]),family = family[index])$coef)
            
            
            
          }
          
          else{
            
            U_update = c(1,U.joint[j,-1])
            U_pre = c(1,U.joint[j,-1]+1)
            l1=1
            
            while(max(abs(U_pre[-1]-U_update[-1]))>1e-1|l1<=10){
              
              # print(l1)
              
              
              mu = MU(Theta.joint[j,]+Theta.inds[j,],n.size[j,])
              DG = dg(mu,n.size[j,])
              U_pre = matrix(U_update,ncol = K+1)
              index = is.na(X_joint[j,])==0
              
              V_trans = t(t(matrix(V_new.joint[-1,index],nrow = K))*sqrt(1/DG)[index])
              X_trans = sqrt(1/DG)[index]*(t(matrix(V_new.joint[-1,index],K))%*%
                                             t(matrix(U_pre[-1],ncol = K))+matrix((X_joint[j,index]-mu[index])*DG[index],ncol = 1))
              
              
              U_update = c(1,t(solve(V_trans%*%t(V_trans))%*%(V_trans%*%X_trans)-
                                 (ncol(U.joint)-1)*lambda*t(matrix(U_pre[-1],ncol = K))))
              # print(U_update)
              
              Theta.joint[j,] = U_update%*%V_new.joint
              loglikgaussian(X_joint[j,1:(P/2)],Theta.joint[j,1:(P/2)])+loglikbinomial1(X_joint[j,(P/2+1):P],Theta.joint[j,(P/2+1):P])
              
              l1 = l1+1
            }
            #  print(j)
            U_update
            
          }
            
            

          
        }
        
        
      },simplify = F))
      
      U.joint = U_new.joint
      V.joint = V_new.joint
      
      Theta.joint = U.joint%*%V.joint
      
      Theta = Theta.joint + Theta.inds
      
      # print(angle(t(svd(scale(Theta.joint,center = T,scale = F),nu = K,nv = K)$v),t(V)))
      
      eval(parse(text = paste("loglik_cur = ",paste(loglikfun,collapse = "+"),sep = "")))
      # loglik = c(loglik,loglik_cur)
    }
    
    
    
    ####################################
    ####### individual structure #######
    ####################################
    
    Theta.inds.temp = list()
    
    for (d in 1:length(data)){
      
      Theta.inds.temp[[d]] = matrix(NA,nrow(data[[d]]),ncol(data[[d]]))
      
      X.ind = data[[d]][rowSums(is.na(data[[d]]))<D[d],]
      N.ind = nrow(X.ind)
      K.ind = ranka[d]
      if(family[d] == "binomial")
      n.size.ind = n.size[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]]
      
      
      if(ranka[d]==0){
        
        Theta.inds.temp[[d]] = matrix(0,N,ncol(data[[d]]))
        U.ind[[d]] = matrix(0,nrow(X.ind),K.ind)
        V.ind[[d]] = matrix(0,K.ind,ncol(X.ind))
        
      }else{
        
        Theta.nna.ind = Theta.joint[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]]
        
        Theta.ind = Theta.inds[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]]
        
        if(family[d]=="binomial1"){
          
          V_new.ind = do.call(cbind,sapply(1:ncol(X.ind),function(i){
            # print(i)
            
            V_update = V.ind[[d]][,i]
            V_pre = V.ind[[d]][,i]+1
            l1=1
            while(max(abs(V_update-V_pre))>5e-1|l1<=3){
              
              # print(l1)
              
              mu = MUbinomial1(Theta.nna.ind[,i]+as.vector(U.ind[[d]]%*%matrix(V_update,nrow = ranka[[d]])))
              DG = dgbinomial1(mu)
              V_pre = matrix(V_update,nrow = ranka[[d]])
              
              U_trans = t(matrix(U.ind[[d]],ncol = ranka[d])*sqrt(1/DG))
              X_trans = sqrt(1/DG)*(matrix(U.ind[[d]],ncol = ranka[d])%*%
                                      matrix(V_pre,nrow = ranka[d])+matrix((X.ind[,i]-mu)*DG,ncol = 1))
              
              
              V_update = t(solve(U_trans%*%t(U_trans))%*%U_trans%*%X_trans)
              
              # print(V_update)
              
              Theta.ind[,i] = U.ind[[d]]%*%matrix(V_update,nrow = ranka[[d]])
              
              # print(loglikbinomial1(X.ind[,i],Theta.nna.ind[,i]+U.ind[[d]]%*%matrix(V_update,nrow = ranka[[d]])))
              # loglik.ind = c(loglik.ind,loglikbinomial1(X.ind[,i],Theta.nna.ind[,i]+U.ind[[d]]%*%matrix(V_update,nrow = ranka[[d]])))
              l1 = l1+1
            }
            #  print(j)
            matrix(V_update,nrow = ranka[d])
            
          },simplify = F))
          
          Theta.ind = U.ind[[d]]%*%V_new.ind
          
          V.ind[[d]] = V_new.ind
          
          U_new.ind = do.call(rbind,sapply(1:nrow(X.ind),function(j){
            # print(j)
            U_update = U.ind[[d]][j,]
            U_pre = U_update+4
            l1=1
            while(max(abs(U_pre-U_update))>1e-1){
              
              # print(l1)
              
              
              mu = MUbinomial1(Theta.nna.ind[j,]+Theta.ind[j,])
              DG = dgbinomial1(mu)
              U_pre = matrix(U_update,ncol = ranka[[d]])
              
              V_trans = t(t(matrix(V_new.ind,nrow = ranka[[d]]))*sqrt(1/DG))
              X_trans = sqrt(1/DG)*(t(matrix(V_new.ind,ranka[d]))%*%
                                      t(matrix(U_pre,ncol = ranka[d]))+matrix((X.ind[j,]-mu)*DG,ncol = 1))
              
              
              U_update = t(solve(V_trans%*%t(V_trans))%*%V_trans%*%X_trans)
              
              # print(U_update)
              
              Theta.ind[j,] = U_update%*%V_new.ind
              
              # print(loglikbinomial1(X.ind[j,],Theta.nna.ind[j,]+U_update%*%V_new.ind))
              
              l1 = l1+1
            }
            #  print(j)
            U_update
            
          },simplify = F))
          
          U.ind[[d]] = U_new.ind
          
          Theta.ind = U.ind[[d]]%*%V.ind[[d]]
          
          
          
        }else{
          
          V_new.ind = sapply(1:ncol(X.ind),function(i){
            # print(i)
            if(family[d]=="binomial"){
              
              glm(cbind(X.ind[,i],n.size.ind[,i]-X.ind[,i])~-1+U.ind[[d]]+offset(Theta.nna.ind[,i]),family = family[d])$coef
              
            }else if(family[d]=="binomial2") {
              glm(X.ind[,i]~-1+U.ind[[d]]+offset(Theta.nna.ind[,i]),family = "binomial")$coef
              
            }else if(family[d]=="gaussian"){
              
              lm(X.ind[,i]~-1+U.ind[[d]]+offset(Theta.nna.ind[,i]))$coef
              
            }else{
              
              glm(X.ind[,i]~-1+U.ind[[d]]+offset(Theta.nna.ind[,i]),family = family[d])$coef
              
            }
            
            
          })
          
          V_new.ind = matrix(V_new.ind,nrow = K.ind)
          
          U_new.ind = do.call(rbind,sapply(1:nrow(X.ind),function(j){
            # print(j)
            if (family[d]=="binomial"){
              
              glm(cbind(X.ind[j,],n.size.ind[j,]-X.ind[j,])~-1+t(V_new.ind)+offset(Theta.nna.ind[j,]),family = family[d])$coef
              
            }else if(family[d]=="binomial2"){
              
              glm(X.ind[j,]~-1+t(V_new.ind)+offset(Theta.nna.ind[j,]),family = "binomial")$coef
              
            }else if(family[d]=="gaussian"){
              
              lm(X.ind[j,]~-1+t(V_new.ind)+offset(Theta.nna.ind[j,]))$coef
              
            } else{
              
              glm(X.ind[j,]~-1+t(V_new.ind)+offset(Theta.nna.ind[j,]),family = family[d])$coef
              
            }
            
          },simplify = F)) 
          
          
        }
        
        
        
        U_new.ind = matrix(U_new.ind,ncol = K.ind)
        
        Theta.ind = U_new.ind%*%V_new.ind
        U.ind[[d]] = svd(Theta.ind,nu = ranka[d],nv = ranka[d])$u
        V.ind[[d]] = diag(svd(Theta.ind,nu = ranka[d],nv = ranka[d])$d,ranka[d])%*%
          t(svd(Theta.ind,nu = ranka[d],nv = ranka[d])$v)
        
        
        
        Theta.inds.temp[[d]][rowSums(is.na(data[[d]]))<D[d],] = Theta.ind
        
        rm(Theta.nna.ind,V_new.ind,U_new.ind,X.ind,N.ind,K.ind,Theta.ind)
        
      }
      
      
      
      
    }
    Theta.inds = do.call(cbind,Theta.inds.temp)
    
    
    # U.ind.orth = list()
    # 
    # for(d in 1:length(D)){
    #   
    #   U.ind.orth.temp = matrix(NA,nrow(data[[d]]),ranka[d])
    #   U.ind.orth.temp[rowSums(is.na(data[[d]]))<D[d],] = U.joint[rowSums(is.na(data[[d]]))<D[d],]
    #   U.ind.orth[[d]] = diag(N)-U.ind.orth.temp%*%t(U.ind.orth.temp)
    #   
    # }
    # 
    # for(d in 1:length(D)){
    #   
    #   U.ind.temp = matrix(NA,nrow(data[[d]]),ranka[d])
    #   U.ind.temp[rowSums(is.na(data[[d]]))<D[d],] = U.ind[[d]]
    #   
    #   for(i in (1:(length(D)-1))){
    #     
    #     U.ind.temp[eval(parse(text = paste(sapply(c((1:length(D))[-d][i],d),function(d) paste("rowSums(is.na(data[[",d,
    #                                                                                           "]]))<D[",d,"]",sep = "")),
    #                                        collapse = "&"))),] = 
    #       U.ind.orth[[(1:length(D))[-d][i]]][eval(parse(text = paste(sapply(c((1:length(D))[-d][i],d),function(d) paste("rowSums(is.na(data[[",d,
    #                                                                                                                     "]]))<D[",d,"]",sep = "")),
    #                                                                  collapse = "&"))),
    #                                          eval(parse(text = paste(sapply(c((1:length(D))[-d][i],d),function(d) paste("rowSums(is.na(data[[",d,
    #                                                                                                                     "]]))<D[",d,"]",sep = "")),
    #                                                                  collapse = "&")))]%*%
    #       U.ind.temp[eval(parse(text = paste(sapply(c((1:length(D))[-d][i],d),function(d) paste("rowSums(is.na(data[[",d,
    #                                                                                             "]]))<D[",d,"]",sep = "")),
    #                                          collapse = "&"))),]
    #     
    #   }
    #   
    #   U.ind[[d]] = U.ind.temp[rowSums(is.na(data[[d]]))<D[d],]
    #   Theta.inds.temp[[d]] = matrix(NA,nrow(data[[d]]),ncol(data[[d]]))
    #   Theta.inds.temp[[d]][rowSums(is.na(data[[d]]))<D[d],] = U.ind[[d]]%*%V.ind[[d]]
    #   
    # }
    # Theta.inds.new = do.call(cbind,Theta.inds.temp)
    # 
    # Theta.joint.new = Theta.joint+Theta.inds-Theta.inds.new
    # Theta.joint.new[which(is.na(Theta.joint.new),arr.ind = T)] = Theta.joint[which(is.na(Theta.joint.new),arr.ind = T)]
    # Theta.joint = Theta.joint.new
    # Theta.inds = Theta.inds.new
    
    U.joint[,-1] = svd(scale(Theta.joint,scale = F),nu = K,nv = K)$u
    project.matrix = diag(N)-U.joint%*%solve(t(U.joint)%*%U.joint)%*%t(U.joint)
    
    for(d in 1:length(D)){
      
      if(ranka[d]==0){
        
        Theta.inds.temp[[d]] = matrix(0,nrow(data[[d]]),ncol(data[[d]]))
        
      }else{
        
        Theta.inds.temp[[d]] = matrix(NA,nrow(data[[d]]),ncol(data[[d]]))
        U.ind.temp = scale(project.matrix[rowSums(is.na(data[[d]]))<D[d],rowSums(is.na(data[[d]]))<D[d]]%*%
                             U.ind[[d]],scale = F)
        U.ind[[d]] = svd(U.ind.temp%*%V.ind[[d]],nv = ranka[d],nu = ranka[d])$u%*%
          diag(svd(U.ind.temp%*%V.ind[[d]],nv = ranka[d],nu = ranka[d])$d[ranka[d]],ranka[d])
        V.ind[[d]] = t(svd(U.ind.temp%*%V.ind[[d]],nv = ranka[d],nu = ranka[d])$v)
        Theta.inds.temp[[d]][rowSums(is.na(data[[d]]))<D[d],] = U.ind[[d]]%*%V.ind[[d]]
        
      }
      
    }
    Theta.inds.new = do.call(cbind,Theta.inds.temp)
    
    
    Theta.joint.new = Theta.joint+Theta.inds-Theta.inds.new
    Theta.joint.new[which(is.na(Theta.joint.new),arr.ind = T)] = Theta.joint[which(is.na(Theta.joint.new),arr.ind = T)]
    V.joint[1,] = colMeans(Theta.joint.new)
    if(rankj>1){
      
      U.joint[,-1] = svd(scale(Theta.joint.new,scale = F),nv = K,nu = K)$u%*%
        diag(svd(scale(Theta.joint.new,scale = F),nv = K,nu = K)$d,K)
      V.joint[-1,] = t(svd(scale(Theta.joint.new,scale = F),nv = K,nu = K)$v)
      
    }
    
    Theta.joint = Theta.joint.new
    Theta.inds = Theta.inds.new
    # if(sum(ranka)>0){
    #   
    #   Vind = eval(parse(text = paste("as.matrix(bdiag(",paste("V.ind[[",1:length(D),"]]",sep = "",
    #                                                           collapse = ","),"))",sep = "")))
    #   project.matrix = diag(P)-t(Vind)%*%Vind
    #   
    #   Theta.joint.proj = U.joint%*%rbind(V.joint[1,],V.joint[-1,]%*%project.matrix)
    # }else{
    #   
    #   Theta.joint.proj = Theta.joint
    # }
    # 
    # 
    # 
    # if(K>0){
    #   
    #   V.joint[1,] = colMeans(Theta.joint.proj)
    #   V.joint[-1,] = t(svd(scale(Theta.joint.proj,scale = F),nu = K,nv = K)$v)
    #   U.joint[,-1] = svd(scale(Theta.joint.proj,scale = F),nu = K,nv = K)$u%*%diag(svd(scale(Theta.joint.proj,
    #                                                                                          scale = T),
    #                                                                                    nu = K,nv = K)$d,rankj)
    # }
    # 
    # 
    # 
    # Theta.ind.new = Theta.joint-Theta.joint.proj+Theta.inds
    # 
    # for (d in 1:length(D)){
    #   if(ranka[d]>0){
    #     
    #     U.ind[[d]] = svd(Theta.ind.new[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]],
    #                      nu = ranka[d],nv = ranka[d])$u%*%
    #       diag(svd(Theta.ind.new[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]],
    #                nu = ranka[d],nv = ranka[d])$d,ranka[d])
    #     V.ind[[d]] = t(svd(Theta.ind.new[rowSums(is.na(data[[d]]))<D[d],(dimention[d]+1):dimention[d+1]],
    #                        nu = ranka[d],nv = ranka[d])$v)
    #     
    #     
    #   }else{
    #     
    #     U.ind[[d]] = matrix(0,nrow(na.omit(data[[d]])),ranka[d])
    #     V.ind[[d]] = matrix(0,ncol(na.omit(data[[d]])),ranka[d])
    #   }
    #   
    # }
    # 
    # Theta.joint = Theta.joint.proj
    # Theta.inds = Theta.ind.new
    
    
    Theta = Theta.joint + Theta.inds
    # Theta[which(is.na(Theta.inds),arr.ind = T)] = Theta.joint[which(is.na(Theta.inds),arr.ind = T)]
    eval(parse(text = paste("loglik_cur = ",paste(loglikfun,collapse = "+"),sep = "")))
    loglik = c(loglik,loglik_cur)
    
    l = l+1
    
    # print(angle(t(svd(scale(Theta.joint,scale = F),nu = K,nv = K)$v),t(V)))
  }
  
  
  
  singj =  sapply(1:length(D),function(d){
    
    round(svd(t(t(Theta.joint)-V.joint[1,])[,(dimention[d]+1):(dimention[d+1])])$d[1:rankj],2)
    
  },simplify = F)
  
  mu = sapply(1:length(D),function(d){
    
    colMeans(Theta.joint)
    
  })
  
  
  singa = sapply(1:length(D),function(d){
    
    round(svd(na.omit(Theta.inds[,(dimention[d]+1):dimention[d+1]]))$d[1:ranka[d]],2)
    
  },simplify = F)
  
  if(rankj>0){
    
    U.joint[,-1] = svd(scale(Theta.joint,scale = F),nu = K,nv = K)$u
    V.joint[-1,] = t(svd(scale(Theta.joint,scale = F),nu = K,nv = K)$v)
    
  }
  
  
  V.ind = sapply(1:length(D),function(d){
    
    if(ranka[d]>0){
      t(svd(na.omit(Theta.inds[,(dimention[d]+1):dimention[d+1]]),nu = ranka[d],nv = ranka[d])$v)
    }else{
      0
    }
    
  },simplify = F)
  
  return(list(Theta.joint = Theta.joint,Theta.inds = Theta.inds,mu = V.joint[1,],V.joint = V.joint, 
              V.ind = V.ind,U.joint = U.joint,U.ind = U.ind,
              loglik = loglik_cur,singj = singj,singa = singa,l =l))
  
}






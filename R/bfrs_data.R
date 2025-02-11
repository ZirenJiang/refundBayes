brfs_data=function(formula,data,family,func_comp,intercept,cens=NULL,func_parameter=NULL){

  trans.mat=list()

  result_code_data=paste0("data{ \n")
  result_code_para=paste0("parameters{ \n")
  result_code_transpara=paste0("transformed parameters { \n   real lprior = 0;\n")
  result=list()

  X_scalar=brfs_scalar_pred(formula$scalar_var,data)

  if(intercept){
    #X_scalar=cbind(rep(1,NROW(X_scalar)),X_scalar)
    result_code_para=paste0(result_code_para, "   //Linear predictor intercept  \n")
    result_code_para=paste0(result_code_para, "   real eta_0;\n")
  }


  result[["N_num"]]=NROW(data)
  result_code_data=paste0(result_code_data, "   //Total number of observations  \n")
  result_code_data=paste0(result_code_data,"   int<lower=1> N_num;\n")

  if(family=="Cox"){
    args=list()
    y=data[[formula$y_var]]
    min_y <- min(y, na.rm = TRUE)
    max_y <- max(y, na.rm = TRUE)
    diff_y <- max_y - min_y
    lower_knot <- max(min_y - diff_y / 50, 0)
    upper_knot <- max_y + diff_y / 50
    args$Boundary.knots <- c(lower_knot, upper_knot)
    args$x <- y

    Mbasis=splines2:::mSpline(x=args$x,Boundary.knots = args$Boundary.knots,df=5L,intercept=TRUE)
    Ibasis=splines2:::iSpline(x=args$x,Boundary.knots = args$Boundary.knots,df=5L)




    result[["real_inter"]]=0
    result_code_data=paste0(result_code_data,"   real real_inter;\n")

    result[["cens"]]=data[[cens]]
    result_code_data=paste0(result_code_data, "   //Indicator variable for censoring  \n")
    result_code_data=paste0(result_code_data,"   array[N_num] int<lower=-1,upper=2> cens;\n")



    result[["L_num"]]=NCOL(Mbasis)
    result_code_data=paste0(result_code_data, "   //Number of baseline hazard function basis  \n")
    result_code_data=paste0(result_code_data,"   int L_num;\n")

    result[["con_sbhaz"]]= rep(1,NCOL(Mbasis))
    result_code_data=paste0(result_code_data,"   vector<lower=0>[L_num] con_sbhaz;\n")


    result[["Mbasis"]]=Mbasis
    result_code_data=paste0(result_code_data, "   //Baseline hazard function basis  \n")
    result_code_data=paste0(result_code_data,"   matrix[N_num, L_num] Mbasis;\n")

    result[["Ibasis"]]=Ibasis
    result_code_data=paste0(result_code_data, "   //Cumulative baseline hazard function basis  \n")
    result_code_data=paste0(result_code_data,"   matrix[N_num, L_num] Ibasis;\n")

    result[["Y"]]=data[[formula$y_var]]
    result_code_data=paste0(result_code_data, "   //Time-to-event outcome variable  \n")
    result_code_data=paste0(result_code_data,"   real Y[N_num];\n")

    result_code_para=paste0(result_code_para, "   //Baseline hazard spline coefficient   \n")
    result_code_para=paste0(result_code_para, "   simplex[L_num] c;\n")




  }
  if(family=="gaussian"){

    result_code_para=paste0(result_code_para, "   real<lower=0> sigma;\n")
    result[["Y"]]=data[[formula$y_var]]
    result_code_data=paste0(result_code_data, "   //Outcome variable   \n")
    result_code_data=paste0(result_code_data,"   real Y[N_num];\n")

    result[["plocation"]]=mean(data[[formula$y_var]])
    result[["pscale"]]=10
    result_code_data=paste0(result_code_data,"   real plocation;\n")
    result_code_data=paste0(result_code_data,"   real pscale;\n")
  }
  if(family=="binomial"){
    #result_code_para=paste0(result_code_para, "   int<lower=0> sigma;\n")

    result[["Y"]]=data[[formula$y_var]]
    result_code_data=paste0(result_code_data, "   //Outcome variable   \n")
    result_code_data=paste0(result_code_data,"   int Y[N_num];\n")

    result[["plocation"]]=log((mean(data[[formula$y_var]]))/(1-mean(data[[formula$y_var]])))
    result[["pscale"]]=10
    result_code_data=paste0(result_code_data,"   real plocation;\n")
    result_code_data=paste0(result_code_data,"   real pscale;\n")

  }

  if(family=="functional"){
    #result_code_para=paste0(result_code_para, "   int<lower=0> sigma;\n")
    result[["T"]]=dim(data[[formula$y_var]])[2]
    result_code_data=paste0(result_code_data,"   //Number of observed time points\n")
    result_code_data=paste0(result_code_data,"   int<lower=1> T;\n")

    result[["Y"]]=data[[formula$y_var]]
    result_code_data=paste0(result_code_data,"   //Functional outcome for FoSR\n")
    result_code_data=paste0(result_code_data,"   matrix[N_num, T] Y;\n")

    freq.fpca=refund::fpca.face(unclass(data[[formula$y_var]]))

    data_temp=data.frame(inx=1:NROW(data))
    data_temp[["yindex.vec"]]=matrix(rep(1:dim(data[[formula$y_var]])[2],NROW(data)),nrow = NROW(data))
    obj_exp=paste0("object=s(yindex.vec, bs = \"",func_parameter[["type"]],"\", k =",func_parameter[["df"]],")")

    eval(parse(text = obj_exp))
    dk = ExtractData(object,data_temp,NULL)
    crspline.cons = smooth.construct(object,dk$data,dk$knots)
    beta.base = crspline.cons$X

    beta.S = crspline.cons$S[[1]]
    maXX = norm(beta.base,type="I")^2
    maS = norm(beta.S)/maXX
    beta.S = beta.S / maS

    result[["ncol_phi"]]=dim(beta.base)[2]
    result_code_data=paste0(result_code_data,"   //Number of FPCA eigenfunctions\n")
    result_code_data=paste0(result_code_data,"   int ncol_phi; \n")

    result[["phi"]]=t(beta.base)
    result_code_data=paste0(result_code_data,"   //Matrix of FPCA eigenfunctions\n")
    result_code_data=paste0(result_code_data,"   matrix[ncol_phi, T] phi; \n")

    # result[["K"]]=NCOL(X_scalar)
    # result_code_data=paste0(result_code_data,"   int<lower=1> K; \n")
    #
    # result[["X"]]=X_scalar
    # result_code_data=paste0(result_code_data,"   matrix[N, K] phi; \n")

    result[["ncol_psi_mat"]]=ifelse(is.null(dim(freq.fpca$efunctions)[2]),1,dim(freq.fpca$efunctions)[2])
    result_code_data=paste0(result_code_data,"   //Number of spline basis\n")
    result_code_data=paste0(result_code_data,"   int ncol_psi_mat; \n")

    result[["psi_mat"]]=t(freq.fpca$efunctions)
    result_code_data=paste0(result_code_data,"   //Matrix of spline basis\n")
    result_code_data=paste0(result_code_data,"   matrix[ncol_psi_mat, T] psi_mat;  \n")

    result[["S_beta"]]=beta.S
    result_code_data=paste0(result_code_data,"   //Penalty matrix\n")
    result_code_data=paste0(result_code_data,"   matrix[ncol_phi, ncol_phi] S_beta;  \n")


    result[["K_num"]]=NCOL(X_scalar)
    result_code_data=paste0(result_code_data,"   //Number of scalar predictors\n")
    result_code_data=paste0(result_code_data,"   int<lower=1> K_num;\n")

    if(!is.null(X_scalar)){
      result[["Z_mat"]]=X_scalar
      result_code_data=paste0(result_code_data,"   //Matrix of scalar predictors\n")
      result_code_data=paste0(result_code_data,"   matrix[N_num,K_num] Z_mat;\n")

      result_code_para=paste0(result_code_para,"   //FPCA scores \n")
      result_code_para=paste0(result_code_para, "   matrix[N_num,ncol_psi_mat] zxi;\n")
      result_code_para=paste0(result_code_para, "   real<lower=0> sigma_epis;\n")

      result_code_para=paste0(result_code_para,"   //Smoothing parameters\n")
      result_code_para=paste0(result_code_para, "   vector<lower=0>[K_num] sigma_b;\n")


      result_code_para=paste0(result_code_para, "   vector<lower=0>[ncol_psi_mat] sigma_k;\n")


      if(FALSE){
      #if(NCOL(X_scalar)==1){
        result_code_para=paste0(result_code_para, "   row_vector[ncol_phi] b;\n")
      }else{
        result_code_para=paste0(result_code_para, "   //Functional effect spline coefficients\n")
        result_code_para=paste0(result_code_para, "   matrix[K_num,ncol_phi] b;\n")
      }


    }



  }else{

    ### For SoFR and Functional Cox

    result[["K_num"]]=NCOL(X_scalar)
    result_code_data=paste0(result_code_data, "   //Number of scalar predictors   \n")
    result_code_data=paste0(result_code_data,"   int<lower=1> K_num;\n")

    if(!is.null(X_scalar)){
      result[["Z_mat"]]=X_scalar
      result_code_data=paste0(result_code_data, "   //Matrix of scalar predictors   \n")
      result_code_data=paste0(result_code_data,"   matrix[N_num,K_num] Z_mat;\n")
      result_code_para=paste0(result_code_para, "   vector[K_num] gamma;\n")

      #result[["K_c"]]=NCOL(Z_matcalar)-1
      #result_code_data=paste0(result_code_data,"   int<lower=1> K_c;\n")
    }
  }







  if(length(formula$func_var)>0){
    if(sum(func_comp)>0){

    }


    for(i in 1:length(func_comp)){
      if(func_comp[i]){

      }

      if(!func_comp[i]){
        func_data=brfs_func_pred_obs(term=formula$func_var[[i]],data = data)

        trans.mat[[i]]=func_data$trans.mat

        result[[paste0("Kr_",i)]]=NCOL(func_data$Xr)
        result_code_data=paste0(result_code_data,"   int<lower=1> ",paste0("Kr_",i),";\n")

        result[[paste0("X_mat_r_",i)]]=func_data$Xr
        result_code_data=paste0(result_code_data,"   matrix[N_num, ",paste0("Kr_",i),"] ",paste0("X_mat_r_",i),";\n")
        result_code_para=paste0(result_code_para, "   vector[",paste0("Kr_",i),"] ",paste0("zbr_",i),";\n")
        result_code_para=paste0(result_code_para, "   real<lower=0>",paste0("sigmabr_",i),";\n")

        result_code_transpara=paste0(result_code_transpara, "   vector[",paste0("Kr_",i),"] ",paste0("br_",i),";\n")
        result_code_transpara=paste0(result_code_transpara, "   ",paste0("br_",i)," = ",paste0("sigmabr_",i)," * ",paste0("zbr_",i),";\n")


        result[[paste0("Kf_",i)]]=NCOL(func_data$Xf)
        result_code_data=paste0(result_code_data,"   int<lower=1> ",paste0("Kf_",i),";\n")

        result[[paste0("X_mat_f_",i)]]=func_data$Xf
        result_code_data=paste0(result_code_data,"   matrix[N_num, ",paste0("Kf_",i),"] ",paste0("X_mat_f_",i),";\n")
        result_code_para=paste0(result_code_para, "   vector[",paste0("Kf_",i),"] ",paste0("bf_",i),";\n")

      }
    }
  }






  result_code_data=paste0(result_code_data,"}")
  result_code_para=paste0(result_code_para,"}")
  result_code_transpara=paste0(result_code_transpara,"}")

  return(list(standata=result,
              stancode_data=result_code_data,
              stancode_para=result_code_para,
              stancode_transpara=result_code_transpara,
              trans.mat=trans.mat,
              X_scalar=X_scalar))
}


brfs_scalar_pred=function(term, data){
  if(length(term)==0){
    return(NULL)
  }else{
    X_use=matrix(1,nrow =NROW(data), ncol = 1)
    for(i in 1:length(term)){

      if(length(term[[i]])>1){
        #X_itm=data[[term[[i]][[1]]]]
        term_now = term[[i]][[1]]
        for(k in 2:length(term[[i]])){
          term_now = paste0(term_now, "*", term[[i]][[k]])
        }
        X_itm=model.matrix(as.formula(paste0("~",term_now)),data = data)

      }else{
        X_itm=model.matrix(as.formula(paste0("~",term[[i]])),data = data)
        term_now=term[[i]]
      }

      colsave = colnames(X_itm)
      X_itm = as.matrix(X_itm[,-1])
      colnames(X_itm) = colsave[-1]
      # if(i == 1){
      #   term_use=term_now
      # }else{
      #   term_use=paste0(term_use,"+",term_now)
      # }

      X_use=cbind(X_use,X_itm)
    }
    #X_use=model.matrix(as.formula(paste0("~-1+",term_use)),data = data)
    X_use=X_use[,-1]

    good_inx=c()
    for(j in 1:ncol(X_use)){
      if(sum(abs(X_use[,j]))==0){
        #X_use=X_use[,-j]
      }else{
        good_inx=c(good_inx,j)
      }
    }
    X_use=X_use[,good_inx]
    return(X_use=X_use)
  }

}

brfs_func_pred_obs=function(term, data){
  if(term[[1]]=="s"){
    obj=format(term)
    #sm = mgcv:::smoothCon(obj,data = data,absorb.cons = TRUE,diagonal.penalty = TRUE)
    eval(parse(text = paste0("sm = mgcv:::smoothCon(",obj,",data = data,absorb.cons = TRUE,diagonal.penalty = TRUE)")))
    #eval(expr = )
    re = mgcv::smooth2random(sm[[1]], names(data), type = 2)
    Xr=re$rand$Xr
    Xf=re$Xf
    trans.mat=re$trans.U

  }else{

  }


  return(list(Xr=Xr, Xf=Xf,trans.mat=trans.mat))
}
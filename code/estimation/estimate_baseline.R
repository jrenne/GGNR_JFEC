# Optional re-estimation. The published and archived estimates are read-only.
source("code/inference/parameterization.R")
input<-Sys.getenv("GGNR_START_FILE","estimates/start_values.rds")
output<-Sys.getenv("GGNR_OUTPUT_FILE","estimates/reestimated.rds")
if(normalizePath(output,mustWork=FALSE)%in%normalizePath(c("estimates/baseline_r2.rds","estimates/start_values.rds"),mustWork=FALSE))
 stop("Choose a new output file; archived estimates are read-only.")
budget<-as.integer(Sys.getenv("GGNR_EVAL_BUDGET","2000"))
method<-match.arg(Sys.getenv("GGNR_OPTIMIZER","nlminb"),c("nlminb","Nelder-Mead","bobyqa"))
seconds<-as.numeric(Sys.getenv("GGNR_ESTIMATION_SECONDS","180"))
liquid<-tolower(Sys.getenv("GGNR_LIQUID_TIPS_ONLY","false"))%in%c("1","true")
data<-load_paper_data();short<-data$effective_federal_funds_rate
omit_start<-Sys.getenv("GGNR_OMIT_START","");omit_end<-Sys.getenv("GGNR_OMIT_END","")
if(nzchar(omit_start)&&nzchar(omit_end)){
 omit<-data$dates>=as.Date(omit_start)&data$dates<=as.Date(omit_end)
 for(n in c("macro","yields_n","yields_r","surv_infexp","surv_gdpexp","surv_tbexp"))data[[n]][omit,]<-NaN
 data$core_inflation[omit]<-NaN;short[omit]<-NaN
 data$any_survey_release[omit]<-data$any_nonmonthly_release[omit]<-FALSE
}
if(liquid){
 d<-data$dates
 keep<-d>=as.Date("2004-01-01")&!(d>=as.Date("2008-01-01")&d<as.Date("2010-01-01"))&!(d>=as.Date("2020-03-01")&d<as.Date("2021-03-01"))
 data$yields_r[!keep,]<-NaN
}
problem<-make_inference_problem(input,data,short,liquid)
theta<-problem$theta;nn<-names(theta);n<-length(theta)
C<-problem$estimate$differentiation_basis
if(is.null(C))C<-diag(n) else C<-C[match(nn,rownames(C)),,drop=FALSE]
B<-t(chol(tcrossprod(C)));dimnames(B)<-list(nn,nn)
B<-sweep(B,2,pmax(1,apply(abs(B),2,max)/.2),'/')
best<-problem$evaluate_fit(theta);best_theta<-theta;initial<-best$fval
calls<-0L;started<-proc.time()["elapsed"]
expand<-function(z)setNames(theta+drop(B%*%z),nn)
objective<-function(z){
 calls<<-calls+1L
 if(calls>budget||proc.time()["elapsed"]-started>seconds)stop("Explicit estimation budget reached.")
 v<-expand(z)
 f<-tryCatch(problem$evaluate_fit(v),error=function(e)NULL)
 if(is.null(f)||!is.finite(f$fval))return(1e12)
 if(f$fval>best$fval){best<<-f;best_theta<<-v}
 -f$fval
}
optimizer<-tryCatch({
 if(method=="nlminb")nlminb(rep(0,n),objective,control=list(eval.max=budget,iter.max=200,rel.tol=1e-12,x.tol=1e-8,step.max=.5))
 else if(method=="Nelder-Mead")optim(rep(0,n),objective,method="Nelder-Mead",control=list(maxit=budget,reltol=1e-10))
 else {
  if(!requireNamespace("nloptr",quietly=TRUE))stop("Install nloptr to use BOBYQA.")
  z<-nloptr::nloptr(rep(0,n),objective,opts=list(algorithm="NLOPT_LN_BOBYQA",maxeval=budget,maxtime=seconds,xtol_rel=1e-9,ftol_abs=1e-9))
  list(convergence=z$status,message=z$message)
 }
},error=function(e)list(convergence=98L,message=conditionMessage(e)))
result<-problem$estimate;result$working_parameters<-best_theta;result$parameters<-best$objects$pars
result$initial_loglik<-initial;result$final_loglik<-best$fval
result$exact_targets<-best$objects$targets;result$exact_moments<-best$objects$exact_moments
result$free_parameter_count<-n;result$liquid_tips_only<-liquid
result$differentiation_basis<-B;result$optimizer<-optimizer
result$stage<-"Optional re-estimation; verify convergence and Hessian before use"
result$timestamp<-Sys.time();dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
saveRDS(result,output)
message(sprintf("Saved %s; %d free coordinates; %d evaluations; %.1f seconds. Optimizer code %s is not a convergence certificate.",output,n,calls,proc.time()["elapsed"]-started,optimizer$convergence))

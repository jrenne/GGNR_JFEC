# Optional refinement of the 30 free coordinates, preserving the published point.
saved<-Sys.getenv(c("GGNR_START_FILE","GGNR_OUTPUT_FILE","GGNR_EVAL_BUDGET","GGNR_ESTIMATION_SECONDS","GGNR_OPTIMIZER"),unset=NA)
Sys.setenv(GGNR_START_FILE=Sys.getenv("GGNR_START_FILE","estimates/baseline_r2.rds"),
 GGNR_OUTPUT_FILE=Sys.getenv("GGNR_REFINE_OUTPUT","estimates/refined.rds"),
 GGNR_EVAL_BUDGET=Sys.getenv("GGNR_REFINE_MAXEVAL","2000"),
 GGNR_ESTIMATION_SECONDS=Sys.getenv("GGNR_REFINE_SECONDS","180"),GGNR_OPTIMIZER="nlminb")
source("code/estimation/estimate_baseline.R",local=new.env())
Sys.unsetenv(names(saved));keep<-!is.na(saved)
if(any(keep))do.call(Sys.setenv,as.list(saved[keep]))

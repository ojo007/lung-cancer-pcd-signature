# Script 49 — TCGA-LUAD out-of-fold survival calibration
# Uses the established 47 genes and reproduces Script 46 repeat 1,
# with an additional training-fold Cox baseline hazard to obtain
# 1- and 3-year survival probabilities. NO test-outcome refitting.
# Outputs are exploratory INTERNAL validation; NOT external calibration
# or a deployable single-patient probability model.

ROOT <- if (file.exists('03_scripts/00_project_config.R')) normalizePath('.') else
  'C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD'
if (!dir.exists(ROOT)) stop('Project directory not found: ', ROOT)
for (package in c('survival', 'glmnet')) {
  if (!requireNamespace(package, quietly = TRUE)) stop('Install the R package: ', package)
}
suppressPackageStartupMessages(library(survival))

input <- function(...) file.path(ROOT, ...)
meta_path <- input('02_processed_data', 'survival', 'TCGA_LUAD_PCD_survival_metadata.csv')
expr_path <- input('02_processed_data', 'survival', 'TCGA_LUAD_PCD_survival_expression.rds')
gene_path <- input('04_results', 'differential_expression',
                   'TCGA_concordant_PCD_candidates_LUAD_LUSC.csv')
required_files <- c(meta_path, expr_path, gene_path)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop('Missing source files:\n', paste(missing, collapse='\n'))

out_dir <- input('04_results', 'mpcds_calibration')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

meta <- read.csv(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
expr <- readRDS(expr_path)
gene_tab <- read.csv(gene_path, check.names=FALSE, stringsAsFactors=FALSE)
genes <- as.character(gene_tab$gene_name)
if (length(genes) != 47L || anyDuplicated(genes) || !all(genes %in% rownames(expr))) {
  stop('The original set of 47 genes does not match the expression matrix.')
}
if (!identical(as.character(colnames(expr)), as.character(meta$submitter_id))) {
  stop('Survival metadata and expression columns are not in identical patient order.')
}
if (anyDuplicated(meta$submitter_id)) stop('Duplicate patient identifiers.')

# Match the clinical eligibility and factor levels from Scripts 18/46.
stage <- as.character(meta$ajcc_pathologic_stage)
stage_group <- rep(NA_character_, length(stage))
stage_group[stage %in% c('Stage I','Stage IA','Stage IB')] <- 'Stage I'
stage_group[stage %in% c('Stage II','Stage IIA','Stage IIB')] <- 'Stage II'
stage_group[stage %in% c('Stage IIIA','Stage IIIB')] <- 'Stage III'
stage_group[stage %in% c('Stage IV')] <- 'Stage IV'

clinical_all <- data.frame(
  patient_id = as.character(meta$submitter_id),
  OS_time = as.numeric(meta$OS_time),
  OS_status = as.integer(meta$OS_status),
  age_10yr = as.numeric(meta$age_at_diagnosis_years) / 10,
  sex_at_birth = factor(as.character(meta$sex_at_birth), levels=c('female','male')),
  stage_group = factor(stage_group,levels=c('Stage I','Stage II','Stage III','Stage IV')),
  stringsAsFactors=FALSE
)
eligible <- complete.cases(clinical_all) & clinical_all$OS_time > 0 &
  clinical_all$OS_status %in% c(0L,1L)
d <- clinical_all[eligible,,drop=FALSE]
x <- t(expr[genes, d$patient_id, drop=FALSE]); mode(x) <- 'numeric'
if (!identical(rownames(x), d$patient_id) || any(!is.finite(x))) {
  stop('Invalid expression alignment or missing expression values.')
}
cat('Eligible patients:', nrow(d), '| deaths:', sum(d$OS_status), '\n')
if (nrow(d) != 475L || sum(d$OS_status) != 171L) {
  stop('Current sample differs from the published Script 46 475-patient/171-event cohort.')
}

# Identical event-stratified folds to Script 46 (repeat_id = 1).
make_stratified_folds <- function(status, k=5L, seed=12345L) {
  set.seed(seed)
  foldid <- integer(length(status))
  for (ev in sort(unique(status))) {
    idx <- which(status == ev)
    idx <- sample(idx)
    foldid[idx] <- rep(seq_len(k), length.out=length(idx))
  }
  foldid
}

# Baseline cumulative hazard H0(t) is estimated only from each TRAINING fold.
# basehaz(centered=TRUE) and predict(reference='sample') share the
# same Cox-model reference, preventing mismatched centering.
absolute_event_risk <- function(fit, newdata, t_days) {
  base <- survival::basehaz(fit, centered=TRUE)
  lp <- as.numeric(predict(fit, newdata=newdata, type='lp', reference='sample'))
  if (!length(base$time) || any(!is.finite(lp))) stop('Baseline hazard/LP unavailable.')
  vapply(t_days, function(h) {
    idx <- which(base$time <= h)
    h0 <- if (length(idx)) base$hazard[max(idx)] else 0
    risk <- 1 - exp(-h0 * exp(lp))
    pmin(pmax(risk, 0), 1)
  }, numeric(nrow(newdata)))
}

horizons <- c(365, 1095)   # days, not months
fold_id <- make_stratified_folds(d$OS_status, k=5L, seed=12345L)
pred_parts <- vector('list', 5L)

for (fold in seq_len(5L)) {
  train_idx <- which(fold_id != fold)
  test_idx <- which(fold_id == fold)
  train <- d[train_idx,,drop=FALSE]; test <- d[test_idx,,drop=FALSE]
  xtr <- x[train_idx,,drop=FALSE]; xte <- x[test_idx,,drop=FALSE]
  inner <- make_stratified_folds(train$OS_status, k=5L, seed=12345L + fold)

  ridge_cv <- glmnet::cv.glmnet(
    x=xtr, y=survival::Surv(train$OS_time,train$OS_status),
    family='cox', alpha=0, foldid=inner,
    type.measure='C', standardize=TRUE, nlambda=100,
    cox.ties='breslow')
  ridge_fit <- glmnet::glmnet(
    x=xtr, y=survival::Surv(train$OS_time,train$OS_status),
    family='cox', alpha=0, lambda=ridge_cv$lambda.min,
    standardize=TRUE, cox.ties='breslow')

  train$MPCDS_cv <- as.numeric(predict(ridge_fit, newx=xtr, type='link'))
  test$MPCDS_cv <- as.numeric(predict(ridge_fit, newx=xte, type='link'))
  mu <- mean(train$MPCDS_cv); sigma <- sd(train$MPCDS_cv)
  if (!is.finite(sigma) || sigma <= 0) stop('Invalid training-score SD.')
  train$MPCDS_cv_z <- (train$MPCDS_cv-mu)/sigma
  test$MPCDS_cv_z <- (test$MPCDS_cv-mu)/sigma

  f_clinical <- survival::coxph(
    survival::Surv(OS_time,OS_status) ~ age_10yr + sex_at_birth + stage_group,
    data=train, ties='efron', x=TRUE)
  f_score <- survival::coxph(
    survival::Surv(OS_time,OS_status) ~ MPCDS_cv_z,
    data=train, ties='efron', x=TRUE)
  f_combined <- survival::coxph(
    survival::Surv(OS_time,OS_status) ~ MPCDS_cv_z + age_10yr +
      sex_at_birth + stage_group,
    data=train, ties='efron', x=TRUE)

  p_cl <- absolute_event_risk(f_clinical,test,horizons)
  p_sc <- absolute_event_risk(f_score,test,horizons)
  p_co <- absolute_event_risk(f_combined,test,horizons)

  pred_parts[[fold]] <- data.frame(
    patient_id=test$patient_id, outer_fold=fold,
    OS_time=test$OS_time, OS_status=test$OS_status,
    clinical_lp=as.numeric(predict(f_clinical,newdata=test,type='lp')),
    combined_lp=as.numeric(predict(f_combined,newdata=test,type='lp')),
    mpcds_lp=as.numeric(predict(f_score,newdata=test,type='lp')),
    clinical_1y=p_cl[,1],clinical_3y=p_cl[,2],
    mpcds_1y=p_sc[,1],mpcds_3y=p_sc[,2],
    combined_1y=p_co[,1],combined_3y=p_co[,2],
    stringsAsFactors=FALSE)
  cat(sprintf('Completed held-out fold %d/5 (%d patients, %d deaths)\n',
              fold,nrow(test),sum(test$OS_status)))
}

pred <- do.call(rbind,pred_parts); rownames(pred) <- NULL
risk_matrix <- as.matrix(pred[,grep('_(1y|3y)$', names(pred)),drop=FALSE])
if (nrow(pred) != 475L || anyDuplicated(pred$patient_id) ||
    any(!is.finite(risk_matrix))) {
  stop('Out-of-fold predictions failed validation checks.')
}

# Optional audit against the archived Script 46 fold predictions, without
# changing its file. Difference should be machine precision if packages
# and source inputs are identical.
audit_path <- input('04_results','mpcds_incremental_uncertainty',
                    'out_of_fold_patient_predictions.csv')
if (file.exists(audit_path)) {
  old <- read.csv(audit_path, stringsAsFactors=FALSE)
  old <- old[old$repeat_id==1L,,drop=FALSE]
  m <- match(pred$patient_id, old$patient_id)
  if (anyNA(m) || any(pred$outer_fold != old$outer_fold[m])) {
    stop('Fold or patient mismatch against Script 46.')
  }
  d_cl <- max(abs(pred$clinical_lp-old$clinical_lp[m]))
  d_co <- max(abs(pred$combined_lp-old$combined_lp[m]))
  cat(sprintf('Script 46 repeat-1 LP reproduction: max differences %.3g, %.3g\n',
              d_cl,d_co))
  if (max(d_cl,d_co)>1e-6) stop('Script 46 predictions did not reproduce.')
}

write.csv(pred,file.path(out_dir,'TCGA_LUAD_calibration_out_of_fold.csv'),row.names=FALSE)

# Kaplan-Meier observed risk in the five groups ranked by held-out
# predicted event risk. Greenwood/log-log intervals are descriptive;
# grouping is not a full smooth-calibration curve.
calibration_group <- function(dat, model, horizon) {
  column <- paste0(model,'_',if (horizon==365) '1y' else '3y')
  risks <- dat[[column]]
  if (!all(is.finite(risks)) || any(risks<0 | risks>1)) stop('Invalid risks.')
  quintile <- ceiling(5*rank(risks,ties.method='first')/length(risks))
  ans <- lapply(seq_len(5L),function(q) {
    v <- dat[quintile==q,,drop=FALSE]
    fit <- survival::survfit(survival::Surv(OS_time,OS_status)~1,
                             data=v,conf.type='log-log')
    s <- summary(fit,times=horizon,extend=FALSE)
    known_beyond <- sum(v$OS_time>horizon)
    events_by <- sum(v$OS_status==1 & v$OS_time<=horizon)
    if (!length(s$surv) || known_beyond<5L || events_by<1L) {
      observed <- low <- high <- NA_real_
    } else {
      observed <- 1-s$surv[1]; low <- 1-s$upper[1]; high <- 1-s$lower[1]
    }
    data.frame(model=model,horizon_days=horizon,quintile=q,
      patients=nrow(v),events_by_horizon=events_by,
      known_followup_beyond_horizon=known_beyond,
      mean_predicted_event_risk=mean(risks[quintile==q]),
      observed_KM_event_risk=observed,observed_lower95=low,
      observed_upper95=high)
  })
  do.call(rbind,ans)
}
cal <- do.call(rbind,lapply(c('clinical','mpcds','combined'),function(model) {
  do.call(rbind,lapply(horizons,function(h) calibration_group(pred,model,h)))
}))
rownames(cal)<-NULL
write.csv(cal,file.path(out_dir,'TCGA_LUAD_calibration_groups.csv'),row.names=FALSE)

# The LP column names must match the calibration model labels.
# Prevent a missing [[column]] from silently becoming a zero-length vector.
expected_lp_cols <- paste0(c('clinical','mpcds','combined'), '_lp')
missing_lp_cols <- setdiff(expected_lp_cols, names(pred))
if (length(missing_lp_cols)) {
  stop('Missing predicted linear-predictor columns: ',
       paste(missing_lp_cols, collapse=', '))
}

# A supplementary relative-risk calibration slope for repeat 1 only.
slopes <- do.call(rbind,lapply(c('clinical','mpcds','combined'),function(model) {
  tmp <- data.frame(time=pred$OS_time,event=pred$OS_status,
                    lp=pred[[paste0(model,'_lp')]],
                    fold=factor(pred$outer_fold))
  fit <- survival::coxph(survival::Surv(time,event)~lp+strata(fold),
                         data=tmp,ties='efron')
  data.frame(model=model,slope=unname(coef(fit)['lp']),
             CI_lower=unname(confint(fit)['lp',1]),
             CI_upper=unname(confint(fit)['lp',2]))
}))
write.csv(slopes,file.path(out_dir,'TCGA_LUAD_calibration_slopes.csv'),row.names=FALSE)

# One publication-readable plot per horizon, with each model shown separately.
for (h in horizons) {
  f <- file.path(out_dir,paste0('TCGA_LUAD_calibration_',h,'days.png'))
  grDevices::png(f,width=3600,height=1300,res=300)
  oldpar <- par(no.readonly=TRUE)
  par(mfrow=c(1,3),mar=c(5,5,3,1),oma=c(0,0,2,0))
  h_cal <- cal[cal$horizon_days==h,,drop=FALSE]
  vals <- c(h_cal$mean_predicted_event_risk,h_cal$observed_upper95)
  axis_max <- min(1,max(0.40,ceiling(max(vals,na.rm=TRUE)*10)/10))
  for (model in c('clinical','mpcds','combined')) {
    dplot <- cal[cal$model==model & cal$horizon_days==h,,drop=FALSE]
    plot(c(0,axis_max),c(0,axis_max),type='n',
         xlim=c(0,axis_max),ylim=c(0,axis_max),
         asp=1,xlab='Predicted event probability',
         ylab='Observed KM event probability',
         main=switch(model,clinical='Clinical',mpcds='MPCDS',combined='Clinical + MPCDS'),
         cex.lab=1.15)
    abline(a=0,b=1,lty=2,col='gray45',lwd=1.8)
    ok <- is.finite(dplot$observed_KM_event_risk)
    if (any(ok)) {
      points(dplot$mean_predicted_event_risk[ok],dplot$observed_KM_event_risk[ok],
             pch=19,col='#006F72',cex=1.5)
      segments(dplot$mean_predicted_event_risk[ok],dplot$observed_lower95[ok],
               dplot$mean_predicted_event_risk[ok],dplot$observed_upper95[ok],
               col='#006F72',lwd=2)
    }
    if (any(!ok)) mtext('Some groups lack follow-up support',side=3,line=0,cex=.65)
  }
  mtext(paste0('TCGA-LUAD | ',h,'-day internal out-of-fold calibration'),
        side=3,outer=TRUE,line=0,cex=1.25,font=2)
  par(oldpar); grDevices::dev.off()
}
cat('\nCALIBRATION SCRIPT COMPLETED. Output folder:\n',out_dir,'\n')

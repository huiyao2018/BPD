#单因素的logit回归
logitUniVar <- function(dat, group, var, digit = 3){
  formu <- as.formula(paste0(group, " ~ ", var))
  dat[[group]] <- as.factor(dat[[group]])
  subgroup <- levels(as.factor(dat[[group]]))
  subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
  fit <- glm(formu, data = dat, family = binomial())
  unisum <- summary(fit)
  OR <- exp(coef(fit))[2]
  OR <- round(OR, digit)
  ci <- exp(confint(fit))[2,]
  ci <- round(ci, digit)
  cito <- paste0(ci[1], " - ", ci[2])
  p <- unisum$coefficients[2, "Pr(>|z|)"]
  #p <- ifelse(p < 0.001, "< 0.001", round(p, 3))
  var1 <- names(exp(coef(fit)))[2]
  result <- c(var1, group,subgroup1, OR, cito, p)
  names(result) <- c("var", "group","subgroup", "OR", "95%CI", "p.val")
  return(result)
}
# 多因素logit回归
logitMultiVar <- function(dat, group, var, adjvar,digit = 3){
  if(length(adjvar) == 1){
    formu <- as.formula(paste0(group, " ~ ", var, "+", adjvar))
  }else{
    formu <- as.formula(paste0(group, " ~ ", var, "+", paste(adjvar, collapse = "+")))
  }
  dat[[group]] <- as.factor(dat[[group]])
  subgroup <- levels(as.factor(dat[[group]]))
  subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
  fit <- glm(formu, data = dat, family = binomial())
  unisum <- summary(fit)
  OR <- exp(coef(fit))[2]
  OR <- round(OR, digit)
  ci <- exp(confint(fit))[2,]
  ci <- round(ci, digit)
  cito <- paste0(ci[1], " - ", ci[2])
  p <- unisum$coefficients[2, "Pr(>|z|)"]
  #p <- ifelse(p < 0.001, "< 0.001", round(p, 3))
  var1 <- names(exp(coef(fit)))[2]
  result <- c(var1, group,subgroup1, OR, cito, p)
  names(result) <- c("var", "group", "subgroup","OR", "95%CI", "p.val")
  return(result)
}
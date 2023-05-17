# Dependencies

library(fs) # file manipulation
library(ggplot2) # plotting 
library(dplyr) # data manipulation
library(tidyr)
library(gridExtra) # plotting
library(rlang)
library(mgcv)
library(Rcpp)
library(ipmr)
library(nlme)
library(graphics)
library(doBy)
library(boot)

get_gg_legend<-function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#sourceCpp(file = 'Ana_Israel_IPM/Cpp/cpp_utils.cpp')

# Top level wrappers for sensitivity and elasticity computations.

sensitivity <- function(K, h, level = c('kernel', 'vr', 'param'),
                        ...) {
  
  switch(level,
         "kernel" = .kernel_sens(K, h),
         "vr"     = .vr_sens(K, h, ...))
}

elasticity <- function(K, h, level = c('kernel', 'vr', 'param'),
                       ...) {
  
  switch(level,
         "kernel" = .kernel_elas(K, h),
         "vr"     = .vr_elas(K, h, ...))
}


# Sensitivity for kernel, vital rates, and parameters
.kernel_sens <- function(K, h) {
  w <- Re(eigen(K)$vectors[ , 1])
  v <- Re(eigen(t(K))$vectors[ , 1])
  
  out <- outer(v, w) / sum(v * w * h)
  
  return(out)
}


# Elasticity for kernel, vital rates, and parameters

.kernel_elas <- function(K, h) {
  
  sens <- sensitivity(K, h, 'kernel')
  
  lambda <- Re(eigen(K)$values[1])
  
  out <- sens * (K / h) / lambda
  
  return(out)
}

all_data<- read.csv("ASESORADOS/Erasmo/IPM2019/Marzo2019/Coyote/Coyote15.csv",
                    stringsAsFactors = FALSE)

# Now we are ready to explore the data a bit! first, fit some
# growth models, then plot the data to see how it looks.
hist(all_data$size, main=,xlab="Size t")
hist(all_data$sizeNext, main=,xlab="Size t+1")

grow_mod_int  <- lm(sizeNext ~ 1,
                    data = all_data[!is.na(all_data$sizeNext) &
                                      !is.na(all_data$size), ])

grow_mod_lin  <- lm(sizeNext ~ size, 
                    data = all_data)
sd_g <- sd(resid(grow_mod_lin))
sd_g <-as.list(sd_g)
names(sd_g) <- c("sd_g")

grow_mod_quad <- lm(sizeNext ~ poly(size, 2), 
                    data = all_data[!is.na(all_data$size), ])

# Now, fit a gam with a few knots to see if the linear is doing well
grow_gam <- gam(sizeNext ~ s(size, k = 8), data = all_data)

xx <- seq(-5, 135, 0.1)
plot(sizeNext ~ size, data = all_data)

abline(grow_mod_lin, col = 'red')
abline(grow_mod_int, col = 'blue')
lines(xx, predict(grow_gam, 
                  newdata = data.frame(size = xx),
                  type = 'response'),
      lty = 2,
      col = 'red')
lines(xx, predict(grow_mod_quad,
                  data.frame(size = xx),
                  type = 'response'),
      col = 'blue',
      lty = 2)

abline(a = 0, b = 1)

legend('topleft',
       legend = c('Linear Model',
                  'Intercept only',
                  'GAM', 
                  'Quadratic fit',
                  '1:1 Line'),
       col = c('red', 'blue', 'red', 'blue', 'black'),
       lty = c(1, 1, 2, 2, 1))

print(summary(grow_mod_int)) 
print(summary(grow_mod_lin))
print(summary(grow_mod_quad))
print(summary(grow_gam))

#par(mfrow = c(2,2))
plot(grow_mod_lin, ask = FALSE)

grow_aic <- AIC(grow_mod_int, grow_mod_lin, grow_mod_quad, grow_gam)

# Based on the residual plots, it looks like we have a slightly decreasing
# variance with increasing size. Next, I'll try fitting a gls with non-constant
# variance and see how that improves the fit.

grow_exp_var <- gls(sizeNext ~ size,
                    data = all_data,
                    weights = varExp(),
                    na.action = na.omit,
                    method = 'ML')


grow_aic <- AIC(grow_mod_int, grow_mod_lin,
                grow_mod_quad, grow_gam,
                grow_exp_var)

#par(mfrow = c(1, 1))

plot(sizeNext ~ size, data = all_data)

abline(grow_mod_lin, col = 'red')

lines(xx, 
      predict(grow_exp_var,
              data.frame(size = xx),
              type = 'response'),
      col = 'blue',
      lty = 2)

abline(a = 0, b = 1)

legend('topleft',
       legend = c('Linear Model - constant variance',
                  'Linear Model - non constant variance',
                  '1:1 Line'),
       col = c('red', 'blue', 'black'),
       lty = c(1, 2, 1))

# Survival models

surv_mod_int  <- glm(surv ~ 1,
                     data   = all_data,
                     family = binomial())

surv_mod_lin  <- glm(surv ~ size, 
                     data   = all_data, 
                     family = binomial())

surv_mod_quad <- glm(surv ~ size + I(size^2), 
                     data   = all_data, 
                     family = binomial())

surv_gam      <- gam(surv ~ s(size, k = 6), 
                     data   = all_data[!is.na(all_data$size), ],
                     family = binomial())

#par(mfrow = c(1,1))
plot(surv ~ size, data = all_data)

lines(xx, 
      rep(1/(1 + exp(-(coef(surv_mod_int)[1]))), length(xx)), 
      col = 'blue')

lines(xx, predict(surv_mod_lin, 
                  data.frame(size = xx),
                  type = 'response'),
      col = 'red')

lines(xx,
      predict(surv_mod_quad, 
              data.frame(size = xx), 
              type = 'response'),
      col = 'blue', 
      lty = 2)

lines(xx,
      predict(surv_gam,
              data.frame(size = xx),
              type = 'response'),
      col = 'red',
      lty = 2)

legend('bottomright',
       legend = c('Linear Model',
                  'Intercept only',
                  'GAM', 
                  'Quadratic fit'),
       col = c('red', 'blue', 'red', 'blue'),
       lty = c(1, 1, 2, 2))

# par(mfrow = c(2,2))
# plot(surv_mod, ask = FALSE)

print(summary(surv_mod_int)) 
print(summary(surv_mod_lin))
print(summary(surv_mod_quad))
print(summary(surv_gam))

surv_aic <- AIC(surv_mod_int,
                surv_mod_lin,
                surv_mod_quad,
                surv_gam)

message('Survival AIC scores\n\n')
print(surv_aic)
message('\n\nGrowth AIC scores\n\n')
print(grow_aic)

# Extract coefficients so we can generate a "data_list" for ipmr

grow_exp_var_coef_list <- coef(grow_mod_lin) %>%
  as.list() %>%
  setNames(c('g_int', 'g_slope'))


surv_coef_list <- coef(surv_mod_lin) %>%
  as.list() %>%
  setNames(c('s_int', 's_slope'))

# Fecundity parameters ----------

# Recruitment is defined as the rate of new plants(T+1) per flower at T.
#Number of seeds produced per fruit:
fec2 <- 512

#New recruits have no size(t), but do have size(t + 1)
recr_data <- subset(all_data, is.na(size))
recr_mu <- mean(recr_data$sizeNext)
recr_sd <- sd(recr_data$sizeNext)

# This data set doesn't include information on germination and establishment.
# Thus, we'll compute the realized recruitment parameter as the number
# of observed recruits divided by the number of seeds per fruits produced in the prior
# year.

recr_n <- length(recr_data$sizeNext)
flow_n <- sum(all_data$fec1, na.rm = TRUE)*fec2
recr_pr <- recr_n / flow_n

#Probability of germination within the year of seed production:
fec3 <- recr_pr

n_flow     <- sum(all_data$fec1, na.rm = TRUE)*fec2
new_plants <- subset(all_data, is.na(size))
n_new <- dim(new_plants)[1]

# Probability of reproducing. f_r gets multiplied by s_z and p_r
# to generate a fecundity kernel
p_r_mod_int <- glm(fec0 ~ 1,        data = all_data, family = binomial())
p_r_mod_lin <- glm(fec0 ~ size, data = all_data, family = binomial())

p_r_aic <- AIC(p_r_mod_int, p_r_mod_lin)

print(summary(p_r_mod_int))
print(summary(p_r_mod_lin))

message('\n\nPr(reproduction) AIC table\n\n')
print(p_r_aic)

#par(mfrow = c(1, 1))
plot(fec0 ~ size, data = all_data)
lines(xx, 
      predict(p_r_mod_lin, 
              data.frame(size = xx),
              type = 'response'),
      col = 'red')

lines(xx, 
      rep(1/(1 + exp(-(coef(p_r_mod_int)[1]))), length(xx)), 
      col = 'blue')

# New plant size distribution --------

f_d_mu <- mean(new_plants$sizeNext, na.rm = TRUE)
f_d_sd <- sd(new_plants$sizeNext, na.rm = TRUE)

# Fruit production ~ size 

f_s_int <- glm(fec1 ~ 1, data = all_data, family = poisson())
f_s_mod <- glm(fec1 ~ size, data = all_data, family = poisson())

# inspect AIC table
f_s_aic <- AIC(f_s_int, f_s_mod)

print(f_s_aic)

# compute approximate dispersion
dev_ratio <- summary(f_s_mod)$deviance / summary(f_s_mod)$df.residual

print(dev_ratio)

f_s_mod_quas <- glm(fec1 ~ size, 
                    data = all_data, 
                    family = quasipoisson())

f_s_confint  <- confint(f_s_mod_quas) 

f_s_mus      <- coef(f_s_mod_quas)
f_s_sds      <- apply(
  f_s_confint,
  1,
  function(x) {
    (x[2] - x[1]) / 3.92 # The full range is almost 4 SDs (remember, + and -!)
  }
) %>%
  as.list() %>%
  setNames(
    c(
      "f_s_int_sd",
      "f_s_slope_sd"
    )
  )

plot(fec1 ~ size, data = all_data)

lines(xx, 
      exp(coef(f_s_mod)[1] + coef(f_s_mod)[2] * xx),
      col = 'red')

print(summary(f_s_mod_quas))

fec_coef_list <- c(coef(p_r_mod_lin),
                   coef(f_s_mod_quas),
                   f_s_sds,
                   n_new,
                   f_d_mu,
                   f_d_sd) %>%
  as.list() %>%
  setNames(c('p_r_int', 'p_r_slope',
             'f_s_int', 'f_s_slope',
             "f_s_int_sd", "f_s_slope_sd",
             'n_new', 'f_d_mu', 'f_d_sd'))

all_param_list  <- c(fec_coef_list, grow_exp_var_coef_list, sd_g, surv_coef_list)

# End parameter estimation!

#IPM
s_z <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

L <- min(c(all_data$size,
           all_data$sizeNext),
         na.rm = TRUE) * 1.2

U <- max(c(all_data$size,
           all_data$sizeNext),
         na.rm = TRUE) * 1.2

n_mesh_p <- 100

carp_ipmr <- init_ipm("simple", "di", "det") %>%
  define_kernel(
    name = "P",
    formula = S * G,
    family = "CC",
    G = dnorm(sa_2, mu_g, sd_g),
    mu_g = g_int + g_slope * sa_1,
    S = inv_logit(s_int, s_slope, sa_1),
    data_list = all_param_list,
    states = list(c('sa')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "G")
  ) %>%
  define_kernel(
    name = "F",
    formula = f_r * f_s * f_d * p_r,
    family = "CC",
    f_s = exp(f_s_int + f_s_slope * sa_1),
    f_d = dnorm(sa_2, f_d_mu, f_d_sd),
    p_r = inv_logit(p_r_int, p_r_slope, sa_1),
    f_r = n_new / sum(f_s),
    data_list = all_param_list,
    states = list(c("sa")),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "f_d")
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule = rep('midpoint', 2),
      state_start = rep('sa', 2),
      state_end = rep('sa', 2) 
    ) 
  ) %>%
  define_domains(
    sa = c(L, U, n_mesh_p)
  )  %>%
  define_pop_state(
    n_sa = rep(1/100, n_mesh_p)
  ) %>%
  make_ipm(usr_funs = list(inv_logit = s_z),
           iterate = TRUE,
           iterations = 1000)

lambda_ipmr <- lambda_exp_var <- lambda(carp_ipmr)
lambda_ipmr

# Option 1: is_conv_to_asymptotic
is_conv_to_asymptotic(carp_ipmr)

# Option 2: generate iteration kernel and compute eigenvalues
K <- make_iter_kernel(carp_ipmr)
lam_eigen <- Re(eigen(K$mega_matrix)$values[1])

# If we've iterated our model enough, this should be approximately 0 (though
# maybe a little off due to floating point errors).
lambda_ipmr - lam_eigen

# Sub-kernels have their own print method to display the range of values
# and some diagnotic information.
carp_ipmr$sub_kernels

#Further analysis
R_nought <- function(ipm_obj) {
  Pm <- ipm_obj$sub_kernels$P
  Fm <- ipm_obj$sub_kernels$F
  I <- diag(dim(Pm)[1])
  N <- solve(I - Pm)
  R <- Fm %*% N
  return(
    Re(eigen(R)$values)[1]
  )
}

gen_time <- function(ipm_obj) {
  lamb <- unname(lambda(ipm_obj))
  r_nought <- R_nought(ipm_obj)
  return(log(r_nought) / log(lamb))
}

R0 <- R_nought(carp_ipmr)
gen_T <- gen_time(carp_ipmr)
R0
#R0 can be defined as the
#average number of offspring that an individual produces over their lifetime.2
#That is, each individual in the present generation contributes R0 individuals,
#on average, to the next generation.

gen_T
#Another measure of generation
#time is the mean age of mothers at offspring production.

make_N <- function(ipm) {
  
  P     <- ipm$sub_kernels$P
  
  I     <- diag(nrow(P))
  N     <- solve(I - P)
  
  return(N)
}

eta_bar <- function(ipm) {
  
  N     <- make_N(ipm)
  out   <- colSums(N)
  
  return(as.vector(out))
  
}

sigma_eta <- function(ipm) {
  
  N     <- make_N(ipm)  
  
  out <- colSums(2 * (N %^% 2L) - N) - colSums(N) ^ 2
  
  return(as.vector(out))
}

mean_l <- eta_bar(carp_ipmr)
var_l  <- sigma_eta(carp_ipmr)

mesh_ps <- int_mesh(carp_ipmr)$sa_1 %>%
  unique()

#par(mfrow = c(1,2))

plot(mesh_ps, mean_l, type = "l", xlab = expression( "Initial size z"[0]))
plot(mesh_ps, var_l, type = "l", xlab = expression( "Initial size z"[0]))

#browseVignettes("ipmr")

v_a <- left_ev(carp_ipmr, iterations = 1000)
w_a <- right_ev(carp_ipmr, iterations = 1000)

#par(mfrow = c(1, 2))

plot(1:100, seq(0, max(unlist(w_a)), length.out = 100), type = "n",
     ylab = expression(paste("w"[a],"(z)")),
     xlab = "Size bin")

for(i in seq_along(w_a)) {
  
  lines(w_a[[i]])
  
}

plot(1:100, 
     seq(0, max(unlist(v_a)), length.out = 100),
     type = "n",
     ylab = expression(paste("v"[a], "(z)")),
     xlab = "Size bin")

for(i in seq_along(v_a)) {
  
  lines(v_a[[i]]) 
  
}


#PLOTS
lab_seq <- round(seq(L, U, length.out = 6), 2)
tick_seq <- c(1, 20, 40, 60, 80, 100)
#par(mfrow = c(2, 2))
# Sub-kernels - ipmr contains plot methods for sub-kernels
plot(carp_ipmr$sub_kernels$P,
     do_contour = TRUE,
     main = "P",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

plot(carp_ipmr$sub_kernels$F,
     do_contour = TRUE,
     main = "F",
     xlab = "size (t)",
     ylab = "size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

# Iterate over meshpoints to make sure matrix size doesn't matter

mesh_p <- seq(100, 500, by = 50)

lambdas <- data.frame(mesh_p = mesh_p,
                      lambdas = c(lambda_ipmr, rep(NA_real_, 8)))

i <- 1
it<- 1

for(i in mesh_p) {
  
  
  carp_ipmr_test <- init_ipm("simple", "di", "det") %>%
    define_kernel(
      name = "P",
      formula = S * G,
      family = "CC",
      G = dnorm(sa_2, mu_g, sd_g),
      mu_g = g_int + g_slope * sa_1,
      S = inv_logit(s_int, s_slope, sa_1),
      data_list = all_param_list,
      states = list(c('sa')),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>%
    define_kernel(
      name = "F",
      formula = f_r * f_s * f_d * p_r,
      family = "CC",
      f_s = exp(f_s_int + f_s_slope * sa_1),
      f_d = dnorm(sa_2, f_d_mu, f_d_sd),
      p_r = inv_logit(p_r_int, p_r_slope, sa_1),
      f_r = n_new / sum(f_s),
      data_list = all_param_list,
      states = list(c("sa")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "f_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('sa', 2),
        state_end = rep('sa', 2) 
      ) 
    ) %>%
    define_domains(
      sa = c(L, U, i) # iterate over different meshpoint numbers
    ) %>%
    define_pop_state(
      n_sa = rep(1/i, i)
    ) %>%
    make_ipm(usr_funs = list(inv_logit = s_z))
  
  lambdas$lambdas[it] <- lambda(carp_ipmr_test)
  
  it <- it + 1
  
}

plot(lambdas ~ mesh_p, data = lambdas)

diff(lambdas$lambdas)

# increasing from 100 to 150 does not really make a difference for lambda, so
# i'm going to keep it at 100 meshpoints

#Kernel-level perturbations
# Sensitivity and elasticity

# First, do kernel sensitivities. Function definitions are 
# in 01_Utils_and_Dependencies

K_exp_var <- do.call(`+`, carp_ipmr$sub_kernels)
P_exp_var <- carp_ipmr$sub_kernels$P
F_exp_var <- carp_ipmr$sub_kernels$F

class(K$mega_matrix) <- c("ipmr_matrix", class(K$mega_matrix))
plot(K_exp_var,
     do_contour = TRUE,
     main = "IPM matrix",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))


bs <- seq(L, U, length.out = n_mesh_p + 1)
mesh <- d1 <- (bs[2:(n_mesh_p + 1)] + bs[1:n_mesh_p]) * 0.5

h <- bs[2] - bs[1]

#Local prospective perturbation analysis (hereafter "perturbation analysis")
#amounts to calculating the partial derivative of a quantity with respect to a
#specific perturbation.
#In the ecological literature, the term sensitivity is used for a partial derivative
#(e.g., the partial derivative of ?? with respect to mean soil moisture). In many
#cases a more meaningful quantity is the elasticity or proportional sensitivity,
#the fractional change in the response relative to the fractional change in the
#quantity being perturbed
#The most common perturbation analyses ask how the asymptotic population
#growth rate, ??, of a density-independent model changes when we perturb the
#kernel, K(z, z). 
k_sens_exp_mat   <- sensitivity(K_exp_var, h, level = 'kernel') 

# Elasticities

k_elas_exp_mat   <- elasticity(K_exp_var, h, level = 'kernel') 


# Sub-kernel elasticities

P_fun_ev <- P_exp_var / h
F_fun    <- F_exp_var / h

# Individual kernel elasticities
P_elas_ev_mat <- P_fun_ev * k_sens_exp_mat / lambda_exp_var

F_elas_ev_mat <- F_fun * k_sens_exp_mat / lambda_exp_var

#To better quantify these differences we can integrate over each elasticity surface separately using sum(P.elas)*h^2 and sum(F.elas)*h^2. These show that the
#relative contributions are 0.971 and 0.029, respectively, for survival-growth and
#reproduction. These sum to one (as they must) and demonstrate the relatively
#small effect that proportional changes to F would make to ??.
sum(P_elas_ev_mat)*h^2 
sum(F_elas_ev_mat)*h^2

# Sensitivity and elasticity
class(k_sens_exp_mat) <- c("ipmr_matrix", class(k_sens_exp_mat))
class(k_elas_exp_mat) <- c("ipmr_matrix", class(k_elas_exp_mat))
plot(k_sens_exp_mat,
     do_contour = TRUE,
     main = "K Sensitivity",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

plot(k_elas_exp_mat,
     do_contour = TRUE,
     main = "K Elasticity",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

#Figure shows the P and F elasticity functions. In contrast to the
#sensitivity surfaces, the elasticity surfaces approach 0 at rare transitions. This
#makes sense because even a very large proportional perturbation of the kernel
#where it is near 0 still generates a very small absolute perturbation. Panel C
#clearly shows the survival-growth elasticity surface dominates the total elasticity
#e(z0 , z0) - it is indistinguishable from panel B. The elasticity surface associated
#with F again shows that transitions from the most common size classes of ewes
#have the largest effect on ??, but that overall it makes only a small contribution
#to e(z0 , z0).

class(P_elas_ev_mat) <- c("ipmr_matrix", class(P_elas_ev_mat))
plot(P_elas_ev_mat,
     do_contour = TRUE,
     main = "P Elasticity",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

class(F_elas_ev_mat) <- c("ipmr_matrix", class(F_elas_ev_mat))
plot(F_elas_ev_mat,
     do_contour = TRUE,
     main = "F Elasticity",
     xlab = "Size (t)",
     ylab = "Size (t + 1)",
     yaxt = "none",
     xaxt = "none")
axis(1, at = tick_seq, labels = as.character(lab_seq))
axis(2, at = tick_seq, labels = as.character(lab_seq))

#Graphs with ggplot2
p_df <- ipm_to_df(carp_ipmr$sub_kernels$P)
f_df <- ipm_to_df(carp_ipmr$sub_kernels$F)
k_df <- ipm_to_df(K$mega_matrix)
sens_df <- ipm_to_df(k_sens_exp_mat)
elas_df <- ipm_to_df(k_elas_exp_mat)

# Create a default theme for our plots
def_theme <- theme(
  panel.background = element_blank(),
  axis.text = element_text(size = 16),
  axis.ticks = element_line(size = 1.5),
  axis.ticks.length = unit(0.08, "in"),
  axis.title.x = element_text(
    size = 20,
    margin = margin(
      t = 10,
      r = 0,
      l = 0,
      b = 2
    )
  ),
  axis.title.y = element_text(
    size = 20,
    margin = margin(
      t = 0,
      r = 10,
      l = 2,
      b = 0
    )
  ),
  legend.text = element_text(size = 16)
)

p_plt <- ggplot(p_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("P kernel")

f_plt <- ggplot(f_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("F kernel")

k_plt <- ggplot(k_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K kernel")

sens_plt <- ggplot(sens_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K Sensitivity")

elas_plt <- ggplot(elas_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 0.7,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "red",
                      high = "yellow") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) +
  ggtitle("K Elasticity")

p_plt
f_plt
k_plt

#Figure shows the calculated kernel sensitivity function. The shape of
#the sensitivity surface is governed by w(z) and v(z). A slice through this surface along a line of constant z yields a curve that is proportional to the stable
#size distribution. A slice along a line of constant z produces a curve that is
#proportional to the relative reproductive value, which is an increasing function
#of body size.
#The resulting surface shows that absolute changes to the kernel
#which alter transitions into large size classes from the most abundant size class
#have the greatest relative impact on population growth rate.
sens_plt

#Figure shows the resulting kernel elasticity function. This suggests
#that proportional changes to the survival-growth component of the kernel will
#have the greatest impact on ??, with transitions from the most common size class
#being most important.
elas_plt


#vital_rate_exprs(carp_ipmr)

# Bootstrapping code

# Split the data into plants that existed at T and new recruits at T+1

exists_t_1  <- filter(all_data, !is.na(size))
recruit_t_2 <- filter(all_data, is.na(size))

# Number of existing plants to resample

n_resamp_existing <- dim(exists_t_1)[1]
n_resamp_recruits <- dim(recruit_t_2)[1]

n_resamp <- 1000


# Set up place to hold boot strapping outputs

temp_exp_var_out <- c(grow_exp_var_coef_list, sd_g,
                      surv_coef_list,
                      fec_coef_list,
                      lambda = lambda_exp_var,
                      p_elas = sum(P_elas_ev_mat) * h ^ 2,
                      f_elas = sum(F_elas_ev_mat) * h ^ 2)

exp_var_out <- lapply(temp_exp_var_out,
                      function(x, n_resamp) {
                        c(x, rep(NA_real_, n_resamp))
                      },
                      n_resamp = n_resamp)

for(i in seq_len(n_resamp)) {
  
  # Resample our data with replacement and create a single data frame
  
  ex_resamp_ind <- sample(1:n_resamp_existing,
                          n_resamp_existing, 
                          replace = TRUE)
  re_resamp_ind <- sample(1:n_resamp_recruits,
                          n_resamp_recruits, 
                          replace = TRUE)
  
  boot_existing <- exists_t_1[ex_resamp_ind, ]
  boot_recruits <- recruit_t_2[re_resamp_ind, ]
  
  boot_data <- rbind(boot_existing, boot_recruits)
  
  # Bootstrap growth model
  
  boot_grow_mod_exp <- lm(sizeNext ~ size, data = boot_data)
  # Bootstrap survival model
  
  boot_surv_mod <- glm(surv ~ size, 
                       data   = all_data, 
                       family = binomial())
  
  # Bootstrap fecundity parameters
  
  n_flow <- sum(boot_data$fec1, na.rm = TRUE)*fec2
  n_new<- length(boot_recruits)
  
  # Not computing a number of new recruits because we've constrained that in our
  # procedure. Thus, only the number of flowers can vary from sample to sample
  
  boot_f_r <- n_new / n_flow
  
  boot_p_r_mod <- glm(fec0 ~ size, data = boot_data, family = binomial())
  
  boot_f_d_mu <- mean(boot_recruits$sizeNext, na.rm = TRUE)
  boot_f_d_sd <- sd(boot_recruits$sizeNext, na.rm = TRUE)
  
  boot_f_s_mod <- glm(fec1 ~ size, 
                      data = boot_data, 
                      family = poisson())
  
  # We already have the implementation details from when we generated the point
  # estimate of lambda, so no need to redefine those.
  
  all_param_list <- list(
    s_int    = coef(boot_surv_mod)[1],
    s_slope  = coef(boot_surv_mod)[2],
    g_int    = coef(boot_grow_mod_exp)[1],
    g_slope  = coef(boot_grow_mod_exp)[2],
    sd_g = sd(resid(boot_grow_mod_exp)),
    n_new =length(boot_recruits),
    L = L,
    U = U,
    p_r_int = coef(boot_p_r_mod)[1],
    p_r_slope = coef(boot_p_r_mod)[2],
    f_s_int = coef(boot_f_s_mod)[1],
    f_s_slope = coef(boot_f_s_mod)[2],
    f_r = boot_f_r,
    f_d_mu = boot_f_d_mu,
    f_d_sd = boot_f_d_sd
  )
  
  
  carp_boot <- init_ipm("simple", "di", "det") %>%
    define_kernel(
      name = "P",
      formula = S * G,
      family = "CC",
      G = dnorm(sa_2, mu_g, sd_g),
      mu_g = g_int + g_slope * sa_1,
      S = inv_logit(s_int, s_slope, sa_1),
      data_list = all_param_list,
      states = list(c('sa')),
      evict_cor =  TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>%
    define_kernel(
      name = "F",
      formula = f_r * f_s * f_d * p_r,
      family = "CC",
      f_s = exp(f_s_int + f_s_slope * sa_1),
      f_d = dnorm(sa_2, f_d_mu, f_d_sd),
      p_r = inv_logit(p_r_int, p_r_slope, sa_1),
      f_r = n_new / sum(f_s),
      data_list = all_param_list,
      states = list(c("sa")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "f_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('sa', 2),
        state_end = rep('sa', 2) 
      ) 
    ) %>%
    define_domains(
      sa = c(L, U, n_mesh_p)
    ) %>% 
    define_pop_state(n_sa = rep(1/100, 100)) %>% 
    make_ipm(usr_funs = list(inv_logit = s_z))
  
  K_exp_var_boot <- do.call(`+`, carp_boot$sub_kernels)
  P_exp_var_boot <- carp_boot$sub_kernels$P
  F_exp_var_boot <- carp_boot$sub_kernels$F
  
  # Store lambdas
  exp_var_out$lambda[(i + 1)]   <- l_ev_boot <- lambda(carp_boot)
  
  k_sens_ev_boot <- sensitivity(K_exp_var_boot, h, level = 'kernel')
  
  P_fun_ev_boot <- P_exp_var_boot / h
  F_fun_boot    <- F_exp_var_boot / h
  
  p_elas_ev_mat_boot <- P_fun_ev_boot * k_sens_ev_boot / l_ev_boot
  f_elas_ev_mat_boot <- F_fun_boot * k_sens_ev_boot / l_ev_boot
  
  exp_var_out$p_elas[(i + 1)] <- sum(p_elas_ev_mat_boot) * h ^ 2
  exp_var_out$f_elas[(i + 1)] <- sum(f_elas_ev_mat_boot) * h ^ 2
  
  # Store exponential variance growth parameters
  exp_var_out$g_int[(i + 1)]       <- coef(boot_grow_mod_exp)[1]
  exp_var_out$g_slope[(i + 1)]     <- coef(boot_grow_mod_exp)[2]  
  exp_var_out$sd_g[(i + 1)] <- sd(resid(boot_grow_mod_exp))
  
  # Survival models are the same across both 
  exp_var_out$s_int[(i + 1)]   <- coef(boot_surv_mod)[1]
  exp_var_out$s_slope[(i + 1)] <- coef(boot_surv_mod)[2] 
  
  # Same with pr(repro) and the other fecundity parameters
  exp_var_out$p_r_int[(i + 1)]   <- coef(boot_p_r_mod)[1]
  exp_var_out$p_r_slope[(i + 1)] <- coef(boot_p_r_mod)[2] 
  
  exp_var_out$f_s_int[(i + 1)]   <- coef(boot_f_s_mod)[1]
  exp_var_out$f_s_slope[(i + 1)] <- coef(boot_f_s_mod)[2]
  
  exp_var_out$f_r[(i + 1)]    <- boot_f_r
  exp_var_out$f_d_mu[(i + 1)] <- boot_f_d_mu
  exp_var_out$f_d_sd[(i + 1)] <- boot_f_d_sd
  exp_var_out$n_new[(i + 1)] <- n_new
  
}


exp_var_temp <- lapply(exp_var_out,
                       function(x) {
                         temp <- sort(x[2:1001])
                         return(c(x[1], temp[25], temp[975]))
                       }) %>%
  as_tibble() %>%
  mutate(boot_obs = c('Observed', 'Lower_CI', 'Upper_CI'))

exp_var_temp$lambda
lambda_ipmr

# Final plots

surv_pred <- data.frame(size = xx, 
                        pred = predict(surv_mod_lin, 
                                       data.frame(size = xx),
                                       type = 'response'))

grow_pred <- data.frame(size = xx,
                        pred = predict(grow_exp_var,
                                       data.frame(size = xx),
                                       type = 'response'))

pr_pred <- data.frame(size = xx,
                      pred = predict(p_r_mod_lin,
                                     data.frame(size = xx),
                                     type = 'response'))
fs_pred <- data.frame(size = xx,
                      pred = predict(f_s_mod,
                                     data.frame(size = xx),
                                     type = 'response'))

yy <- seq(min(recruit_t_2$sizeNext, na.rm = TRUE) - 10,
          max(recruit_t_2$sizeNext, na.rm = TRUE) + 10,
          length.out = 400)

recr_pred <- data.frame(sizeNext = yy,
                        density = dnorm(yy, 
                                        f_d_mu,
                                        f_d_sd))


all_vr <- exp_var_temp %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)

lam_pred <- filter(all_vr, vital_rate == 'lambda')

elas_plot <- filter(all_vr, vital_rate %in% c('p_elas', 'f_elas'))


# Set up default themes for contour plots and line plots

theme_contour <- theme(
  panel.background = element_blank(),
  axis.text        = element_text(size   = 8),
  axis.title.x     = element_text(size   = 10,
                                  margin = margin(
                                    t = 5,
                                    r = 0, 
                                    l = 0, 
                                    b = 1
                                  )
  ),
  axis.title.y     = element_text(size   = 8,
                                  margin = margin(
                                    t = 0,
                                    r = 5,
                                    l = 1,
                                    b = 0
                                  )
  )
)

theme_linerange <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_blank(), # Remove x-axis text + title
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 7), # make y-axis text + title bigger
    axis.title.y      = element_text(size = 9,
                                     margin = margin(
                                       t = 0,
                                       l = 2,
                                       r = 5,
                                       b = 0
                                     )),
    strip.text        = element_text(size = 10), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7)
  )

theme_vr <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_text(size = 7), 
    axis.text.y       = element_text(size = 7), # make y-axis text + title bigger
    axis.title.x      = element_text(size   = 8,
                                     margin = margin(
                                       t = 5,
                                       r = 0, 
                                       l = 0, 
                                       b = 7
                                     )
    ),
    axis.title.y     = element_text(size   = 8,
                                    margin = margin(
                                      t = 3,
                                      r = 10,
                                      l = 1,
                                      b = 0
                                    )
    ),
    strip.text        = element_text(size = 8), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7)
  )

# Now, make the figure panel

grow_plot <- ggplot(all_data, aes(x = size,
                                  y = sizeNext)) +
  geom_point(color = 'black', 
             size = 1.25) + 
  geom_line(data = grow_pred,
            aes(x = size,
                y = pred),
            size = 1.25,
            color = 'grey70',
            show.legend = FALSE) + 
  #geom_abline(intercept = 0,
  #            slope = 1,
  #            color = 'grey70',
  #            show.legend = FALSE,
  #            size = 1.25) +
  theme_vr +
  scale_x_continuous('Size t',
                     limits = c(0, 100)) +
  scale_y_continuous('Size t + 1',
                     limits = c(0, 100))

grow_plot

surv_plot <- ggplot(all_data, aes(x = size,
                                  y = surv)) +
  geom_point(color = 'black',
             size = 1.75,
             width = 0,
             height = 0.05) + 
  geom_line(data = surv_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'gray70') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 100)) +
  scale_y_continuous('Survival (t + 1)',
                     limits = c(0, 1),
                     breaks = c(0, 1))

surv_plot

pr_plot <- ggplot(all_data, aes(x = size, 
                                y = fec0)) +
  geom_point(color = 'black',
             size = 1.75,
             width = 0,
             height = 0.05) + 
  geom_line(data = pr_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'gray70') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 100)) +
  scale_y_continuous('Pr(Reproductive, t)',
                     limits = c(0, 1.1),
                     breaks = c(0, 1))

pr_plot

fs_plot <- ggplot(all_data, aes(x = size, 
                                y = fec1)) +
  geom_point(color = 'black',
             size = 1.75) + 
  geom_line(data = fs_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'gray70') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 100)) +
  scale_y_continuous('Number of fruits',
                     limits = c(0, 40),
                     breaks = seq(0, 70, by = 10))

fs_plot
#END


#Graphs with ggplot2
p_df <- ipm_to_df(carp_ipmr$sub_kernels$P)
f_df <- ipm_to_df(carp_ipmr$sub_kernels$F)
k_df <- ipm_to_df(K$mega_matrix)
sens_df <- ipm_to_df(k_sens_exp_mat)
elas_df <- ipm_to_df(k_elas_exp_mat)

# Create a default theme for our plots
def_theme <- theme(
  panel.background = element_blank(),
  axis.text = element_text(size = 26),
  axis.ticks = element_line(size = 1.5),
  axis.ticks.length = unit(0.08, "in"),
  axis.title.x = element_text(
    size = 30,
    margin = margin(
      t = 10,
      r = 0,
      l = 0,
      b = 2
    )
  ),
  axis.title.y = element_text(
    size = 30,
    margin = margin(
      t = 0,
      r = 10,
      l = 2,
      b = 0
    )
  ),
  legend.text = element_text(size = 20)
)

p_plt <- ggplot(p_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 1,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "blue",
                      high = "magenta") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank())

f_plt <- ggplot(f_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 1,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "blue",
                      high = "magenta") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank()) 

k_plt <- ggplot(k_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 1,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "blue",
                      high = "magenta") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank())

sens_plt <- ggplot(sens_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 1,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "blue",
                      high = "magenta") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank())

elas_plt <- ggplot(elas_df) +
  geom_tile(aes(x = t,
                y = t_1,
                fill = value)) +
  geom_contour(aes(x = t,
                   y = t_1,
                   z = value),
               color = "black",
               size = 1,
               bins = 5) +
  scale_fill_gradient("Value",
                      low = "blue",
                      high = "magenta") +
  scale_x_continuous(name = "Size (t)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  scale_y_continuous(name = "Size (t + 1)",
                     labels = lab_seq,
                     breaks = tick_seq) +
  def_theme +
  theme(legend.title = element_blank())

p_plt
f_plt
k_plt

#Figure shows the calculated kernel sensitivity function. The shape of
#the sensitivity surface is governed by w(z) and v(z). A slice through this surface along a line of constant z yields a curve that is proportional to the stable
#size distribution. A slice along a line of constant z produces a curve that is
#proportional to the relative reproductive value, which is an increasing function
#of body size.
#The resulting surface shows that absolute changes to the kernel
#which alter transitions into large size classes from the most abundant size class
#have the greatest relative impact on population growth rate.
sens_plt

#Figure shows the resulting kernel elasticity function. This suggests
#that proportional changes to the survival-growth component of the kernel will
#have the greatest impact on ??, with transitions from the most common size class
#being most important.
elas_plt

###LAST PLOTS
surv_pred <- data.frame(size = xx, 
                        pred = predict(surv_mod_lin, 
                                       data.frame(size = xx),
                                       type = 'response'))

grow_pred <- data.frame(size = xx,
                        pred = predict(grow_exp_var,
                                       data.frame(size = xx),
                                       type = 'response'))

pr_pred <- data.frame(size = xx,
                      pred = predict(p_r_mod_lin,
                                     data.frame(size = xx),
                                     type = 'response'))
fs_pred <- data.frame(size = xx,
                      pred = predict(f_s_mod,
                                     data.frame(size = xx),
                                     type = 'response'))

yy <- seq(min(recruit_t_2$sizeNext, na.rm = TRUE) - 10,
          max(recruit_t_2$sizeNext, na.rm = TRUE) + 10,
          length.out = 400)

recr_pred <- data.frame(sizeNext = yy,
                        density = dnorm(yy, 
                                        f_d_mu,
                                        f_d_sd))


all_vr <- exp_var_temp %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)

lam_pred <- filter(all_vr, vital_rate == 'lambda')

elas_plot <- filter(all_vr, vital_rate %in% c('p_elas', 'f_elas'))

# Set up default themes for contour plots and line plots

theme_contour <- theme(
  panel.background = element_blank(),
  axis.text        = element_text(size   = 26),
  axis.title.x     = element_text(size   = 30,
                                  margin = margin(
                                    t = 5,
                                    r = 0, 
                                    l = 0, 
                                    b = 1
                                  )
  ),
  axis.title.y     = element_text(size   = 30,
                                  margin = margin(
                                    t = 0,
                                    r = 5,
                                    l = 1,
                                    b = 0
                                  )
  )
)

theme_linerange <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_blank(), # Remove x-axis text + title
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 26), # make y-axis text + title bigger
    axis.title.y      = element_text(size = 30,
                                     margin = margin(
                                       t = 0,
                                       l = 2,
                                       r = 5,
                                       b = 0
                                     )),
    strip.text        = element_text(size = 25), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 26),
    legend.title      = element_text(size = 27)
  )

theme_vr <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_text(size = 26), 
    axis.text.y       = element_text(size = 26), # make y-axis text + title bigger
    axis.title.x      = element_text(size   = 30,
                                     margin = margin(
                                       t = 5,
                                       r = 0, 
                                       l = 0, 
                                       b = 7
                                     )
    ),
    axis.title.y     = element_text(size   = 30,
                                    margin = margin(
                                      t = 3,
                                      r = 10,
                                      l = 1,
                                      b = 0
                                    )
    ),
    strip.text        = element_text(size = 28), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 26),
    legend.title      = element_text(size = 27)
  )

# Now, make the figure panel

grow_plot <- ggplot(all_data, aes(x = size,
                                  y = sizeNext)) +
  geom_point(color = 'black', 
             size = 4) + 
  geom_line(data = grow_pred,
            aes(x = size,
                y = pred),
            size = 1.25,
            color = 'blue',
            show.legend = FALSE) + 
  #geom_abline(intercept = 0,
  #            slope = 1,
  #            color = 'grey70',
  #            show.legend = FALSE,
  #            size = 1.25) +
  theme_vr +
  scale_x_continuous('Size t',
                     limits = c(0, 123)) +
  scale_y_continuous('Size t + 1',
                     limits = c(0, 123))

grow_plot

surv_plot <- ggplot(all_data, aes(x = size,
                                  y = surv)) +
  geom_point(color = 'black',
             size = 3,
             width = 0,
             height = 0.05) + 
  geom_line(data = surv_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'blue') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 125)) +
  scale_y_continuous('Survival (t + 1)',
                     limits = c(0, 1),
                     breaks = c(0, 1))

surv_plot

pr_plot <- ggplot(all_data, aes(x = size, 
                                y = fec0)) +
  geom_point(color = 'black',
             size = 3,
             width = 0,
             height = 0.05) + 
  geom_line(data = pr_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'blue') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 125)) +
  scale_y_continuous('Pr(Flowering)',
                     limits = c(0, 1.1),
                     breaks = c(0, 1))

pr_plot

fs_plot <- ggplot(all_data, aes(x = size, 
                                y = fec1)) +
  geom_point(color = 'black',
             size = 3) + 
  geom_line(data = fs_pred,
            aes(x = size,
                y = pred),
            #linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'blue') +
  theme_vr + 
  scale_x_continuous('Size t', 
                     limits = c(0, 125)) +
  scale_y_continuous('Fruits number',
                     limits = c(0, 55),
                     breaks = seq(0, 70, by = 10))

fs_plot

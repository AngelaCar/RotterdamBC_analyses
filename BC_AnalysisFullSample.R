# Preliminaries -----------------------------------------------------------

# install.packages("TwoTimeScales")
library("TwoTimeScales")
library("fields")       # for plotting functions
library("survival")     # for dataset
library("tictoc")       # for tic() toc() functions to measure estimation time

# Breast cancer dataset -- preliminary operations -------------------------
d <- rotterdam

# Address concern about the 43 individuals with censoring time for recurrence
# smaller than death times
# Censor at the time of recurrence 
d[d$rtime < d$dtime & d$recur == 0, ]$death <- 0
d[d$rtime < d$dtime & d$recur == 0, ]$dtime <-  d[d$rtime < d$dtime & d$recur == 0, ]$rtime

# Time variables
d$rtimey <- d$rtime/365.25        # time at recurrence in years
d$dtimey <- d$dtime/365.25        # time from surgery to death in years

# Age variables
d$rage <- d$age + d$rtimey        # age at recurrence
d$dage <- d$age + d$dtimey        # age at death

d$grade <- as.factor(d$grade)     

# For competing risks analysis
d$fetimey <- pmin(d$rtimey, d$dtimey)  # time from surgery to first event in years

# Type of first event
d$first_event <- ifelse(d$recur == 1, "recurrence", NA)
d$first_event <- ifelse(d$recur == 0 & d$death == 1, "death", d$first_event)
d[is.na(d$first_event),]$first_event <- "censored"

# Example 1: All deaths (with or without recurrence) ----------------------
# ----- Prepare the data for estimation -----------------------------------
death_cov <- prepare_data(data = d,
                          u = "age",
                          s_out = "dtimey",
                          events = "death",
                          ds = .5,
                          du = 1,
                          min_u = 24,
                          min_s = 0,
                          individual = TRUE,
                          covs = c("grade"))
print(death_cov)

# ----- Estimation --------------------------------------------------------
tic()
mod_deathcov <- fit2ts(data2ts = death_cov,
                       Bbases_spec = list(bdeg = 3,
                                          nseg_u = 12,
                                          min_u = 24,
                                          nseg_s = 7,
                                          min_s = 0),
                       pord = 2,
                       lrho = c(1,-2),
                       optim_method = "ucminf",
                       optim_criterion = "bic"
)
toc()
summary(mod_deathcov)
#save(mod_deathcov, file = "BC_AllDeathsMod_BIC.Rda")

# ----- Visualization ----------------------------------------------------

# Reproduces Figure 3a
# Baseline hazard in (u,s)-plane
#pdf("Haz-us_BIC.pdf", width = 6.5, height = 6) #- comment off to save pdf file
par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2)
plot(mod_deathcov,
     plot_grid = list(c(umin = 24, umax = 90, du = .2),
                      c(smin = 0, smax = 19.5, du = .1)),
     plot_options = list(n_shades = 100,
                         tmax = 90,
                         cut_extrapolated = T,
                         main = "Baseline hazard of death",
                         xlab = "Age at surgery (years)",
                         ylab = "Time since surgery (years)",
                         xlim = c(24, 90)))

#dev.off()

# Reproduced Figure 3b
# Associated SEs to the baseline hazard
#pdf("SE_Haz-us_BIC.pdf", width = 6.5, height = 6) #- comment off to save pdf file
par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2)
plot(mod_deathcov,
     which_plot = "SE",
     plot_grid = list(c(umin = 24, umax = 90, du = .2),
                      c(smin = 0, smax = 19.5, du = .1)),
     plot_options = list(n_shades = 100,
                         tmax = 90,
                         cut_extrapolated = T,
                         main = "Standard Errors of the baseline hazard of death",
                         xlab = "Age at surgery (years)",
                         ylab = "Time since surgery (years)",
                         xlim = c(24, 90)))

#dev.off()

# ----- Predictions -------------------------------------------------------

newdata <- as.data.frame(expand.grid("age" = c(40, 50, 60), 
                                     "dtimey" = seq(0, 15, by = .1),
                                     "grade_3" = 1))
newdata$t <- newdata$age + newdata$dtimey
newdata <- subset(newdata, t <= 70)
predicted_haz <- predict(mod_deathcov_BIC,
                          newdata = newdata,
                          u = "age", s = "dtimey",
                          z = "grade_3")
round(head(predicted_haz),3)

# Reproduces Figure 4 left
#pdf("PredictedHazardEx1.pdf", width = 6, height = 5)
par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2,
    mar = c(8, 4, 4, 1),
    oma = c(0, 0, 0, 0)+.2)

plot(1,1, type = "n",
     xlim = c(0, 15),
     ylim = c(0, max(predicted_haz$hazard)),
     main = "Predicted hazard \n Grade = 3",
     xlab = "Time since surgery",
     ylab = "Hazard"
     )
with(subset(predicted_haz, age == 40), 
     lines(dtimey, hazard,
           lwd = 2,
           lty = 1, 
           col = viridis(3)[3]))
with(subset(predicted_haz, age == 50), 
     lines(dtimey, hazard,
           lwd = 2,
           lty = 2, 
           col = viridis(3)[2]))
with(subset(predicted_haz, age == 60), 
     lines(dtimey, hazard,
           lwd = 2,
           lty = 3, 
           col = viridis(3)[1]))
legend("bottom", horiz = TRUE,
       cex = 1.2,
  xpd = TRUE,
       inset = c(0,-0.6),
       legend = c(40,50,60),
       title = "Age at surgery",
       lwd = 2,
       lty = 1:3,
       bty = "n",
       col = rev(viridis(3)))
#dev.off()

# Prediction for age at surgery = 50 with confidence intervals
predicted_haz$lci <- predicted_haz$hazard - 1.96 * predicted_haz$se_hazard
predicted_haz$uci <- predicted_haz$hazard + 1.96 * predicted_haz$se_hazard

# Reproduces Figure 5 right
#pdf("PredictedHazardEx1_Age50withCI.pdf",  width = 6, height = 5)

par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2,
    mar = c(8, 4, 4, 1),
    oma = c(0, 0, 0, 0)+.2)

plot(1,1, type = "n",
     xlim = c(0, 15),
     ylim = c(0, max(predicted_haz$uci)),
     main = "Predicted hazard \n Grade = 3",
     xlab = "Time since surgery",
     ylab = "Hazard"
)
with(subset(predicted_haz, age == 50), 
     lines(dtimey, hazard,
           lwd = 2,
           lty = 1, 
           col = viridis(3)[2]))
with(subset(predicted_haz, age == 50), 
     lines(dtimey, lci,
           lwd = 1,
           lty = 1, 
           col = viridis(3)[2]))
with(subset(predicted_haz, age == 50), 
     lines(dtimey, uci,
           lwd = 1,
           lty = 1, 
           col = viridis(3)[2]))
with(subset(predicted_haz, age == 50), 
     polygon(c(rev(dtimey), dtimey), c(rev(lci), uci),
        col = adjustcolor(viridis(3)[2], alpha = .5), 
        border = NA
)
)
legend("bottom", horiz = TRUE,
       cex = 1.2,
       xpd = TRUE,
       inset = c(0,-0.6),
       legend = c(50),
       title = "Age at surgery",
       lwd = 2,
       lty = 1,
       bty = "n",
       col = viridis(3)[2])

#dev.off()

# Alternative models (with AIC and LMMsolver) ----------------------------
tic()
mod_deathcov_LMM <- fit2ts(data2ts = death_cov,
                           Bbases_spec = list(bdeg = 3,
                                              nseg_u = 12,
                                              min_u = 24,
                                              nseg_s = 7,
                                              min_s = 0),
                           pord = 2,
                           optim_method = "LMMsolver",
                           control_algorithm = list(monitor_ev = TRUE)
)
toc()
summary(mod_deathcov_LMM)

tic()
mod_deathcov_AIC <- fit2ts(data2ts = death_cov,
                           Bbases_spec = list(bdeg = 3,
                                              nseg_u = 12,
                                              min_u = 24,
                                              nseg_s = 7,
                                              min_s = 0),
                           pord = 2,
                           lrho = c(1,-2),
                           optim_criterion = "aic",
                           control_algorithm = list(monitor_ev = TRUE)
)
toc()
summary(mod_deathcov_AIC)

# Reproduces Figure 6 (Appendix), left and right panel
#pdf("Haz-us_AIC.pdf", width = 6.5, height = 6) #- comment off to save pdf file
par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2)
plot(mod_deathcov,
     plot_grid = list(c(umin = 24, umax = 90, du = .2),
                      c(smin = 0, smax = 19.5, du = .1)),
     plot_options = list(n_shades = 100,
                         tmax = 90,
                         cut_extrapolated = T,
                         main = "Hazard of death\n optim: AIC",
                         xlab = "Age at surgery (years)",
                         ylab = "Time since surgery (years)",
                         xlim = c(24, 90)))

#dev.off()

#pdf("Haz-us_LMM.pdf", width = 6.5, height = 6) #- comment off to save pdf file
par(font.main = 1,
    cex.main = 1.8, cex.lab = 1.2)
plot(mod_deathcov_LMM,
     plot_grid = list(c(umin = 24, umax = 90, du = .2),
                      c(smin = 0, smax = 19.5, du = .1)),
     plot_options = list(n_shades = 100,
                         tmax = 90,
                         cut_extrapolated = T,
                         main = "Hazard of death\n optim: LMMsolver",
                         xlab = "Age at surgery (years)",
                         ylab = "Time since surgery (years)",
                         xlim = c(24, 90)))

#dev.off()

# Example 2: Competing risks after surgery --------------------------------
# ---- Data preparation ---------------------------------------------------
# The data have to be prepared separately for the two different event types
rec2ts <- prepare_data(data = d,
                       u = "age",
                       s_out = "fetimey",
                       event = "recur",
                       min_s = 0,
                       min_u = 24,
                       du = 1, ds = .5)
print(rec2ts)

death2ts <- prepare_data(data = d,
                         u = "age",
                         s_out = "fetimey",
                         event = "death",
                         min_s = 0,
                         min_u = 24,
                         du = 1, ds = .5)
print(death2ts)

# ---- Estimation ---------------------------------------------------------
# Model for recurrence
mod_rec <- fit2ts(rec2ts,
                  Bbases_spec = list(bdeg = 3,
                                     nseg_u = 12,
                                     min_u = 24,
                                     nseg_s = 7,
                                     min_s = 0),
                  optim_criterion = "bic")

summary(mod_rec)

# Model for death
mod_death <- fit2ts(death2ts,
                    Bbases_spec = list(bdeg = 3,
                                       nseg_u = 12,
                                       min_u = 24,
                                       nseg_s = 7,
                                       min_s = 0),
                    optim_criterion = "bic")

summary(mod_death)

# ---- Competing risks quantities -----------------------------------------
# Calculate the competing risks quantities 
# Cumulative hazards - separately by event type
H_rec <- cumhaz2ts(mod_rec,
                   plot_grid = list(c(umin = 24, umax = 90, du = .2),
                                    c(smin = 0,  smax = 19.5, ds = .1)),
                   cause = "recurrence")
H_death <- cumhaz2ts(mod_death,
                     plot_grid = list(c(umin = 24, umax = 90, du = .2),
                                      c(smin = 0,  smax = 19.5, ds = .1)),
                     cause = "death")
# Cumulative incidence functions - both obtained from one call to cuminc2ts
cif <- cuminc2ts(haz = list(H_rec$Haz$hazard, H_death$Haz$hazard),
                 ds = .1,
                 cause = c("recurrence", "death"))

# ---- Plot of the Cumulative Incidence Functions -------------------------
# Reproduces Figure 5 
# Here we opt for a different plotting strategy, we use the image() plotting function
# together with contour() to plot the two cif surfaces on the same plotting window

# CIF surface of recurrence
pl1 <- cif[[1]]
# CIF surface of death (without recurrence)
pl2 <- cif[[2]]

max(max(pl1), max(pl2))                   # identify common range for cif levels
brk <- seq(0, 0.75, length = 101)         # color legend breaks
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples")[1:7]) # color palette
mypurp <- "#2C0158"                       # color for contour lines

x <- H_rec$Haz$new_plot_grid$intu         
y <- H_rec$Haz$new_plot_grid$ints

#pdf("BC-CIF.pdf", width = 10, height = 4.5)
par(font.main = 1,
    oma = c(0,0,0,1.1),
    cex.lab = 1.2, cex.main = 1.5)

split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
split.screen(c(1,2), screen=1)-> ind
# first image
screen(ind[1])
image(x, y, pl1,
      col = pal(100),
      breaks = brk,
      main = "Recurrence",
      xlab = "Age at surgery (years)",
      ylab = "Time since surgery (years)",
      xlim = c(24, 90))
contour(x, y, pl1,
        col = mypurp,
        add = TRUE)
# second image
screen(ind[2])
image(x, y, pl2,
      col = pal(100),
      breaks = brk,
      main = "Death before recurrence",
      xlab = "Age at surgery (years)",
      ylab = "Time since surgery (years)",
      xlim = c(24, 90))
contour(x, y, pl2,
        col = mypurp,
        add = TRUE)
# draw the legend strip
screen(2)
par(pin = c(0, 3.91))
image.plot(zlim = c(0, 0.75), legend.only = TRUE, smallplot = c(.0, .15, .1, .9),
           col = pal(100))
#dev.off()

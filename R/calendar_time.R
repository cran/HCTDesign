#' Calculate Calendar Times for Interim Analysis Looks
#'
#' @param k Numeric vector of event fractions for each look (e.g., c(0.5, 1))
#' @param d2max maximum number of events in the experimental group calculated from the design function.
#' @param htime Historical survival times
#' @param hevent Historical event indicators (0/1)
#' @param delta Hazard ratio (experimental vs. historical)
#' @param ta Accrual time
#' @param tf Follow-up time
#' @param method "exponential", "log_normal" or "KM" for survival function estimation(Default is KM)
#' @param event_ind Event indicator value (default = 1)
#' @author Tushar Patni, Yimei Li and Jianrong Wu.
#'
#' @return
#'   Data frame with columns: Distribution, look, time fraction, events, calendar time
#' @importFrom survival survfit Surv
#' @importFrom flexsurv flexsurvreg
#' @import stats
#' @references
#' \insertRef{doi:10.1002/pst.1756}{HCTDesign}
#' @references
#' \insertRef{doi:10.1080/10543406.2019.1684305}{HCTDesign}
#' @importFrom Rdpack reprompt
#' @importFrom diversitree set.defaults
#' @examples
#' time <- c(20, 65, 12, 50, 58, 65, 45, 44)
#' event <- c(1, 0, 0, 0, 1, 1, 1, 1)
#' gg<-calendar_time(k=c(0.5, 1),d2max=46,htime=time,hevent=event,delta=0.7,ta=5,tf=4,method="KM")
#' @export
calendar_time <- function(k, d2max, htime, hevent, delta, ta, tf, method = c("KM"), event_ind = 1) {
  method <- match.arg(method, c("KM", "exponential", "log_normal"))

  # browser()

  # Calculate n if not provided

  n <- as.numeric(SM(
    time = htime, event = hevent, d2max = d2max, opt = method,
    event_ind = event_ind, ta = ta, tf = tf, delta = delta
  ))

  # Calculate accrual rate
  r <- n / ta
  if (method == "exponential") {
    # Exponential fit
    total_time_h <- sum(htime)
    lambda_h <- sum(hevent) / total_time_h
    lambda_e <- delta * lambda_h
    # Expected events function for exponential
    expected_events_exp <- function(t, r_val = r, tacc = ta, lambda = lambda_e) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        int_0_t <- (1 - exp(-lambda * t)) / lambda
        exposure <- t - int_0_t
      } else {
        int_ta_t <- (exp(-lambda * t) / lambda) * (exp(lambda * tacc) - 1)
        exposure <- tacc - int_ta_t
      }
      return(r_val * exposure)
    }
  } else if (method == "KM") {
    # KM fit
    SurvObj <- Surv(htime, hevent)
    model <- survfit(SurvObj ~ 1, conf.type = "log-log")
    p0 <- c(1, summary(model)$surv)
    t0 <- c(0, summary(model)$time)
    outKM <- data.frame(t0 = t0, p0 = p0)
    KM <- function(t) {
      k_len <- length(outKM$t0)
      if (t < 0) {
        return(0)
      }
      if (t >= outKM$t0[k_len]) {
        return(outKM$p0[k_len])
      }
      for (i in 1:(k_len - 1)) {
        if (t >= outKM$t0[i] & t < outKM$t0[i + 1]) {
          return(outKM$p0[i])
        }
      }
      return(1)
    }
    S1 <- function(t) KM(t)^delta
    int_S_e <- function(a, b) {
      mid <- (a + b) / 2
      result <- ((b - a) / 6) * (S1(a) + (4 * S1(mid)) + S1(b))
      return(result)
    }
    expected_events_km <- function(t, r_val = r, tacc = ta) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        return(r_val * (t - int_S_e(0, t)))
      } else {
        return(r_val * (tacc - int_S_e(t - tacc, t)))
      }
    }
  } else {
    s1 <- Surv(htime, hevent == event_ind)
    model <- flexsurvreg(s1 ~ 1, dist = "lognormal") # fit the log-normal model #
    summary(model)
    meanlog <- model$res[1, 1]
    sdlog <- model$res[2, 1]
    lognormal <- function(t) {
      1 - pnorm((log(t) - meanlog) / sdlog)
    }
    S2 <- function(t) {
      lognormal(t)^delta
    }
    int_S_e <- function(a, b) {
      integrate(S2, lower = a, upper = b)$value
    }
    expected_events_logn <- function(t, r_val = r, tacc = ta) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        return(r_val * (t - int_S_e(0, t)))
      } else {
        return(r_val * (tacc - int_S_e(t - tacc, t)))
      }
    }
  }
  expected_events<-switch(method,
    "exponential" = expected_events_exp,
    "KM" = expected_events_km,
    "log_normal" = expected_events_logn
  )

 # browser()
  # Calculate expected events function for delta = 1 (H0)
  if (method == "exponential") {
    lambda_e_h0 <- 1 * lambda_h # delta = 1
    expected_events_h0_exp <- function(t, r_val = r, tacc = ta, lambda = lambda_e_h0) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        int_0_t <- (1 - exp(-lambda * t)) / lambda
        exposure <- t - int_0_t
      } else {
        int_ta_t <- (exp(-lambda * t) / lambda) * (exp(lambda * tacc) - 1)
        exposure <- tacc - int_ta_t
      }
      return(r_val * exposure)
    }
  } else if (method == "KM") {
    # KM fit for H0
    S1_h0 <- function(t) KM(t)^1 # delta = 1
    int_S_e_h0 <- function(a, b) {
      mid <- (a + b) / 2
      result <- ((b - a) / 6) * (S1_h0(a) + (4 * S1_h0(mid)) + S1_h0(b))
      return(result)
    }
    expected_events_h0_km <- function(t, r_val = r, tacc = ta) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        return(r_val * (t - int_S_e_h0(0, t)))
      } else {
        return(r_val * (tacc - int_S_e_h0(t - tacc, t)))
      }
    }
  } else {
    s1 <- Surv(htime, hevent == event_ind)
    model <- flexsurvreg(s1 ~ 1, dist = "lognormal") # fit the log-normal model #
    summary(model)
    meanlog <- model$res[1, 1]
    sdlog <- model$res[2, 1]
    lognormal <- function(t) {
      1 - pnorm((log(t) - meanlog) / sdlog)
    }
    S2_h0 <- function(t) {
      lognormal(t)^1
    }
    int_S_e_h0 <- function(a, b) {
      integrate(S2_h0, lower = a, upper = b)$value
    }
    expected_events_h0_logn <- function(t, r_val = r, tacc = ta) {
      if (t <= 0) {
        return(0)
      }
      if (t <= tacc) {
        return(r_val * (t - int_S_e_h0(0, t)))
      } else {
        return(r_val * (tacc - int_S_e_h0(t - tacc, t)))
      }
    }
  }
  expected_events_h0<-switch(method,
    "exponential" = expected_events_h0_exp,
    "KM" = expected_events_h0_km,
    "log_normal" = expected_events_h0_logn
  )

  # Calculate calendar times for each look
  results <- data.frame(
    "Distribution" = method,
    "Look" = seq_along(k),
    `Time fraction` = k,
    "Events" = k * d2max, check.names = FALSE
  )
  max_time <- ta + tf + (ta + tf) # Reasonable upper bound
  for (i in seq_along(k)) {
    d_look <- k[i] * d2max
    # Calendar time under Ha (original delta)
    tryCatch(
      {
        sol <- uniroot(
          function(t) expected_events(t) - d_look,
          interval = c(0, max_time),
          tol = 1e-6
        )
        results$`Calendar time under H1`[i] <- sol$root
      },
      error = function(e) {
        warning(sprintf("Could not solve for look %d (Ha): %s", i, e$message))
      }
    )
    # Calendar time under H0 (delta = 1)
    tryCatch(
      {
        sol_h0 <- uniroot(
          function(t) expected_events_h0(t) - d_look,
          interval = c(0, max_time),
          tol = 1e-6
        )
        results$`Calendar time under H0`[i] <- sol_h0$root
      },
      error = function(e) {
        warning(sprintf("Could not solve for look %d (H0): %s", i, e$message))
      }
    )
  }
  return(results)
}

#calendar_time<-set.defaults(calendar_time,method="KM")
##############################

#' Compute rejection target function
#'
#' @param theta parameters in weights and G
#' @param weights vector of the original weights
#' @param mat_g transition matrix G
#' @param p_values vector of p-values for the elementary hypotheses
#' @param alpha overall alpha level
#' @param utility rejection utility of each elementary hypothesis
#' @param null_h null hypothesis index
#'
#' @return vector of the following values: 1. rejection utility, by default,
#'     this is the average number of rejected hypothesis 2. reject any null
#'     hypothesis rate, i.e., familywise false discovery rate 3. average
#'     rejection rate for each hypothesis
#'
#' @export
#'
om_rejection <- function(weights, mat_g, p_values, alpha = 0.05,
                         utility = 1, null_h = NULL) {

    rej <- apply(p_values, 1, function(x) {
        c_mtp(x,
              mat_g  = mat_g,
              alphas = alpha * weights)
    })
    rej <- t(rej)

    ## rejection rate separately
    rej_rates <- apply(rej, 2, mean)

    ## rejection any null
    if (0 == length(null_h)) {
        rej_any <- 0
    } else {
        rej_any <- apply(rej[, null_h, drop = FALSE],
                         1, sum)
        rej_any <- mean(rej_any > 0)
    }

    ## rejection utility
    inx_h <- seq_len(length(weights))
    if (0 < length(null_h))
        inx_h <- inx_h[-null_h]

    if (0 == length(inx_h)) {
        rej_uti <- 0
    } else {
        utility <- rep(utility, length.out = length(inx_h))
        rej_uti <- sum(utility * rej_rates[inx_h])
    }

    ## return
    c(rej_uti     = rej_uti,
      rej_anynull = rej_any,
      rej_rates)
}


#' Fill parameters to the G and weights
#'
#' @inheritParams om_rejection
#'
#' @export
#'
om_theta_fill <- function(theta = NULL, weights, mat_g, tot = 1) {
    fill_vec <- function(vec, ta) {
        if (is.null(ta))
            return(vec)

        inx_na <- which(is.na(vec))
        n_na   <- length(inx_na)

        if (1 == n_na) {
            vec[inx_na] <- tot - sum(vec[-inx_na])
        } else if (1 < n_na) {
            i_ta            <- seq_len(n_na - 1)
            vec[inx_na[-1]] <- ta[i_ta]
            vec[inx_na[1]]  <- tot - sum(vec[-inx_na[1]])
            ta              <- ta[-i_ta]
        }

        list(vec, ta)
    }

    n_h <- length(weights)
    stopifnot(n_h == ncol(mat_g))

    ## fill in weights
    fw       <- fill_vec(weights, theta)
    weights  <- fw[[1]]
    theta    <- fw[[2]]

    ## fill G
    rst_g <- NULL
    for (i in seq_len(nrow(mat_g))) {
        fw    <- fill_vec(mat_g[i, ], theta)
        rst_g <- rbind(rst_g, fw[[1]])
        theta <- fw[[2]]
    }

    ## return
    list(weights = weights,
         mat_g   = rst_g)
}

#' Fill parameter then calculate rejection
#'
#' @inheritParams om_rejection
#'
#' @export
#'
om_theta_rejection <- function(theta, weights, mat_g, ...) {
    fill_theta <- om_theta_fill(theta, weights, mat_g)
    om_rejection(fill_theta$weights, fill_theta$mat_g, ...)
}


#' Get optim constraints
#'
#' @inheritParams om_rejection
#'
#' @export
#'
om_theta_constraints <- function(mat, tot = 1) {
    n_theta <- NULL
    row_sum <- NULL
    for (i in seq_len(nrow(mat))) {
        vec     <- mat[i, ]
        inx_na  <- which(is.na(vec))
        n_na    <- length(inx_na)
        if (n_na < 2)
            next

        n_theta <- c(n_theta, n_na - 1)
        row_sum <- c(row_sum, sum(vec, na.rm = TRUE))
    }

    if (is.null(n_theta)) {
        return(list(tot_theta = 0))
    }

    tot_theta <- sum(n_theta)
    tot_rows  <- length(n_theta)

    ## all parameters >= 0
    rst_ci <- rep(0, tot_theta)
    rst_ui <- diag(tot_theta)

    ## all sum < 1
    rst_ci <- c(rst_ci, -tot + row_sum)

    for (i in seq_len(tot_rows)) {
        if (1 == i) {
            inx_1 <- 1
        } else {
            inx_1 <- 1 + sum(n_theta[1:(i-1)])
        }

        inx_2 <- inx_1 + n_theta[i] - 1

        c_ui              <- rep(0, tot_theta)
        c_ui[inx_1:inx_2] <- -1

        ## append
        rst_ui <- rbind(rst_ui, c_ui)
    }

    ## return
    return(list(
        ui = rst_ui,
        ci = rst_ci,
        tot_theta = tot_theta
    ))
}

#' Optimize weights and G
#'
#' @inheritParams om_rejection
#'
#' @param par_optim options for constrOptim function
#' @param tot restrict on the sum of parameters
#'
#' @export
#'
om_rejection_optim <- function(weights, mat_g, ...,
                               tot = 1,
                               par_optim = list()) {

    f_opt <- function(theta) {
        - om_theta_rejection(theta, weights, mat_g, ...)[1]
    }

    ## get constraints
    ui_ci <- om_theta_constraints(rbind(weights, mat_g, tot = tot))

    if (0 == ui_ci$tot_theta) {
        theta <- NULL
        value <- - f_opt(NULL)
    } else {
        init_theta <- rep(0.000001, rep = ui_ci$tot_theta)
        if (ui_ci$tot_theta > 1) {
            rst <- do.call(constrOptim, c(
                                            list(
                                                f     = f_opt,
                                                grad  = NULL,
                                                theta = init_theta,
                                                ui = ui_ci$ui,
                                                ci = ui_ci$ci
                                            ),
                                            par_optim
                                        ))
        } else {
            rst <- do.call(optim, c(
                                      list(
                                          par    = init_theta,
                                          fn     = f_opt,
                                          gr     = NULL,
                                          method = "Brent",
                                          lower  = 0,
                                          upper  = max(abs(ui_ci$ci))
                                      ),
                                      par_optim
                                  ))

        }

        if (0 != rst$convergence) {
            theta <- NA
            value <- NA
        } else {
            theta <- rst$par
            value <- rst$value
        }
    }

    ## return
    list(theta   = theta,
         rej_uti = value)
}

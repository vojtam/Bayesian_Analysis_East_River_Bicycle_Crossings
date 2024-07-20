box::use(
  tibble[tibble],
  dplyr[...],
  tidyr[...],
  ggplot2[ggplot, geom_line, aes, geom_contour, labs, ggtitle, geom_point],
  stats[dgamma, qgamma, dpois],
  HDInterval[hdi],
)


# POISSON -----------------------------------------------------------------


#' get_poisson_prior_tab
#'
#' vytvori dataframe s apriornim rozdelenim pravdepodobnosti
#'
#' @param theta_start levy konec intervalu parametrickeho prostoru
#' @param theta_end pravy konec intervalu parametrickeho prostoru
#' @param a apriorni parametr gamma rozdeleni == celkovy pocet udalosti ve
#'  fiktivnim apriornim vyber
#' @param b apriorni parametr gamma rozdeleni == rozsah fiktivniho apriorniho
#'  vyberu (\strong{degree of belief})
#' @param theta_res "by" parametr pro \emph{seq} rikajici, jak huste chceme
#' nasamplovat parametr \theta. Default: 0.01
#' @param model_name jmeno modelu, ktere se vlozi do sloupce \emph{model}.
#' Default: "prior"
#' @return dataframe se sloupci: n | y | p | theta | model
#' @export
#'
#' @examples
#' get_poisson_prior_tab(0, 5, 2, 1)
get_poisson_prior_tab <-
  function(theta_start,
           theta_end,
           a,
           b,
           theta_res = 0.1,
           model_name = "prior",
           is_jeffreys = F
           ) {
    theta <- seq(theta_start, theta_end, by = theta_res)

    if (is_jeffreys) {
      d.prior <- data.frame(theta = theta) |> 
        mutate(
          n = NA, 
          y = NA, 
          p = 1 / sqrt(theta), 
          model = model_name
        )
      return(d.prior)
    }
    d.prior <- data.frame(theta = theta) |> 
      mutate(
        n = NA, 
        y = NA, 
        p = dgamma(theta, a, b), 
        model = model_name
      )
    return(d.prior)
  }



#'
#' Funkce vytvori a vrati DataFrame s posteriornim rozdelenim pravdepodobnosti
#' pro poissonuv model
#'
#' @param a apriorni parametr gamma rozdeleni == celkovy pocet udalosti ve
#'  fiktivnim apriornim vyber
#' @param b apriorni parametr gamma rozdeleni == rozsah fiktivniho apriorniho
#'  vyberu (\strong{degree of belief})
#' @param n = rozsah nahodneho vyberu v experimentu
#' @param y = celkovy pocet udalosti pozorovanych v experimentu
#' @param theta vektor nebo list nasamplovanych hodnot parametru theta
#' @param model_name jmeno modelu, ktere se vlozi do sloupce \emph{model}.
#' Default: "posterior"
#'
#' @return dataframe se sloupci: n | y | p | theta | model
#' @export
#'
#' @examples
get_poisson_posterior_tab <-
  function(a, b, n, y, theta, model_name = "posterior") {
    posterior_tab <- tibble(n = n, y = y, theta = theta) |>
      mutate(p = dgamma(theta, a + y, b + n),
             model = model_name)
  }



#' get_poisson_likelihood_tab
#'
#' Vypocita verohodnostni funkci (likelihood) podle:
#' (Y \mid \theta) \sim Po(n \theta)
#' f(y \mid \theta) = \frac{(n\theta)^y}{y!} e^{-n\theta}
#'
#'
#' @param n = rozsah nahodneho vyberu v experimentu
#' @param y = celkovy pocet udalosti pozorovanych v experimentu
#' @param theta vektor nebo list nasamplovanych hodnot parametru theta
#' @param likelihood_multiplier konstanta, kterou vynasobime likelihood pro
#' lepsi vizualizaci. Vynasobeni konstantou se na poloze verohodnostni funkce
#' neprojevi
#' @param model_name jmeno modelu, ktere se vlozi do sloupce \emph{model}.
#' Default: "likelihood"
#'
#' @return dataframe se sloupci: n | y | p | theta | model
#' @export
#'
#' @examples
get_poisson_likelihood_tab <-
  function(n,
           y,
           theta,
           likelihood_multiplier = 1,
           model_name = "likelihood") {
    likelihood_tab <- tibble(
      theta = theta,
      model = model_name,
      n = n,
      y = y
    ) |>
      mutate(p = dpois(y, n * theta) * likelihood_multiplier)
  }




#' get_distribution plot
#'
#' @param full_tab dataframe, ktery obsahuje sloupce n | y | p | theta | model
#' a kde model = c("prior", "posterior", "likelihood")
#' @param linewidth ggplot line witdth
#'
#' @return ggplot2 plot posteriorni distribuce, apriorni distribuce a likelihood
#' @export
#'
#' @examples
get_distribution_plot <- function(full_tab, linewidth = 1) {
  plot <-
    full_tab |> ggplot(aes(
      x = theta,
      y = p,
      color = model,
      group = model
    )) +
    geom_line(linewidth = linewidth)
}


#' get_prior_poisson_characteristics
#'
#' @param prior_a apriorni parametr gamma rozdeleni == celkovy pocet udalosti ve
#'  fiktivnim apriornim vyber
#' @param prior_b apriorni parametr gamma rozdeleni == rozsah fiktivniho apriorniho
#'  vyberu (\strong{degree of belief})
#'
#' @return dataframe se sloupci: mean | SD | mode | lwr | upr | HDP.lower | HDP.upper
#' @export
#'
#' @examples
get_prior_poisson_characteristics <- function(prior_a, prior_b, model_name = "prior") {
  prior <- list(a = prior_a, b = prior_b)
  
  stat.prior <- with(
    prior,
    data.frame(
      model = model_name,
      MLE = NA,
      MAP = NA,
      mean = a / b,
      SD = sqrt(a / b ^ 2),
      ETCI.lwr = qgamma(0.025, a, b),
      ETCI.upr = qgamma(0.975, a, b),
      HPD = t(hdi(
        qgamma,
        credMass = 0.95,
        shape = a,
        rate = b
      ))
    )
  )
}



#' get_posterior_poisson_characteristics
#'
#' @param prior_a apriorni parametr gamma rozdeleni == celkovy pocet udalosti ve
#'  fiktivnim apriornim vyber
#' @param prior_b apriorni parametr gamma rozdeleni == rozsah fiktivniho apriorniho
#'  vyberu (\strong{degree of belief})
#' @param n rozsah nahodneho vyberu v experimentu
#' @param y celkovy pocet udalosti pozorovanych v experimentu
#' @param full_tab dataframe, ktery obsahuje sloupce n | y | p | theta | model
#' a kde model = c("prior", "posterior", "likelihood")
#'
#' @return Dataframe se sloupci: mean | SD | mode | lwr | upr | MLE | HDP.lower | HDP.upper
#' @export
#'
#' @examples
get_posterior_poisson_characteristics <-
  function(prior_a,
           prior_b,
           n,
           y,
           full_tab,
           model_name = "posterior") {
    posterior_hyperparameters <- list(a = prior_a + y, b = prior_b + n)
    MaxL <- slice(.data = full_tab, which.max(full_tab$p))
    
    stats_posterior <- with(
      posterior_hyperparameters,
      data.frame(
        model = model_name,
        MLE = MaxL$theta,
        MAP = (a - 1) / b,
        mean = a / b,
        SD = sqrt(a / (b ^ 2)),
        ETCI.lwr = qgamma(0.025, a, b),
        ETCI.upr = qgamma(0.975, a, b),
        HPD = t(hdi(\(x) qgamma(x, a, b), credMass = 0.95))
      )
    )
  }



#' calculate_posterior_mean
#'
#' @param mean.prior
#' @param belief == b
#' @param y
#' @param n
#'
#' @return
#' @export
#'
#' @examples
calculate_posterior_mean <- function(mean.prior, belief, y, n) {
  (belief * mean.prior + y) / (belief + n)
}


#' predict_neg_bin_prior
#'
#' @param pred_start
#' @param pred_end
#' @param prior_a
#' @param prior_b
#' @param model_name
#'
#' @return
#' @export
#'
#' @examples
predict_neg_bin_prior <-
  function(pred_start = 0,
           pred_end = 10,
           prior_a,
           prior_b,
           model_name = "marginal") {
    tibble(
      y.pred = seq(pred_start, pred_end, by = 1),
      f = dnbinom(y.pred, size = prior_a, prob = prior_b / (prior_b + 1)),
      model = model_name
    )
  }


#' predict_neg_bin_marginal_characteristics
#'
#' @param prior_a
#' @param prior_b
#' @param model_name
#'
#' @return
#' @export
#'
#' @examples
predict_neg_bin_marginal_characteristics <-
  function(prior_a, prior_b, model_name = "marginal") {
    marginal <- list(a = prior_a, b = prior_b)
    with(
      marginal,
      data.frame(
        model = model_name,
        mean = a / b,
        SD = sqrt(a * (b + 1) / b ^ 2),
        lower = qnbinom(0.025, a, b / (b + 1)),
        upper = qnbinom(0.975, a, b / (b + 1))
      )
    )
  }


#' predict_neg_bin_posterior_characteristics
#'
#' @param prior_a
#' @param prior_b
#' @param y
#' @param n
#' @param model_name
#'
#' @return
#' @export
#'
#' @examples
predict_neg_bin_posterior_characteristics <-
  function(prior_a, prior_b, y, n, model_name = "posterior") {
    posterior1 <- list(a = prior_a + y, b = prior_b + n)
    with(
      posterior1,
      data.frame(
        model = model_name,
        mean = a / b,
        SD = sqrt(a * (b + 1) / b ^ 2),
        lower = qnbinom(0.025, a, b / (b + 1)),
        upper = qnbinom(0.975, a, b / (b + 1))
      )
    )
  }


#' predict_neg_bin_posterior
#'
#' @param pred_start
#' @param pred_end
#' @param prior_a
#' @param prior_b
#' @param y
#' @param n
#' @param model_name
#'
#' @return
#' @export
#'
#' @examples
predict_neg_bin_posterior <-
  function(pred_start = 0,
           pred_end = 10,
           prior_a,
           prior_b,
           y,
           n,
           how_many = 1,
           model_name = "posterior") {
    tibble(
      y.pred = seq(pred_start, pred_end, by = 1),
      f = dnbinom(
        y.pred,
        size = prior_a + y,
        prob = (prior_b + n) / (prior_b + n + how_many)
      ),
      model = model_name
    )
  }


#' get_predict_dist_plots
#'
#' @param predict_tab tabulka se sloupci: y.pred | f | model
#' @param lab_x popisek osy x
#' @param lab_y popisek osy y
#'
#' @return ggplot2 plot prediktivnich rozdelni
#' @export
#'
#' @examples
get_predict_dist_plots <-
  function(predict_tab, lab_x, lab_y = "prediktivni pravdepodobnosti") {
    predict_tab |> ggplot(aes(
      x = y.pred,
      y = f,
      group = model,
      colour = model
    )) +
      geom_line(linewidth = 0.3, lty = 2) +
      geom_point(size = 1.5, pch = 19) +
      labs(x = lab_x, y = lab_y)
  }


#' get_poisson_sensitivity_plot
#'
#' @param a apriorni parametr gamma rozdeleni == celkovy pocet udalosti ve
#'  fiktivnim apriornim vyber
#' @param b apriorni parametr gamma rozdeleni == rozsah fiktivniho apriorniho
#'  vyberu (\strong{degree of belief})
#' @param n = rozsah nahodneho vyberu v experimentu
#' @param y = celkovy pocet udalosti pozorovanych v experimentu
#' @param mean_prior_left zacatek intervalu pro samplovani apriorni stredni hodnoty
#' @param mean_prior_right konec intervalu pro samplovani apriorni stredni hodnoty
#' @param belief_left zacatek intervalu pro samplovani degree of belief (b)
#' @param belief_right
#'
#' @return
#' @export
#'
#' @examples
get_poisson_sensitivity_plot <-
  function(a,
           b,
           n,
           y,
           mean_prior_left,
           mean_prior_right,
           belief_left,
           belief_right) {
    da <- expand_grid(
      mean.prior = seq(mean_prior_left, mean_prior_right, by = 1),
      belief = seq(belief_left, belief_right, by = 1),
    ) |> mutate(mean.posterior = calculate_posterior_mean(mean.prior, belief, y, n))
    # vrstevnicovy graf s barevnymi vrstevnicemi
    da |> ggplot(aes(x = belief, y = mean.prior, z = mean.posterior)) +
      geom_contour(binwidth = 0.05, aes(colour = after_stat(level))) +
      labs(x = "belief", y = "apriorni stredni hodnota", color = "apost.\nstredni\nhodnota") +
      ggtitle("Aposteriorni stredni hodnota") +
      geom_point(
        x = b,
        y = a / b,
        color = "red",
        show.legend = FALSE
      )
  }



# Priklady ----------------------------------------------------------------


# 
# 
# 
# 
# # 1. zeny
# 
# load(file = "../../04/deti1.RData", verbose = TRUE)
# data <- Y
# 
# tab <- data |>
#   filter(FEMALE == 1) |>
#   filter((YEAR >= 1990 & YEAR < 2000) & (AGE == 40)) |>
#   select(CHILDS, DEGREE)
# 
# tab <- drop_na(tab)
# 
# Y1 <- tab |> filter(DEGREE <= 2) |> mutate(edu = "<BC")
# n1 <- nrow(Y1)
# y1 <- sum(Y1$CHILDS)
# 
# Y2 <- tab |> filter(DEGREE >= 3) |> mutate(edu = "BC+")
# n2 <- nrow(Y2)
# y2 <- sum(Y2$CHILDS)
# 
# Y <- bind_rows(Y1, Y2)
# 
# a <- 2
# b <- 1
# 
# prior_tab_low_ed <- get_poisson_prior_tab(0, 5, a, b, model_name = "BC< prior")
# prior_tab_high_ed <- get_poisson_prior_tab(0, 5, a, b, model_name = "BC+ prior")
# 
# theta <- prior_tab_high_ed$theta
# prior <- bind_rows(prior_tab_low_ed, prior_tab_high_ed)
# 
# 
# posterior_tab_low_ed <- get_poisson_posterior_tab(a, b, n1, y1, theta, "BC< posterior")
# posterior_tab_high_ed <- get_poisson_posterior_tab(a, b, n2, y2, theta, "BC+ posterior")
# 
# posterior <- bind_rows(posterior_tab_low_ed, posterior_tab_high_ed)
# 
# likelihood_tab_low_ed <- get_poisson_likelihood_tab(n1, y1, theta, 10, "BC< likelihood")
# likelihood_tab_high_ed <- get_poisson_likelihood_tab(n2, y2, theta, 10, "BC+ likelihood")
# 
# likelihood <- bind_rows(likelihood_tab_low_ed, likelihood_tab_high_ed)
# 
# full_tab <- bind_rows(prior, posterior, likelihood)
# 
# low_full_tab <- bind_rows(prior_tab_low_ed, posterior_tab_low_ed, likelihood_tab_low_ed)
# high_full_tab <- bind_rows(prior_tab_high_ed, posterior_tab_high_ed, likelihood_tab_high_ed)
# 
# 
# dist_plot <- get_distribution_plot(full_tab)
# dist_plot
# 
# characteristics_low <- get_posterior_poisson_characteristics(a, b, n1, y1, low_full_tab)
# characteristics_high <- get_posterior_poisson_characteristics(a, b, n2, y2, high_full_tab)
# 
# print(rbind(characteristics_low, characteristics_high), digits = 3)
# 
# sens_plt_low <- get_poisson_sensitivity_plot(a, b, n1, y1,  0, 4, 0, 20)
# sens_plt_high <- get_poisson_sensitivity_plot(a, b, n2, y2,  0, 4, 0, 20)
# 
# gridExtra::grid.arrange(sens_plt_low, sens_plt_high, ncol = 2)
# 
# prior_pred <- predict_neg_bin_prior(0, 8, a, b, "marginal")
# post_pred_low <- predict_neg_bin_posterior(0, 8, a, b, y1, n1, "posterior <BC")
# post_pred_high <- predict_neg_bin_posterior(0, 8, a, b, y2, n2, "posterior BC+")
# 
# pred_full <- rbind(prior_pred, post_pred_low, post_pred_high)
# 
# predict_plot <- get_predict_dist_plots(pred_full, "predikovany pocet deti")
# 
# # 2. Nemocnice ------------------------------------------------------------
# 
# data <- read.csv2(file = "srdce1.csv") |>
#   tibble::column_to_rownames(var = "nemocnice")
# print(data)





#' Calculate basal area from diameter at breast height
#'
#' @param dbh Diameter in centimeters
#' @return Basal area in square meters
calc_basal_area <- function(dbh) {
  ((dbh / 200)^2) * pi
}

#' Quadratic mean diameter
#'
#' @param dbh Numeric vector of diameters
#' @return QMD value
calc_qmd <- function(dbh) {
  sqrt(mean(dbh^2, na.rm = TRUE))
}

#' Estimate biomass based on dbh, height and merchantable volume
#'
#' Formulas adapted from compilation scripts.
#' @param species Species code
#' @param dbh Diameter at breast height
#' @param height Tree height
#' @param vol_1307 Merchantable volume (1307 definition)
#' @return Estimated biomass in kg
calc_biomass <- function(species, dbh, height, vol_1307) {
  m <- ifelse(dbh < 15 | vol_1307 == 0,
    dplyr::case_when(
      species == "AW" ~ 0.26738 + 0.01917 * dbh^2 * height,
      species == "BW" ~ 2.47035 + 0.02454 * dbh^2 * height,
      species == "FA" ~ 7.03447 + 0.01477 * dbh^2 * height,
      species %in% c("FB", "FD") ~ 7.94339 + 0.01465 * dbh^2 * height,
      species %in% c("LT", "LA", "LW") ~ 4.48372 + 0.01768 * dbh^2 * height,
      species == "PB" ~ 10.74706 + 0.01350 * dbh^2 * height,
      species == "PJ" ~ 2.91931 + 0.01678 * dbh^2 * height,
      species %in% c("PL", "PW", "PF", "PX") ~ 8.18267 + 0.01597 * dbh^2 * height,
      species == "SB" ~ 2.79552 + 0.01698 * dbh^2 * height,
      species %in% c("SW", "SE", "SX") ~ 6.03377 + 0.01500 * dbh^2 * height,
      TRUE ~ NA_real_
    ),
    dplyr::case_when(
      species == "AW" ~ 499.508 * vol_1307 ^ 0.980765,
      species == "BW" ~ 703.360 * vol_1307 ^ 0.946751,
      species == "FA" ~ 434.694 * vol_1307 ^ 0.903315,
      species %in% c("FB", "FD") ~ 444.532 * vol_1307 ^ 0.873007,
      species == "PJ" ~ 477.288 * vol_1307 ^ 0.983019,
      species %in% c("PL", "PW", "PF", "PX") ~ 436.564 * vol_1307 ^ 0.962308,
      species %in% c("LT", "LA", "LW") ~ 530.347 * vol_1307 ^ 0.922289,
      species == "SB" ~ 516.226 * vol_1307 ^ 1.001660,
      species %in% c("SW", "SE", "SX") ~ 451.544 * vol_1307 ^ 0.958852,
      TRUE ~ NA_real_
    )
  )
  m
}

# internal helper functions for taper and volume
hr_iter <- function(df, topdib, max_iter = 1000) {
  df[, diff := ifelse(dbh > topdib, 1, NA_real_)]
  df[, r0 := ifelse(dbh > topdib, 1, NA_real_)]
  df[, r1 := ifelse(dbh > topdib, 1, NA_real_)]
  hr1 <- function(b1, b2, b3, b4, b5, dbh, ht, r0) {
    b1 * r0^2 + b2 * log(r0 + 0.001) + b3 * sqrt(r0) + b4 * exp(r0) + b5 * (dbh / ht)
  }
  hr2 <- function(a0, a1, a2, b1, b2, b3, b4, b5, dbh, ht, r0, topdib) {
    (1 - ((topdib / (a0 * dbh^a1 * a2^dbh))^(1 / hr1(b1, b2, b3, b4, b5, dbh, ht, r0))) * (1 - sqrt(0.2250)))^2
  }
  n <- 0
  while (df[, max(diff, na.rm = TRUE)] > 1e-08) {
    df[dbh > topdib & diff > 1e-08, r1 := hr2(a0, a1, a2, b1, b2, b3, b4, b5, dbh, ht, r0, topdib)]
    df[dbh > topdib & diff > 1e-08, r0 := (r1 + r0) / 2]
    df[dbh > topdib & diff > 1e-08, diff := abs(r1 - r0)]
    n <- n + 1
    if (n > max_iter) {
      break
    }
  }
  invisible(NULL)
}

dib <- function(a0, a1, a2, b1, b2, b3, b4, b5, dbh, dheight, ht, r = (dheight / ht)) {
  (a0 * dbh^a1) * (a2^dbh) * ((1 - sqrt(r)) / (1 - sqrt(0.2250)))^(b1 * (r)^2 + b2 * log(r + 0.001) + b3 * sqrt(r) + b4 * exp(r) + b5 * dbh / ht)
}

segment_volume <- function(sl, stumpH, dbh, ht, a0, a1, a2, b1, b2, b3, b4, b5) {
  svol <- function(sl, d0, d1, d2) (sl / 6) * ((1 / 200)^2 * pi) * (d0^2 + 4 * d1^2 + d2^2)
  v <- data.table()
  h0 <- stumpH
  d0 <- dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, stumpH, ht)
  h1 <- h0 + sl / 2
  d1 <- dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, h1, ht)
  h2 <- h1 + sl / 2
  d2 <- dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, h2, ht)
  v[, v1 := svol(sl, d0, d1, d2)]
  for (n in 2:20) {
    h0 <- h2
    d0 <- d2
    h1 <- h0 + sl / 2
    d1 <- dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, h1, ht)
    h2 <- h1 + sl / 2
    d2 <- dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, h2, ht)
    v[, paste0("v", n) := svol(sl, d0, d1, d2)]
  }
  rowSums(v)
}

#' Calculate stem volume
#'
#' @param table Data table of tree measurements
#' @param merch Merchantability rule: "total", "1307", or "1510"
#' @param max_iter Maximum iterations for taper solution
#' @return Numeric vector of volumes (m^3)
calc_volume <- function(table, merch = "total", max_iter = 1000) {
  df <- data.table::copy(data.table::as.data.table(table))
  if (!all(c("dbh", "height", "species", "natural_subregion") %in% names(df))) {
    stop("Missing required columns: dbh, height, species, natural_subregion")
  }
  if (any(is.na(df$height))) stop("Table contains missing heights.")
  stumpH <- 0.3
  if (merch == "total") {
    topdib <- 0.0125
  } else if (merch == "1307") {
    topdib <- 7
    stumpD <- 13
    minlen <- 3.66
  } else if (merch == "1510") {
    topdib <- 10
    stumpD <- 15
    minlen <- 3.66
  } else {
    stop("Invalid merch parameter")
  }
  params <- data.table::fread("GYPSY/GYPSY data/lookup/taper.csv")
  df <- df[params, on = c(species = "species", natural_subregion = "natsub")]
  hr_iter(df, topdib, max_iter)
  df[dbh > topdib, dibs := dib(a0, a1, a2, b1, b2, b3, b4, b5, dbh, stumpH, height)]
  df[dbh > topdib, dobs := k7 + k8 * dibs]
  if (merch == "total") {
    df[, mh := height * r0]
    df[, ml := mh - stumpH]
    df[, sl := ml / 20]
    df[, vol := segment_volume(sl, stumpH, dbh, height, a0, a1, a2, b1, b2, b3, b4, b5)]
    df[, tvol := vol + pi * (topdib/200)^2 * (height - mh)/3 + pi * (dibs/200)^2 * stumpH]
    df[is.na(tvol), tvol := 0]
    df$tvol
  } else {
    df[dobs > stumpD, mh := height * r0]
    df[dobs > stumpD, ml := mh - stumpH]
    df[ml > minlen, sl := ml / 20]
    df[ml > minlen, vol := segment_volume(sl, stumpH, dbh, height, a0, a1, a2, b1, b2, b3, b4, b5)]
    df[is.na(vol), vol := 0]
    df$vol
  }
}

#' Predict diameter from height
#'
#' @param height Height in meters
#' @param species Species code
#' @param natsub Natural subregion
#' @return Predicted DBH in centimeters
predict_dbh <- function(height, species, natsub) {
  params <- data.table::fread("GYPSY/GYPSY data/lookup/ht_to_dbh_2016r1.csv")
  dt <- data.table::data.table(height = height, species = species, natsub = natsub)
  dt <- dt[params, on = c("species", "natsub")]
  dt[, b1 * (height - 1.3)^b2 * exp(-b3 * (height - 1.3))]
}

#' Predict height from diameter
#'
#' @param dbh Diameter at breast height (cm)
#' @param species Species code
#' @param natsub Natural subregion
#' @return Predicted height in meters
predict_height <- function(dbh, species, natsub) {
  params <- data.table::fread("GYPSY/GYPSY data/lookup/dbh_to_ht_2016r4.csv")
  dt <- data.table::data.table(dbh = dbh, species = species, natsub = natsub)
  dt <- dt[params, on = c("species", "natsub")]
  dt[, ifelse(model == 1,
              1.3 + b1 * (1 - exp(-1 * b2 * dbh))^b3,
              1.3 + (b1 / (1 + exp(b2 + b3 * log(dbh + 0.01)))))]
}

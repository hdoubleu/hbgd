get_int_var_in_range <- function(par_a, par_b, par_c, par_d, par_e) {
  # print(cbind(par_a, par_c, par_d, par_e))
  ret_par_d = NA
  observed = NA
  if (any(par_c == 0, na.rm = TRUE)) {
    ret_par_d = par_d[par_c == 0]
    observed = par_a[par_c == 0]
  } else {
    if ("first" == tolower(par_e)) {
      ret_par_d = head(par_d, n = 1L)
      observed = head(par_a, n = 1L)
    } else {
      if ("last" == tolower(par_e)) {
        ret_par_d = tail(par_d, n = 1L)
        observed = tail(par_a, n = 1L)
      } else {
        if ("nearest" == tolower(par_e)) {
          # par_c_par_d <- cbind.data.frame(par_c,par_d)
          a = cbind.data.frame(par_a, par_c, par_d)[order(par_c), ][1, ]
          observed = a[[1]]
          ret_par_d = a[[3]]
        } else {
          if ("weightedavg" == tolower(par_e)) {
            par_c_rec = 1/par_c
            par_c_rec_sum = sum(par_c_rec)
            ret_par_d = par_d %*% (par_c_rec/par_c_rec_sum)
            if (1 < length(par_c)) {
              observed = NA
            } else {
              observed = head(par_a, n = 1L)
            }
          } else {
            observed = NA
            ret = NA
          }
        }
      }
    }
  }
  # print(c(observed, ret_par_d))
  return(c(observed, ret_par_d))
}

#' Type checking on the column names.
#'
#' @param cols A vector of strings (vector of characters), that
#'   consists of the column names of a _long_ format data frame.
#' @return Logical FALSE if the type checking fails.
long_colnames_correct_type <- function(cols) {
  if (is.null(cols)) {
    return(FALSE)
  }
  if (!is.vector(cols)) {
    return(FALSE)
  }
  if (is.list(cols)) {
    return(FALSE)
  }
  if (0 == length(cols)) {
    return(FALSE)
  }
  if (!is.character(cols)) {
    return(FALSE)
  }
  if (anyNA(cols)) {
    return(FALSE)
  }

  return(TRUE)
}

#' Checks whether all the _time invariant_ co-variates are
#' present in the column name.
#' 
#' @param vec_timeinv Vector of strings (vector of characters)
#'   that contains the time-invariant column names.
#' @param df_colnames Vector of strings (vector of characters) that
#'   contains the column names of _long_ format data frame.
#' @return logical FALSE when invalied `vec_timeinv` is supplied or
#'   not all the variables in `vec_timeinv` are present in
#'   `df_colnames`
check_time_invar_colnames <- function(vec_timeinv, df_colnames) {
  if (missing(vec_timeinv)) {
    return(TRUE)
  } else {
    if (0 == length(vec_timeinv)) {
      return(TRUE)
    } else {
      if (!is.vector(vec_timeinv)) {
        print("Time invariant variables are not in the form of vector.")
        return(FALSE)
      } else {
        if (!is.character(vec_timeinv)) {
          print("Time invariant variables are not character strings.")
          return(FALSE)
        } else {
          return(all(vec_timeinv %in% df_colnames))
        }
      }
    }
  }
}


# Input is a list, with each element being
# [1]: specified the column name in the output (mandatory)
# [2]: the column name in the input (mandatory)
# [3]: the target (in 'agedays') (optional, default = 100)
# [4]: the tolerance (optional, default = 15)
# [5]: the strategy, one of 'first', 'last', 'nearest' and 'weightedavg' (optional, default = 'nearest')
check_one_time_variant <- function(input, df_col_names) {
  if (missing(input) || !is.list(input)) {
    print("input missing or not a list of var - target - tolerance - strategy")
    return(FALSE)
  }

  len <- length(input)
  if (2 > len || 5 < len) {
    print("input not in acceptable length")
    return(FALSE)
  }

  if (!input[[2]] %in% df_col_names) {
    print("input not in the given data frame")
    return(FALSE)
  }

  if (3 <= len && !is.numeric(input[[3]])) {
    print("input has no valid target")
    return(FALSE)
  }

  if (4 <= len && !is.numeric(input[[4]])) {
    print("input has no valid tolerance range")
    return(FALSE)
  }

  if (5 <= len && !input[[5]] %in% c("first", "last", "nearest", "weightedavg")) {
    print("input has no supported strategy")
    return(FALSE)
  }

  return(TRUE)
}


check_list_time_variant <- function(list_timev, df_colnames) {
  if (missing(list_timev)) {
    return(TRUE)
  } else {
    if (!is.list(list_timev)) {
      print("Time variant variables are not in the form of list.")
      return(FALSE)
    } else {
      if (0 == length(list_timev)) {
        return(TRUE)
      } else {
        for (timev in list_timev) {
          if (!check_one_time_variant(timev, df_colnames)) {
            return(FALSE)
          }
        }
        return(TRUE)
      }
    }
  }
}


# Inputs:
# df- Input data frame
# col - a character string being the name of the column to be filled
fill_col <- function(df, col_regexp) {
  allcols <- colnames(df)
  cols_indices <- grep(col_regexp, allcols)
  for (cindex in cols_indices) {
    a <- df[cindex]
    setorder(a, na.last = T)
    df[cindex] <- rep(a[1, 1], nrow(df))
  }
  return(df)
}

#' Converts _long_ to _wide_ format for one variable.
#' Assumes the input \code{data.frame} is already sorted as `df[with(df, order(get(key_var), get(mea_var))),]`
#'
#' @param df \code{data.frame} to convert.
#' @param oname
#' @param target_val
#' @param tol_ran
#' @param int_var
#' @param str_val String (vector of characters) indicating the name of
#'   strategy to use when converting the data. Valid strategies are
#'   'first', 'last', 'nearest' and 'weightedavg'. Default is `nearest`.
#' @param key_var String (vector of characters) indicating the subject level
#'  variable. Defaults to `SUBJID`
#' @param mea_var String (vector of characters) indicating the measurement variable.
#'   Defaults to `AGEDAYS`
#'
#' @return
long_to_wide_one_var <- function(df, oname, target_val, tol_ran, int_var, str_val = "nearest", key_var = "SUBJID",
                                 mea_var = "AGEDAYS") {
  dt1 <- setDT(df)[is.finite(get(int_var)), `:=`(offset, abs(get(mea_var) - target_val))]
  dt1 <- dt1[offset <= tol_ran, `:=`(new1, .(list(get_int_var_in_range(get(mea_var), target_val, offset,
                                                                       get(int_var), str_val)))), by = key_var]
  dt1 <- dt1[exists("new1"), `:=`(output_observed_on = sapply(new1, `[`, 1), output_value = sapply(new1,
                                                                                                   `[`, 2))]
  dt1 <- dt1[!exists("output_observed_on"), `:=`(output_observed_on = NA, output_value = NA)]
  is.na(dt1) <- dt1 == "NULL"
  dt1 <- dt1[, `:=`(output_varname = oname, input_varname = int_var, target = target_val, tolerance = tol_ran,
                    strategy = str_val)]
  return(dt1[!is.na(output_value), .SD[1], .SDcols = c("output_varname", "output_value", "output_observed_on",
                                                       "input_varname", "target", "tolerance", "strategy"), by = key_var])
}


#' Converts data from _long_ to _wide_ format
#'
#' @param df_long \code{\link{data.frame}} that contains data in _long_ format.
#'   The co-variates, which are columns, can be either _time variant_ (e.g., `MAGE`) or
#'   _time invariant_ (e.g., `HAZ`). The data.frame _must_ contain the
#'   columns `SUBJID` and `AGEDAYS`.
#' @param vector_timeinvar A vector of _time invariant_ co-variates to extract
#'    and include.
#' @param list_timevar A list of _time variant_ co-variates to extract and
#'    include.
#'    (a) a list of entities
#'    (b) each represents a request of summary of one time variant covariate; the time variant covariate needs to exist
#'    in df_long
#'    (c) each request is in the form of a list, see the documentation of check_one_time_variant() for its content.
#' @return
long_to_wide_detailed <- function(df_long, vector_timeinvar = c(), list_timevar = list()) {

  gen_log_msg <- gen_log_msg_fn("long_to_wide_detailed")

  if (missing(df_long) || !is.data.frame(df_long))
    stop(gen_log_msg("Input data not given or not in the form of data.frame"))

  df_long_col_names <- colnames(df_long)

  colname_subjid <- grep("SUBJID", df_long_col_names, ignore.case = T, value = T)
  if ( 0 == length(colname_subjid) )
    stop(gen_log_msg("Input data frame does not have the column `SUBJID`"))

  colname_agedays <- grep("AGEDAYS", df_long_col_names, ignore.case = T, value = T)
  if ( 0 == length(colname_agedays) )
    stop(gen_log_msg("Input data frame does not have the column agedays."))

  df_long <- df_long[with(df_long, order( get(colname_subjid), get(colname_agedays) ) ), ]

  if (!long_colnames_correct_type(df_long_col_names)) {
    stop(gen_log_msg("Column names of the input data frame have incorrect type."))
  }

  if (!check_time_invar_colnames(vector_timeinvar, df_long_col_names)) {
    stop(gen_log_msg("Not all requested time invariant columns are present in the input data frame."))
  }

  if (!check_list_time_variant(list_timevar, df_long_col_names)) {
    stop(gen_log_msg("Not all requested time variant columns are present in the input data frame."))
  }

  df_wide <- data.frame(sort(unique(df_long[,colname_subjid])))
  colnames(df_wide) <- colname_subjid
  num_var <- 0

  for (v in vector_timeinvar) {
    df <- long_to_wide_one_var(df_long, v, 0, Inf, v, "first", colname_subjid, colname_agedays)

    # if ( 0 == nrow(df) || 0 == ncol(df) ) next Commented out this line to make the all-NA returned data
    # frame still appear in the final outcome.

    num_var <- num_var + 1
    df_ncol <- ncol(df)
    colnames(df)[2:df_ncol] <- paste0(colnames(df)[2:df_ncol], "_", num_var)
    df_wide <- merge(setDT(df_wide), setDT(df), by = colname_subjid, all = T)
    is.na(df_wide) <- df_wide == "NULL"
  }

  for (v in list_timevar) {
    len <- length(v)

    if (5 > len) {
      strategy <- "nearest"
    } else {
      strategy <- v[[5]]
      if (4 > len) {
        tolerance <- 15
      } else {
        tolerance <- v[[4]]
        if (3 > len) {
          target <- 100
        } else {
          target <- v[[3]]
        }
      }
    }
    iname <- v[[2]]
    oname <- v[[1]]

    df <- long_to_wide_one_var(df_long, oname, target, tolerance, iname
                               , strategy, colname_subjid, colname_agedays)

    # if ( 0 == nrow(df) || 0 == ncol(df) ) next Commented out this line to make the all-NA returned data
    # frame still appear in the final outcome.

    num_var <- num_var + 1
    df_ncol <- ncol(df)
    colnames(df)[2:df_ncol] <- paste0(colnames(df)[2:df_ncol], "_", num_var)
    df_wide <- merge(setDT(df_wide), setDT(df), by = colname_subjid, all = T)
    is.na(df_wide) <- df_wide == "NULL"
  }

  df_wide <- fill_col(data.table::setDF(df_wide), "output_varname_*|input_varname_*|target_*|tolerance_*|strategy_*")

  return(df_wide)
}

#'
#' Note: This function coherces the \code{\link{data.frame}} in the parameter
#' \code{dw_detailed} into a \code{\link[data.table]{data.table}} by
#' _reference_ (as opposed to making a copy). As a result, beware when calling
#' this function in parallelised function while making changes to
#' \code{dw_detailed}.
#'
#' @param dw_detailed \code{\link{long_to_wide_detailed}} function.
long_to_wide_analysis <- function(dw_detailed) {
  df <- data.table::setDF(dw_detailed)  # Cohereces by reference.
  dw_analysis <- df[1] # Assuming "SUBJID" is on column 1.
  col_iterator <- 0
  va_index <- 0
  num_col <- ncol(dw_detailed)
  while (col_iterator < num_col) {
    va_index <- va_index + 1
    varange <- grep(paste0("_", va_index, "$"), colnames(dw_detailed))
    col_iterator <- tail(varange, 1)
    dw_analysis <- cbind(dw_analysis, df[varange[2]])
    colnames(dw_analysis)[1 + va_index] <- df[1, varange[1]]
  }
  return(dw_analysis)
}

#' Converts \code{\link{data.frame}} from _long_ format to _wide_ format.
#'
#' @param df_long
#' @param vector_timeinvar
#' @param list_timevar
#' @param detailed
#' @return
#'
#' @export
long_to_wide <- function(df_long, vector_timeinvar = c(), list_timevar = list(), detailed = FALSE) {

  dw_detailed <- long_to_wide_detailed(df_long, vector_timeinvar, list_timevar)

  if (detailed)
    return(dw_detailed)
  else
    return(long_to_wide_analysis(dw_detailed))
}


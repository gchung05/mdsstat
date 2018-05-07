#' Check Input Parameters
#'
#' Verifies correct class and, optionally, verifies existence in a data frame
#'
#' @param x Single object name, character, or vector of characters
#' @param check_class Vector of correct classes in character format
#' @param check_names Default: \code{NULL}. Optional data frame to check if
#' \code{x} exists in the data frame.
#' @param null_ok Default: \code{T}. Allow \code{x} to be \code{NULL}. If
#' \code{x} is \code{NULL}, bypass checking.
#' @param exclusions Default: \code{NULL}. Optional values that are excluded
#' from checking in \code{check_names} (if it is not \code{NULL}).
#' @param max_length Default: \code{NULL}. Optional maximum length of \code{x}.
#' If \code{NULL}, no max length check will occur.
#'
#' @return Stop error, if an error is found. Else nothing.
#' @keywords internal
input_param_checker <- function(
  x,
  check_class=NULL,
  check_names=NULL,
  null_ok=T,
  exclusions=NULL,
  max_length=NULL
){
  if (!is.null(exclusions)){
    x <- x[-which(x %in% exclusions)]
  }
  if (length(x) > 0){
    if (length(check_class) > 0){
      if (length(check_names) > 0){
        if (class(x) == "character"){
          this_class <- lapply(check_names[x], class)
        } else this_class <- class(x)
      } else this_class <- class(x)
      if (!any(sapply(this_class, function(x) x %in% check_class))){
        stop(paste0(ifelse(class(x) == "character", paste(x, collapse=" "), quote(x)),
                    " must be of class '", paste(check_class, collapse=" "), "'"))
      }
    }
    if (length(check_names) > 0){
      if (!all(x %in% names(check_names))){
        stop(paste0(x, " names not found in ", check_names))
      }
    }
    if (!is.null(max_length)){
      if (length(x) > max_length){
        stop(paste0(x, " cannot exceed maximum length of ", max_length))
      }
    }
  }
}

#' Convert to Acceptable Date
#'
#' Converts a \code{Date} vector into its equivalent daily, monthly, yearly,
#' etc... \code{Date} vector.
#'
#' @param x Input vector of class \code{Date}.
#' @param convert_type Default: \code{"months"}. String value indicating the
#' date type. Possible values are \code{"days"}, \code{"months"}.
#' @param convert_to_n Default: \code{1}. Numeric value indicating the number of
#' \code{convert_type}s to convert to.
#'
#' Example: \code{convert_type="months"} and \code{convert_to_n=3} indicates a
#' quarterly conversion.
#'
#' @return Converted \code{Date} vector of class \code{mdpms.Date} with
#' attributes \code{adder} (function that adds units),
#' \code{convert_to_n}, and \code{convert_type}.
#' @keywords internal
convert_date <- function(
  x,
  convert_type="months",
  convert_to_n=1
){
  input_param_checker(convert_type, check_class="character",
                      check_names=char_to_df(c("months", "days")), max_length=1)
  if (convert_to_n < 1 | (convert_to_n %% 1 != 0)){
    stop("convert_to_n must be positive integer")
  }
  if (convert_type == "days"){
    # Days
    this <- lubridate::floor_date(x, lubridate::days(convert_to_n))
    adder <- function(t, n){
      lubridate::ymd(t) + lubridate::ddays(convert_to_n * n)
    }
  } else if (convert_type == "months"){
    # Months
    this <- lubridate::floor_date(x, months(convert_to_n))
    adder <- function(t, n){
      lubridate::add_with_rollback(lubridate::ymd(t), months(convert_to_n * n))
    }
  }
  # Save the output class
  out <- structure(this,
                   convert_to_n=convert_to_n,
                   convert_type=convert_type,
                   adder=adder)
  class(out) <- append("mds_Date", class(out))
  return(out)
}

#' Character Vector to Header of Empty Data Frame
#' Converts a character vector into the column names of an empty data frame.
#' Used in conjunction with \code{check_names} parameter of
#' \code{input_param_checker()}
#'
#' @param x Input character vector.
#' @return Empty data frame where columns are named \code{x}.
#' @keywords internal
char_to_df <- function(
  x
){
  out <- as.data.frame(t(matrix(rep("", length(x)))), stringsAsFactors=F)
  names(out) <- x
  return(out)
}

#' Min/Max With All NA's Allowed
#' Min and Max functions that allows the vector \code{x} to be \code{NA}, where
#' \code{f} is \code{min} or \code{max}
#' @param x Numeric input vector
#' @param f Either \code{min} or \code{max}
#' @return Numeric value corresponding to \code{min} or \code{max}, or \code{NA}
#' @keywords internal
fNA <- function(x, f){
  if (all(is.na(x))){
    NA
  } else{
    f(x, na.rm=T)
  }
}

#' Return next level up device
#' Returns the variable name of the next level in the device hierarchy
#' @param x String input of device name, such as \code{"device_1"},
#' \code{"device_2"}, etc.
#' @return String representation of the next device level variable
#' @keywords internal
next_dev <- function(
  x
){
  paste0("device_", as.numeric(gsub("device_", "", x)) + 1)
}

#' Return next level up event
#' Returns the variable name of the next level in the event hierarchy
#' @param x String input of event name, such as \code{"event_1"},
#' \code{"event_2"}, etc.
#' @return String representation of the next event level variable
#' @keywords internal
next_ev <- function(
  x
){
  paste0("event_", as.numeric(gsub("event_", "", x)) + 1)
}

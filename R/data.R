#' Bone Cement MAUDE Events in 2017
#'
#' A dataset containing 535 events reported into the FDA MAUDE database on bone
#' cement in 2017. Data were obtained via the openFDA API
#' (\url{https://open.fda.gov/api/}).
#'
#' @format A data frame with 535 rows and 15 variables. Full variable
#' descriptions may be found on the FDA Device Reference Guide
#' (\url{https://open.fda.gov/device/event/reference/}). Note that \code{region}
#' is a simulated variable not present in MAUDE. Descriptions as follows:
#' \describe{
#'   \item{report_number}{Identifying number for the adverse event report.}
#'   \item{event_type}{Outcomes associated with the adverse event.}
#'   \item{date_received}{Date the report was received by the FDA.}
#'   \item{product_problem_flag}{Indicates whether or not a report was about the
#'   quality, performance or safety of a device.}
#'   \item{adverse_event_flag}{Whether the report is about an incident where the
#'   use of the device is suspected to have resulted in an adverse outcome in a
#'   patient.}
#'   \item{report_source_code}{Source of the adverse event report.}
#'   \item{lot_number}{The lot number found on the label or packaging material.}
#'   \item{model_number}{The exact model number found on the device label or
#'   accompanying packaging.}
#'   \item{manufacturer_d_name}{Device manufacturer name.}
#'   \item{manufacturer_d_country}{Device manufacturer country.}
#'   \item{brand_name}{The trade or proprietary name of the suspect medical
#'   device as used in product labeling or in the catalog.}
#'   \item{device_name}{This is the proprietary name, or trade name, of the
#'   cleared device.}
#'   \item{medical_specialty_description}{Regulation Medical Specialty is
#'   assigned based on the regulation (e.g. 21 CFR Part 888 is Orthopedic
#'   Devices).}
#'   \item{device_class}{A risk based classification system for all medical
#'   devices ((Federal Food, Drug, and Cosmetic Act, section 513)}
#'   \item{region}{A simulated, randomly assigned geographical region for
#'   package example purposes.}
#' }
#' @source \url{https://open.fda.gov/data/maude/}
"maude"

#' Simulated Device Sales Data
#'
#' A dataset containing simulated monthly sales by device and country for
#' devices reported in the \code{maude} dataset. For package usage examples,
#' this data serves as a proxy for exposures. Data were generated using a
#' random normal distribution weighted by the number of reported events by
#' device and country.
#'
#' @format A data frame with 360 rows and 4 variables:
#' \describe{
#'   \item{device_name}{Name of the device mapped from the \code{maude} dataset.}
#'   \item{region}{Geographical region mapped from the \code{maude} dataset.}
#'   \item{sales_month}{The month of sales.}
#'   \item{sales_volume}{The volume of sales.}
#' }
#' @source Random normal distribution using \code{rnorm()}.
"sales"

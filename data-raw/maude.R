# Query MAUDE
maude <- data.frame()
myapikey <- Sys.getenv("FDAapiKey")
initst <- 20170101
initend <- 20171231
searchstring <- "knee+prosthesis"
searchstring <- "bone+cement"
while(initst <= initend){
  query <- paste0(
    "https://api.fda.gov/device/event.json?",
    "api_key=", myapikey,
    "&search=device.generic_name:\"",
    searchstring,
    "\"+AND+date_received:[", initst, "+TO+", initend, "]",
    "&limit=100")
  this <- jsonlite::fromJSON(query)$results
  maude <- plyr::rbind.fill(maude, this)
  range(as.numeric(maude$date_received))
  initst <- as.numeric(max(this$date_received)) + 1
}

# Select certain columns to reduce size
cols <- c("report_number", "event_type", "date_received", "product_problem_flag",
          "adverse_event_flag", "report_source_code")
build <- maude[, cols]
cols_dev <- c("lot_number", "model_number", "manufacturer_d_name",
              "manufacturer_d_country", "brand_name")
for (i in cols_dev){
  build[[i]] <- sapply(maude$device, function(x) x[[i]][1])
}
cols_dev_openfda <- c("device_name", "medical_specialty_description",
                      "device_class")
for (i in cols_dev_openfda){
  build[[i]] <- sapply(maude$device, function(x) x$openfda[[i]][1])
}
# Simulate a reporting geographical region
build$region <- as.character(factor(as.numeric(cut(runif(nrow(build)), 3)),
                                    labels=c("West", "Central", "East")))
# Done
maude <- build

# Format as ASCII
for (i in 1:ncol(maude)){
  maude[, i] <- iconv(maude[, i], "UTF-8", "ASCII")
}

usethis::use_data(maude, overwrite=T)

# Get date range for exposure simulation, based on MAUDE
dt_range <- convert_date(range(as.Date(parsedate::parse_date(maude$date_received))))
dt_seq <- seq(dt_range[1], dt_range[2], "months")
# Simulate exposures by device_name and region
sales <- data.frame()
wts <- table(maude$device_name, maude$region)
for (i in 1:nrow(wts)){
  for (j in 1:ncol(wts)){
    ij <- data.frame(device_name=rownames(wts)[i],
                     region=colnames(wts)[j],
                     sales_month=dt_seq,
                     sales_volume=round(abs(rnorm(length(dt_seq),
                                                  1e2 * (wts[i, j] + 1),
                                                  (1e2 * (wts[i, j] + 1)) / 4))),
                     stringsAsFactors=F)
    sales <- rbind(sales, ij)
  }
}

devtools::use_data(sales, overwrite=T)

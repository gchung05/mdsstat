# Time series
Pde <- mds::deviceevent(
  maude,
  time="date_received",
  device_hierarchy=c("device_name", "device_class"),
  event_hierarchy=c("event_type", "medical_specialty_description"),
  key="report_number",
  covariates="region",
  descriptors="_all_")
Pexp <- mds::exposure(
  sales,
  time="sales_month",
  device_hierarchy="device_name",
  match_levels="region",
  count="sales_volume"
)
Pda <- mds::define_analyses(
  Pde, "device_name",
  exposure=Pexp,
  covariates="region")
Pts <- mds::time_series(Pda[sample(c(1:length(Pda)), 3)], Pde, Pexp)
mds_ts <- Pts

devtools::use_data(mds_ts, overwrite=T)

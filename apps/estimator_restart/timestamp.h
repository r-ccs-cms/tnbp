uint64_t make_timestamp() {
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);
  std::tm tm = *std::localtime(&in_time_t);
  auto micros = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()) % 1000000;
  uint64_t stamp =
    (static_cast<int64_t>(tm.tm_year + 1900) * 1000000000000000LL) +
    (static_cast<int64_t>(tm.tm_mon + 1)    * 10000000000000LL) +
    (static_cast<int64_t>(tm.tm_mday)       * 100000000000LL) +
    (static_cast<int64_t>(tm.tm_hour)       * 1000000000LL) +
    (static_cast<int64_t>(tm.tm_min)        * 10000000LL) +
    (static_cast<int64_t>(tm.tm_sec)        * 100000LL) +
    micros.count();
  return stamp;
}

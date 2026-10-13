// This file is a part of r-ccs-cms/tnbp
/**
 * @file tnbp/framework/timestamp.h
 * @brief Numeric local-time timestamps for application logs.
 */
#pragma once
#include <chrono>
#include <cstdint>
#include <ctime>
#include <stdexcept>

namespace tnbp {

/**
 * @brief Return the current local time as YYMMDDhhmmssffffff.
 * @return A decimal timestamp with six microsecond digits below the seconds.
 *
 * The two-digit year is modulo 100. This is a wall-clock log label, not a
 * monotonic timer or a unique identifier. Integer output omits leading zeros.
 * @throws std::runtime_error If conversion to local calendar time fails.
 */
inline std::uint64_t make_timestamp() {
  const auto now = std::chrono::system_clock::now();
  const auto seconds = std::chrono::floor<std::chrono::seconds>(now);
  const auto calendar_time = std::chrono::system_clock::to_time_t(seconds);
  const auto* local = std::localtime(&calendar_time);
  if(!local) throw std::runtime_error("cannot convert timestamp to local time");
  const std::tm tm = *local;
  const auto micros = std::chrono::duration_cast<std::chrono::microseconds>(now-seconds).count();
  std::uint64_t stamp = static_cast<std::uint64_t>((tm.tm_year + 1900) % 100);
  stamp = stamp * 100 + static_cast<std::uint64_t>(tm.tm_mon + 1);
  stamp = stamp * 100 + static_cast<std::uint64_t>(tm.tm_mday);
  stamp = stamp * 100 + static_cast<std::uint64_t>(tm.tm_hour);
  stamp = stamp * 100 + static_cast<std::uint64_t>(tm.tm_min);
  stamp = stamp * 100 + static_cast<std::uint64_t>(tm.tm_sec);
  return stamp * 1000000 + static_cast<std::uint64_t>(micros);
}

} // namespace tnbp

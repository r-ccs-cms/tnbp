#pragma once
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace tcapi::detail::io {
inline std::string path_from(std::string_view value) {
    if (value.empty() || value.find('\0') != std::string_view::npos)
        throw std::invalid_argument("min-tcapi: empty path or embedded NUL");
    return std::string(value);
}
inline std::string path_from(const std::string& value) {
    return path_from(std::string_view(value));
}
inline std::string path_from(const char* value) {
    if (!value) throw std::invalid_argument("min-tcapi: null path");
    return path_from(std::string_view(value));
}
// Optional storage types specialize this adapter without changing core I/O.
template<class T> struct path_adapter;
template<class S>
auto path_from(const S& value) -> decltype(path_adapter<std::decay_t<S>>::convert(value)) {
    return path_adapter<std::decay_t<S>>::convert(value);
}
} // namespace tcapi::detail::io

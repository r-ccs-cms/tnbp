#pragma once
#include "io_path.h"
#include "verbose.h"
#include <filesystem>

namespace tcapi::detail::io {
template<> struct path_adapter<std::filesystem::path> {
    static std::string convert(const std::filesystem::path& value) {
        return path_from(value.native());
    }
};
} // namespace tcapi::detail::io

namespace tcapi::detail::verbose {
template<> struct formatter<std::filesystem::path> {
    static void write(std::ostream& os, const std::filesystem::path& value) {
        quoted(os, value.string());
    }
};
} // namespace tcapi::detail::verbose

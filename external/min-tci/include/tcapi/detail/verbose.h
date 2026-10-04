#pragma once
#include "../types.h"
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <filesystem>
#include <iomanip>
#include <locale>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <unordered_map>
#include <vector>

namespace tcapi::detail::verbose {
// One initialization across translation units. Environment changes afterwards
// are deliberately ignored; launch a new process to select another level.
inline int level() noexcept {
    static const int value = [] {
        const char* text = std::getenv("TCAPI_VERBOSE");
        if (text && text[0] && !text[1]) {
            if (text[0] == '1') return 1;
            if (text[0] == '2') return 2;
        }
        return 0;
    }();
    return value;
}
inline thread_local unsigned depth = 0;
inline std::mutex output_mutex;

template<class T> constexpr const char* scalar_name() noexcept {
    if constexpr (std::is_same_v<T, float>) return "float32";
    else if constexpr (std::is_same_v<T, double>) return "float64";
    else if constexpr (std::is_same_v<T, std::complex<float>>) return "complex64";
    else if constexpr (std::is_same_v<T, std::complex<double>>) return "complex128";
    else return "unknown";
}
inline void quoted(std::ostream& os, std::string_view text) {
    constexpr char hex[] = "0123456789abcdef";
    os << '"';
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') os << '\\' << static_cast<char>(c);
        else if (c < 32 || c == 127) os << "\\x" << hex[c >> 4] << hex[c & 15];
        else os << static_cast<char>(c);
    }
    os << '"';
}
template<class T, class Enable = void> struct formatter {
    static void write(std::ostream& os, const T&) { os << "<opaque>"; }
};
template<class T> struct formatter<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static void write(std::ostream& os, T value) { os << value; }
};
template<class T> struct formatter<std::complex<T>> {
    static void write(std::ostream& os, const std::complex<T>& value) { os << value; }
};
template<> struct formatter<std::string> {
    static void write(std::ostream& os, const std::string& value) { quoted(os, value); }
};
template<> struct formatter<std::string_view> {
    static void write(std::ostream& os, std::string_view value) { quoted(os, value); }
};
template<> struct formatter<const char*> {
    static void write(std::ostream& os, const char* value) {
        if (value) quoted(os, value); else os << "null";
    }
};
template<> struct formatter<char*> : formatter<const char*> {};
template<class T> void write(std::ostream& os, const T& value) {
    formatter<std::decay_t<decltype(value)>>::write(os, value);
}
template<class A, class B> struct formatter<std::pair<A, B>> {
    static void write(std::ostream& os, const std::pair<A, B>& value) {
        os << '('; verbose::write(os, value.first); os << ',';
        verbose::write(os, value.second); os << ')';
    }
};
template<class T> struct formatter<std::reference_wrapper<T>> {
    static void write(std::ostream& os, const std::reference_wrapper<T>& value) {
        verbose::write(os, value.get());
    }
};
template<class Range> void sequence(std::ostream& os, const Range& values) {
    os << '[';
    bool first = true;
    for (const auto& value : values) {
        if (!first) os << ',';
        first = false;
        verbose::write(os, value);
    }
    os << ']';
}
template<class T, class A> struct formatter<std::vector<T, A>> {
    static void write(std::ostream& os, const std::vector<T, A>& value) { sequence(os, value); }
};
template<class K, class V, class C, class A> struct formatter<std::map<K, V, C, A>> {
    static void write(std::ostream& os, const std::map<K, V, C, A>& value) { sequence(os, value); }
};
template<class K, class V, class H, class E, class A>
struct formatter<std::unordered_map<K, V, H, E, A>> {
    static void write(std::ostream& os, const std::unordered_map<K, V, H, E, A>& value) {
        sequence(os, value);
    }
};
template<> struct formatter<std::filesystem::path> {
    static void write(std::ostream& os, const std::filesystem::path& value) {
        quoted(os, value.string());
    }
};
template<class T> void field(std::ostream& os, const char* name, const T& value) {
    os << ' ' << name << '=';
    verbose::write(os, value);
}


// Tensor metadata only; even allocated-but-uninitialized tensors are safe.
template<class E> struct formatter<gqten::tensor<E>> {
    static void write(std::ostream& os, const gqten::tensor<E>& a) {
        os << "{dtype=" << scalar_name<E>() << ",shape=";
        sequence(os,a.Shape());os << '}';
    }
};

// Declared first in each public API: metadata is captured before mutation,
// and the destructor runs after the API's local resources have been released.
// Metadata only: no tensor elements are read or numerical operations performed.
class call {
    using clock = std::chrono::steady_clock;
    int level_;
    bool outer_ = false;
    int exceptions_ = 0;
    const char* name_;
    std::string inputs_;
    clock::time_point start_{};
public:
    template<class Describe> call(const char* name, Describe&& describe) noexcept
        : level_(level()), name_(name) {
        if (!level_) return;
        outer_ = depth++ == 0;
        if (!outer_) return;
        exceptions_ = std::uncaught_exceptions();
        try {
            std::ostringstream os;
            os.imbue(std::locale::classic());
            os << std::setprecision(17);
            describe(os);
            inputs_ = os.str();
        } catch (...) {
            // Logging must not change computation or mask its exception.
        }
        if (level_ == 2) start_ = clock::now();
    }
    ~call() noexcept {
        if (!level_) return;
        --depth;
        if (!outer_) return;
        const auto stop = level_ == 2 ? clock::now() : clock::time_point{};
        try {
            std::ostringstream os;
            os.imbue(std::locale::classic());
            os << "TCAPI " << name_ << inputs_
               << " status=" << (std::uncaught_exceptions() > exceptions_ ? "exception" : "ok");
            if (level_ == 2)
                os << " elapsed_ms=" << std::setprecision(9)
                   << std::chrono::duration<double, std::milli>(stop - start_).count();
            os << '\n';
            const auto line = os.str();
            std::lock_guard<std::mutex> lock(output_mutex);
            (void)std::fwrite(line.data(), 1, line.size(), stderr);
        } catch (...) {
            // Best-effort diagnostics, including during stack unwinding.
        }
    }
    call(const call&) = delete;
    call& operator=(const call&) = delete;
};

// User callbacks form a new public-call boundary. Restore nesting even when
// a callback throws. Internal calls made by TCAPI itself remain suppressed.
class callback_boundary {
    unsigned saved_;
public:
    callback_boundary() noexcept : saved_(depth) { depth = 0; }
    ~callback_boundary() noexcept { depth = saved_; }
    callback_boundary(const callback_boundary&) = delete;
    callback_boundary& operator=(const callback_boundary&) = delete;
};
template<class F, class... Args> decltype(auto) invoke_callback(F&& f, Args&&... args) {
    if (!level()) return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    callback_boundary boundary;
    return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
}
} // namespace tcapi::detail::verbose

#include <fstream> // Some libc++ releases include filesystem internally.
#define TCAPI_TEST_NO_FILESYSTEM 1
#include "../common/check.h"
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
struct temporary_directory {
    std::string path;
    temporary_directory() {
        char pattern[] = "/tmp/min-tcapi-string-XXXXXX";
        const auto result = ::mkdtemp(pattern);
        if (!result) throw std::runtime_error("mkdtemp failed");
        path = result;
    }
    ~temporary_directory() {
        std::remove((path + "/tensor.bin").c_str());
        ::rmdir(path.c_str());
    }
};
template<class E> void test() {
    using T = gqten::tensor<E>; using R = tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    temporary_directory dir;
    const std::string name = dir.path + "/tensor.bin";
    auto a = tcapi::fill<T>(ctx, {2,3}, value<E>(3,1));
    auto roundtrip = [&](auto&& storage) {
        tcapi::save(ctx, a, storage);
        check(tcapi::close(ctx, a, tcapi::load<T>(ctx, storage), R(0)), "string path round trip");
    };
    roundtrip(name); roundtrip(name.c_str());
    const std::string suffixed = name + "ignored";
    roundtrip(std::string_view(suffixed.data(), name.size()));
    char writable[256]; std::snprintf(writable, sizeof(writable), "%s", name.c_str());
    roundtrip(writable);
    std::stringstream stream(std::ios::in|std::ios::out|std::ios::binary);
    tcapi::save(ctx, a, stream);
    check(tcapi::close(ctx, a, tcapi::load<T>(ctx, stream), R(0)), "stream round trip");
    auto bad_path = [&](auto&& storage) {
        throws<std::invalid_argument>([&] { tcapi::save(ctx, a, storage); });
        throws<std::invalid_argument>([&] { tcapi::load<T>(ctx, storage); });
    };
    bad_path(""); bad_path(static_cast<const char*>(nullptr));
    bad_path(std::string("a\0b",3)); bad_path(std::string_view("a\0b",3));
    tcapi::destroy_context(ctx);
}
RUN_FOUR_TYPES(test)

#include "tcapi/input_and_output_filesystem.h"
#include "tcapi/tcapi.h"
#include "../common/check.h"
#include <sstream>
void check_filesystem_header_before() {
    const std::filesystem::path path("tensor.bin");
    check(tcapi::detail::io::path_from(path)=="tensor.bin", "filesystem adapter before tcapi.h");
    std::ostringstream log;
    tcapi::detail::verbose::field(log, "path", path);
    check(log.str()==" path=\"tensor.bin\"", "filesystem verbose formatter");
}

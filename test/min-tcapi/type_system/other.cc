#include "tcapi/tcapi.h"
int other_translation_unit() {
    tcapi::gqten_handle ctx;
    tcapi::create_context(ctx);
    auto a = tcapi::eye<gqten::tensor<double>>(ctx, 1);
    return static_cast<int>(tcapi::get_elem(ctx, a, {0,0}));
}

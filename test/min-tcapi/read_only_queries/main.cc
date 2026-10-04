#include "../common/check.h"
template<class E> void test() {
    using T = gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    for (const auto& dims : {tcapi::shape_t<T>{}, {1}, {2,3}}) {
        auto a = tcapi::fill<T>(ctx, dims, value<E>(2, 3));
        const std::size_t count = dims.size() == 2 ? 6 : 1;
        check(tcapi::shape(ctx,a) == dims, "shape");
        check(tcapi::order(ctx,a) == static_cast<int>(dims.size()), "order");
        check(tcapi::size(ctx,a) == count, "logical size");
        check(tcapi::size_bytes(ctx,a) == count*sizeof(E), "numeric storage bytes");
        tcapi::elem_coors_t<T> coors(dims.size(),0);
        check(tcapi::get_elem(ctx,a,coors) == value<E>(2,3), "logical value");
    }
    auto a = tcapi::fill<T>(ctx,{2,3},value<E>(2,1));
    a.Scale(E(2));
    check(tcapi::get_elem(ctx,a,{1,2}) == value<E>(4,2), "lazy scale read");
    throws<std::invalid_argument>([&] { tcapi::get_elem(ctx,a,{0}); });
    throws<std::out_of_range>([&] { tcapi::get_elem(ctx,a,{-1,0}); });
    throws<std::out_of_range>([&] { tcapi::get_elem(ctx,a,{0,3}); });
}
RUN_FOUR_TYPES(test)

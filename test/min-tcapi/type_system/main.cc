#include "../common/check.h"
#include "tci/tci.h" // The two public entry points must coexist during migration.
int other_translation_unit();
template<class E> void test() {
    using T = gqten::tensor<E>;
    static_assert(std::is_same_v<tcapi::elem_t<T>, E>);
    static_assert(std::is_same_v<tcapi::ten_t<T>, T>);
    static_assert(std::is_same_v<tcapi::shape_t<T>, tcapi::List<tcapi::bond_dim_t<T>>>);
    static_assert(std::is_same_v<tcapi::elem_coors_t<T>, tcapi::List<tcapi::elem_coor_t<T>>>);
    static_assert(std::is_same_v<tcapi::CRef<T>, std::reference_wrapper<const T>>);
    static_assert(std::is_same_v<tcapi::Map<int,int>, std::unordered_map<int,int>>);
    static_assert(std::is_same_v<tcapi::Pair<int,int>, std::pair<int,int>>);
    static_assert(std::is_same_v<tcapi::real_ten_t<T>, gqten::tensor<tcapi::real_t<T>>>);
    static_assert(std::is_same_v<tcapi::cplx_ten_t<T>, gqten::tensor<tcapi::cplx_t<T>>>);
    tcapi::context_handle_t<T> ctx;
    throws<std::logic_error>([&] { (void)tcapi::zeros<T>(ctx, {}); });
    tcapi::create_context(ctx);
    auto a = tcapi::fill<T>(ctx, {}, value<E>(3, 2));
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&] { (void)tcapi::get_elem(ctx, a, {}); });
    tcapi::create_context(ctx);
    check(tcapi::get_elem(ctx, a, {}) == value<E>(3,2), "context reinitialization");
    check(other_translation_unit() == 1, "public header multi-TU linkage");
}
RUN_FOUR_TYPES(test)

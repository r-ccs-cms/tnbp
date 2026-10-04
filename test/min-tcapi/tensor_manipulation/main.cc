#include "../common/check.h"
#include "shapes.h"
#include "components.h"
#include "traversal.h"
#include "regions.h"
#include "joins.h"
template<class E> void test() {
    test_regions<E>();
    test_joins<E>();
    test_shapes<E>();
    test_components<E>();
    test_traversal<E>();
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    auto a=tcapi::fill<T>(ctx,{2,3},value<E>(3,1));
    a.Scale(E(0));
    tcapi::set_elem(ctx,a,{1,2},value<E>(5,2));
    check(tcapi::get_elem(ctx,a,{1,2})==value<E>(5,2),"write through zero lazy scale");
    check(tcapi::get_elem(ctx,a,{0,0})==E(0),"other zero elements preserved");
    a.Scale(E(2));
    tcapi::set_elem(ctx,a,{0,0},value<E>(2,4));
    check(tcapi::get_elem(ctx,a,{0,0})==value<E>(2,4),"scaled write");
    check(tcapi::get_elem(ctx,a,{1,2})==value<E>(10,4),"other scaled values preserved");
    throws<std::out_of_range>([&] { tcapi::set_elem(ctx,a,{2,0},E(1)); });
    throws<std::invalid_argument>([&] { tcapi::set_elem(ctx,a,{},E(1)); });
}
RUN_FOUR_TYPES(test)

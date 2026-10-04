#include "../common/check.h"
template<class E> void test() {
    using T = gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    for (const auto& dims : {tcapi::shape_t<T>{}, {1}, {2,3}}) {
        auto z = tcapi::zeros<T>(ctx,dims);
        check(tcapi::get_elem(ctx,z,tcapi::elem_coors_t<T>(dims.size(),0)) == E{}, "zeros");
        auto a = tcapi::allocate<T>(ctx,dims);
        tcapi::set_elem(ctx,a,tcapi::elem_coors_t<T>(dims.size(),0),value<E>(7,2));
        check(tcapi::get_elem(ctx,a,tcapi::elem_coors_t<T>(dims.size(),0)) == value<E>(7,2), "allocate then set");
        int calls = 0;
        auto gen = [&] { ++calls; return value<E>(calls, -calls); };
        auto r = tcapi::random<T>(ctx,dims,gen);
        check(calls == static_cast<int>(tcapi::size(ctx,r)), "generator call count");
        check(tcapi::get_elem(ctx,r,tcapi::elem_coors_t<T>(dims.size(),0)) == value<E>(1,-1), "generator values unchanged");
        auto duplicate = tcapi::copy(ctx,r);
        tcapi::set_elem(ctx,duplicate,tcapi::elem_coors_t<T>(dims.size(),0),E(99));
        check(tcapi::get_elem(ctx,r,tcapi::elem_coors_t<T>(dims.size(),0)) == value<E>(1,-1), "deep copy independence");
        const auto* before = r.GetRaw();
        auto moved = tcapi::move(ctx,r);
        check(moved.GetRaw() == before, "move transfers storage");
        check(tcapi::shape(ctx,r).empty() && r.GetRaw() == nullptr, "move restores default source");
        check(tcapi::get_elem(ctx,r,{}) == E(1), "default gqten scalar");
        tcapi::clear(ctx,moved);
        check(moved.GetRaw() == nullptr && tcapi::shape(ctx,moved).empty(), "clear releases dense storage");
        tcapi::clear(ctx,moved);
    }
    auto identity = tcapi::eye<T>(ctx,3);
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        check(tcapi::get_elem(ctx,identity,{i,j}) == E(i==j), "identity");
    auto scaled = tcapi::fill<T>(ctx,{2},value<E>(1,2));
    scaled.SetLabels({7});
    scaled.Scale(E(3));
    auto duplicate = tcapi::copy(ctx,scaled);
    check(tcapi::get_elem(ctx,duplicate,{1}) == value<E>(3,6), "copy lazy scale");
    check(duplicate.GetLabels() == scaled.GetLabels(), "copy backend metadata");
    auto uninitialized = tcapi::allocate<T>(ctx,{2});
    auto uninitialized_copy = tcapi::copy(ctx,uninitialized);
    tcapi::set_elem(ctx,uninitialized_copy,{0},E(3));
    check(tcapi::get_elem(ctx,uninitialized_copy,{0}) == E(3), "copy then initialize");
    const auto max = std::numeric_limits<tcapi::bond_dim_t<T>>::max();
    throws<std::invalid_argument>([&] { tcapi::zeros<T>(ctx,{0,3}); });
    throws<std::invalid_argument>([&] { tcapi::allocate<T>(ctx,{-1}); });
    throws<std::invalid_argument>([&] { tcapi::eye<T>(ctx,0); });
    throws<std::overflow_error>([&] { tcapi::allocate<T>(ctx,{max,max,max}); });
    throws<std::overflow_error>([&] { tcapi::allocate<T>(ctx,{max,max,2}); });
    int calls=0;
    auto gen=[&]() -> E { if (++calls==3) throw std::runtime_error("generator failure"); return E(1); };
    throws<std::runtime_error>([&] { tcapi::random<T>(ctx,{4},gen); });
}
RUN_FOUR_TYPES(test)

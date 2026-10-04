#pragma once
#include "../common/check.h"
#include <memory>
#include <algorithm>
template<class E> void test_traversal() {
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    for(const auto& dims : {tcapi::shape_t<T>{},{2,3}}) {
        auto a=tcapi::fill<T>(ctx,dims,value<E>(2,3));
        a.Scale(E(2));
        const int count=dims.empty()?1:6;
        int calls=0;
        E sum{};
        const T& input=a;
        const auto* before=a.GetRaw();
        auto readonly=[&](const E& e) { ++calls; sum+=e; };
        static_assert(std::is_void_v<decltype(tcapi::for_each(ctx,input,readonly))>);
        tcapi::for_each(ctx,input,readonly);
        check(calls==count && sum==E(count)*value<E>(4,6),"const traversal count and logical sum");
        check(a.GetRaw()==before,"const traversal does not replace tensor");
        tcapi::for_each(ctx,a,[state=std::make_unique<int>(1)](E& e) { e+=E(*state); });
        calls=0;
        tcapi::for_each(ctx,input,[&](const auto& e) {
            static_assert(std::is_const_v<std::remove_reference_t<decltype(e)>>);
            ++calls; check(e==value<E>(5,6),"mutable traversal result");
        });
        check(calls==count,"const generic callback count");
        std::vector<int> visits(count,0);
        tcapi::for_each_with_coors(ctx,a,[&](E& e,const auto& c) {
            static_assert(std::is_const_v<std::remove_reference_t<decltype(c)>>);
            const int index=c.empty()?0:c[0]+2*c[1];
            ++visits[index]; e+=E(index);
        });
        for(int v:visits) check(v==1,"coordinate traversal exactly once");
        std::fill(visits.begin(),visits.end(),0);
        auto coordinate_read=[&](const auto& e,const auto& c) {
            static_assert(std::is_const_v<std::remove_reference_t<decltype(e)>>);
            static_assert(std::is_const_v<std::remove_reference_t<decltype(c)>>);
            const int index=c.empty()?0:c[0]+2*c[1];
            ++visits[index];check(e==value<E>(5+index,6),"coordinate traversal values");
        };
        static_assert(std::is_void_v<decltype(tcapi::for_each_with_coors(ctx,input,coordinate_read))>);
        tcapi::for_each_with_coors(ctx,input,coordinate_read);
        for(int v:visits) check(v==1,"const coordinate traversal exactly once");
        if(dims.empty()) check(a.GetRaw()==nullptr,"scalar traversal remains inline");

        auto saved=tcapi::copy(ctx,a);
        int step=0;
        throws<std::runtime_error>([&] {
            tcapi::for_each(ctx,a,[&](E& e) { e=E(100); if(++step==count) throw std::runtime_error("callback"); });
        });
        step=0;
        throws<std::runtime_error>([&] {
            tcapi::for_each_with_coors(ctx,a,[&](E& e,const auto&) { e=E(200); if(++step==count) throw std::runtime_error("callback"); });
        });
        tcapi::for_each_with_coors(ctx,input,[&](const E& e,const auto& c) {
            check(e==tcapi::get_elem(ctx,saved,c),"throwing traversal preserves input");
        });
        a.Scale(E(0));
        tcapi::for_each(ctx,a,[](E& e) { e+=value<E>(3,2); });
        tcapi::for_each(ctx,input,[](const E& e) { check(e==value<E>(3,2),"zero-scale mutable traversal"); });
    }
    auto a=tcapi::zeros<T>(ctx,{2});
    a.SetLabels({7});
    tcapi::for_each(ctx,a,[](E& e) { e=E(1); });
    tcapi::for_each_with_coors(ctx,a,[](E& e,const auto&) { e+=E(1); });
    check(a.GetLabels()==std::decay_t<decltype(a.GetLabels())>({7}),"traversal retains labels");
    tcapi::destroy_context(ctx);
    const T& input=a;
    throws<std::logic_error>([&] { tcapi::for_each(ctx,a,[](E&) {}); });
    throws<std::logic_error>([&] { tcapi::for_each(ctx,input,[](const E&) {}); });
    throws<std::logic_error>([&] { tcapi::for_each_with_coors(ctx,a,[](E&,const auto&) {}); });
    throws<std::logic_error>([&] { tcapi::for_each_with_coors(ctx,input,[](const E&,const auto&) {}); });
}

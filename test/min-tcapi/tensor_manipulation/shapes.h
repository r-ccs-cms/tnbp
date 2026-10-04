#pragma once
#include "../common/check.h"
template<class E> void test_shapes() {
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    std::vector<E> values;
    for (int i=0; i<24; ++i) values.push_back(value<E>(i+1, -i-2));
    auto a=tcapi::assign_from_range<T>(ctx,{2,3,4},values.begin(),[](const auto& c) { return c[0]+2*c[1]+6*c[2]; });
    a.Scale(value<E>(2,1));
    auto original=tcapi::copy(ctx,a);
    const auto* storage=a.GetRaw();
    tcapi::reshape(ctx,a,{6,4});
    check(a.GetRaw()==storage,"dense reshape reuses payload");
    for(int j=0;j<4;++j) for(int i=0;i<6;++i)
        check(tcapi::get_elem(ctx,a,{i,j})==value<E>(2,1)*values[i+6*j],"reshape preserves column-major values");
    T out;
    tcapi::reshape(ctx,a,{4,6},out);
    tcapi::set_elem(ctx,out,{0,0},E(99));
    check(tcapi::get_elem(ctx,a,{0,0})==value<E>(2,1)*values[0],"out reshape independent");
    tcapi::reshape(ctx,a,{2,3,4},a);
    check(tcapi::shape(ctx,a)==tcapi::shape_t<T>({2,3,4}),"alias reshape shape");
    throws<std::invalid_argument>([&] { tcapi::reshape(ctx,a,{23}); });
    throws<std::invalid_argument>([&] { tcapi::reshape(ctx,a,{0,24}); });
    throws<std::invalid_argument>([&] { tcapi::reshape(ctx,a,{-2,-12},out); });
    check(tcapi::get_elem(ctx,out,{0,0})==E(99),"invalid out reshape preserves output");
    check(tcapi::shape(ctx,a)==tcapi::shape_t<T>({2,3,4}),"invalid in-place reshape preserves shape");

    a.SetLabels({10,20,30});
    tcapi::transpose(ctx,a,{2,0,1},out);
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({4,2,3}),"transpose shape");
    check(out.GetLabels()==std::decay_t<decltype(out.GetLabels())>({30,10,20}),"transpose labels");
    for(int k=0;k<4;++k) for(int j=0;j<3;++j) for(int i=0;i<2;++i)
        check(tcapi::get_elem(ctx,out,{k,i,j})==tcapi::get_elem(ctx,original,{i,j,k}),"transpose coordinates");
    check(tcapi::shape(ctx,a)==tcapi::shape_t<T>({2,3,4}),"out transpose leaves input unchanged");
    tcapi::transpose(ctx,out,{1,2,0});
    tcapi::transpose(ctx,a,{2,0,1},a);
    tcapi::transpose(ctx,a,{1,2,0});
    for(int k=0;k<4;++k) for(int j=0;j<3;++j) for(int i=0;i<2;++i) {
        check(tcapi::get_elem(ctx,out,{i,j,k})==tcapi::get_elem(ctx,original,{i,j,k}),"transpose inverse");
        check(tcapi::get_elem(ctx,a,{i,j,k})==tcapi::get_elem(ctx,original,{i,j,k}),"alias transpose inverse");
    }
    throws<std::invalid_argument>([&] { tcapi::transpose(ctx,a,{0,0,2}); });
    throws<std::invalid_argument>([&] { tcapi::transpose(ctx,a,{0,1}); });
    throws<std::out_of_range>([&] { tcapi::transpose(ctx,a,{-1,1,2}); });
    throws<std::out_of_range>([&] { tcapi::transpose(ctx,a,{0,1,3},out); });
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({2,3,4}),"invalid transpose preserves output");

    auto scalar=tcapi::fill<T>(ctx,{},value<E>(7,-3));
    for(const auto& dims : {tcapi::shape_t<T>{1},{1,1},{}}) {
        tcapi::reshape(ctx,scalar,dims);
        check(tcapi::size(ctx,scalar)==1,"reshape scalar logical size");
        check(tcapi::get_elem(ctx,scalar,tcapi::elem_coors_t<T>(dims.size(),0))==value<E>(7,-3),"reshape scalar value");
    }
    check(scalar.GetRaw()==nullptr,"dense-to-scalar reshape releases payload");
    tcapi::reshape(ctx,scalar,{},out);
    check(out.GetRaw()==nullptr && tcapi::get_elem(ctx,out,{})==value<E>(7,-3),"scalar replaces dense destination");
    tcapi::transpose(ctx,scalar,{},scalar);
    check(tcapi::get_elem(ctx,scalar,{})==value<E>(7,-3),"scalar alias transpose");
    tcapi::transpose(ctx,scalar,{},out);
    throws<std::invalid_argument>([&] { tcapi::transpose(ctx,scalar,{0}); });
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&] { tcapi::reshape(ctx,scalar,{}); });
    throws<std::logic_error>([&] { tcapi::transpose(ctx,scalar,{}); });
}

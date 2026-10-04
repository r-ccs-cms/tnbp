#pragma once
#include "../common/check.h"
template<class E> void test_regions() {
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx; tcapi::create_context(ctx);
    auto a=tcapi::zeros<T>(ctx,{2,3});
    for(int j=0;j<3;++j) for(int i=0;i<2;++i) tcapi::set_elem(ctx,a,{i,j},value<E>(1+i+2*j,2-i+j));
    a.Scale(value<E>(2,1));
    T out;
    tcapi::expand(ctx,a,{{0,1},{1,2}},out);
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({3,5}),"expanded shape");
    for(int j=0;j<5;++j) for(int i=0;i<3;++i)
        check(tcapi::get_elem(ctx,out,{i,j})==(i<2&&j<3?tcapi::get_elem(ctx,a,{i,j}):E{}),"expanded region and padding");
    tcapi::shrink(ctx,out,{{0,{0,2}},{1,{0,3}}});
    for(int j=0;j<3;++j) for(int i=0;i<2;++i) check(tcapi::get_elem(ctx,out,{i,j})==tcapi::get_elem(ctx,a,{i,j}),"expand/shrink inverse");
    tcapi::shrink(ctx,a,{{1,{1,3}}},out);
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({2,2}),"shrink retains unlisted axes");
    for(int j=0;j<2;++j) for(int i=0;i<2;++i) check(tcapi::get_elem(ctx,out,{i,j})==tcapi::get_elem(ctx,a,{i,j+1}),"slice offset");
    tcapi::extract_sub(ctx,a,{{1,2},{0,2}},out);
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({1,2}),"extract shape");
    check(tcapi::get_elem(ctx,out,{0,1})==tcapi::get_elem(ctx,a,{1,1}),"extract coordinates");
    tcapi::extract_sub(ctx,out,{{0,1},{1,2}},out);
    check(tcapi::get_elem(ctx,out,{0,0})==tcapi::get_elem(ctx,a,{1,1}),"extract alias");
    auto sub=tcapi::fill<T>(ctx,{1,2},value<E>(9,-2)); sub.Scale(E(2));
    tcapi::replace_sub(ctx,a,sub,{1,1},out);
    for(int j=0;j<3;++j) for(int i=0;i<2;++i)
        check(tcapi::get_elem(ctx,out,{i,j})==(i==1&&j>=1?value<E>(18,-4):tcapi::get_elem(ctx,a,{i,j})),"replacement region");
    tcapi::replace_sub(ctx,a,sub,{1,1},sub); // output aliases the subregion input
    check(tcapi::shape(ctx,sub)==tcapi::shape(ctx,a),"replacement aliases sub output");
    check(tcapi::get_elem(ctx,sub,{1,2})==value<E>(18,-4),"aliased sub value");
    tcapi::replace_sub(ctx,sub,sub,{0,0});
    check(tcapi::get_elem(ctx,sub,{1,2})==value<E>(18,-4),"self replacement");
    auto small=tcapi::fill<T>(ctx,{1,1},value<E>(7,3));
    tcapi::replace_sub(ctx,sub,small,{0,0});
    check(tcapi::get_elem(ctx,sub,{0,0})==value<E>(7,3),"in-place replacement");
    auto alias=tcapi::copy(ctx,a);
    tcapi::expand(ctx,alias,{{0,1}},alias);
    tcapi::shrink(ctx,alias,{{0,{0,2}}},alias);
    tcapi::extract_sub(ctx,alias,{{0,2},{0,3}});
    check(tcapi::get_elem(ctx,alias,{1,2})==tcapi::get_elem(ctx,a,{1,2}),"region alias round trip");
    tcapi::expand(ctx,alias,{{1,0}});
    check(tcapi::shape(ctx,alias)==tcapi::shape(ctx,a),"zero increment allowed");
    tcapi::set_elem(ctx,alias,{0,0},E(99));
    check(tcapi::get_elem(ctx,a,{0,0})!=E(99),"region copies independent");
    const auto sentinel=tcapi::get_elem(ctx,out,{0,0});
    throws<std::invalid_argument>([&]{tcapi::expand(ctx,a,{{0,-1}},out);});
    throws<std::out_of_range>([&]{tcapi::expand(ctx,a,{{-1,1}},out);});
    throws<std::overflow_error>([&]{tcapi::expand(ctx,a,{{0,std::numeric_limits<int32_t>::max()}},out);});
    throws<std::out_of_range>([&]{tcapi::shrink(ctx,a,{{2,{0,1}}},out);});
    throws<std::out_of_range>([&]{tcapi::shrink(ctx,a,{{0,{-1,1}}},out);});
    throws<std::out_of_range>([&]{tcapi::shrink(ctx,a,{{0,{1,0}}},out);});
    throws<std::out_of_range>([&]{tcapi::shrink(ctx,a,{{1,{0,4}}},out);});
    throws<std::invalid_argument>([&]{tcapi::shrink(ctx,a,{{1,{1,1}}},out);});
    throws<std::invalid_argument>([&]{tcapi::extract_sub(ctx,a,{{0,2}},out);});
    throws<std::invalid_argument>([&]{tcapi::replace_sub(ctx,a,small,{0},out);});
    throws<std::out_of_range>([&]{tcapi::replace_sub(ctx,a,small,{-1,0},out);});
    throws<std::out_of_range>([&]{tcapi::replace_sub(ctx,a,small,{2,0},out);});
    check(tcapi::get_elem(ctx,out,{0,0})==sentinel,"invalid region operations preserve output");
    auto scalar=tcapi::fill<T>(ctx,{},value<E>(4,-1));
    tcapi::expand(ctx,scalar,{},out);
    check(out.GetRaw()==nullptr,"scalar expansion replaces dense storage");
    tcapi::shrink(ctx,out,{}); tcapi::extract_sub(ctx,out,{});
    check(tcapi::get_elem(ctx,out,{})==value<E>(4,-1),"scalar region identity");
    tcapi::replace_sub(ctx,out,scalar,{});
    check(tcapi::get_elem(ctx,out,{})==value<E>(4,-1),"scalar replacement");
    throws<std::out_of_range>([&]{tcapi::expand(ctx,scalar,{{0,1}});});
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::expand(ctx,a,{});});
    throws<std::logic_error>([&]{tcapi::shrink(ctx,a,{});});
    throws<std::logic_error>([&]{tcapi::extract_sub(ctx,a,{});});
    throws<std::logic_error>([&]{tcapi::replace_sub(ctx,a,a,{0,0});});
}

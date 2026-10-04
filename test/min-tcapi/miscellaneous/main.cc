#include "../common/check.h"
#include "close_convert.h"
template<class E> void test() {
    test_close<E>();
    test_convert_from<E>();
    if constexpr(std::is_same_v<E,double>) test_conversion_specials();
    using T = gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    std::vector<E> input;
    for (int i=0;i<6;++i) input.push_back(value<E>(i+1, 2*i+1));
    auto map=[](const tcapi::elem_coors_t<T>& c) -> std::ptrdiff_t { return 3*c[0]+c[1]; };
    auto a=tcapi::assign_from_range<T>(ctx,{2,3},input.begin(),map);
    for(int i=0;i<2;++i) for(int j=0;j<3;++j)
        check(tcapi::get_elem(ctx,a,{i,j}) == input[3*i+j], "mapped input coordinates");
    std::vector<E> output(6);
    tcapi::to_range(ctx,a,output.begin(),map);
    check(input==output,"mapped range round trip");
    int calls=0;
    auto scalar_map=[&](const tcapi::elem_coors_t<T>& c) -> std::ptrdiff_t {
        ++calls; check(c.empty(),"scalar coordinates"); return 0;
    };
    auto scalar=tcapi::assign_from_range<T>(ctx,{},input.begin(),scalar_map);
    tcapi::to_range(ctx,scalar,output.begin(),scalar_map);
    check(calls==2 && output[0]==input[0],"scalar mapping each called once");
    a.Scale(E(2));
    tcapi::to_range(ctx,a,output.begin(),map);
    for(int i=0;i<6;++i) check(output[i]==E(2)*input[i],"scaled range export");
    throws<std::out_of_range>([&] { tcapi::to_range(ctx,a,output.begin(),[](const auto&) { return -1; }); });
    throws<std::out_of_range>([&] { tcapi::to_range(ctx,a,output.begin(),[](const auto&) { return 6; }); });
    throws<std::out_of_range>([&] { tcapi::assign_from_range<T>(ctx,{},input.begin(),[](const auto&) { return -1; }); });
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&] { tcapi::to_range(ctx,a,output.begin(),map); });
}
RUN_FOUR_TYPES(test)

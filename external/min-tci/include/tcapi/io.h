#pragma once
#include "detail/core.h"
#include <filesystem>
#include <fstream>
#include <string_view>

namespace tcapi {
namespace detail {
inline std::streamsize io_size(std::size_t bytes) {
    if (bytes>static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
        throw std::overflow_error("min-tcapi: I/O byte count exceeds streamsize");
    return static_cast<std::streamsize>(bytes);
}
inline void read_bytes(std::istream& in, void* data, std::size_t bytes) {
    in.read(static_cast<char*>(data),io_size(bytes));
    if (!in) throw std::ios_base::failure("min-tcapi: incomplete tensor stream");
}
template<class TenT> TenT read_tensor(std::istream& in) {
    // Preserve gqten's native binary format; labels are not part of that format.
    std::size_t rank=0;read_bytes(in,&rank,sizeof(rank));
    if (rank>static_cast<std::size_t>(std::numeric_limits<order_t<TenT>>::max()))
        throw std::overflow_error("min-tcapi: serialized rank exceeds order limit");
    shape_t<TenT> shape;
    // Read dimensions incrementally so a truncated rank header cannot trigger
    // a huge metadata allocation before the first dimension is read.
    for (std::size_t i=0;i<rank;++i) {
        bond_dim_t<TenT> dim;read_bytes(in,&dim,sizeof(dim));
        if (dim<=0) throw std::invalid_argument("min-tcapi: serialized dimensions must be positive");
        shape.push_back(dim);
    }
    const auto count=checked_size<TenT>(shape);
    using E=elem_t<TenT>;
    E scale;read_bytes(in,&scale,sizeof(E));
    if (shape.empty()) return construct<TenT>(shape,[&](std::size_t){return scale;});
    const auto bytes=count*sizeof(E);io_size(bytes);
    std::unique_ptr<E,decltype(&std::free)> data(static_cast<E*>(std::malloc(bytes)),&std::free);
    if (!data) throw std::bad_alloc();
    read_bytes(in,data.get(),bytes);
    TenT result(shape,data.get(),scale);data.release();
    return result;
}
template<class TenT> void write_tensor(const TenT& a,std::ostream& out) {
    const auto count=checked_size<TenT>(a.Shape());
    io_size(a.Rank()*sizeof(bond_dim_t<TenT>));
    if (a.Rank()) {
        io_size(count*sizeof(elem_t<TenT>));
        if (!a.DataMemInitialized()) throw std::invalid_argument("min-tcapi: tensor has no allocated data");
    }
    a.StreamWrite(out);
    if (!out) throw std::ios_base::failure("min-tcapi: tensor write failed");
}
} // namespace detail

template<class TenT,class Storage>
TenT load(context_handle_t<TenT>& ctx, Storage&& storage) {
    detail::verbose::call diagnostic("load", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "storage", storage);
    });
    detail::require_context(ctx);
    if constexpr(std::is_base_of_v<std::istream,std::decay_t<Storage>>) {
        return detail::read_tensor<TenT>(storage);
    } else {
        const std::filesystem::path path(std::forward<Storage>(storage));
        std::ifstream in(path,std::ios::binary);
        if (!in) throw std::ios_base::failure("min-tcapi: cannot open tensor file for reading");
        return detail::read_tensor<TenT>(in);
    }
}
template<class TenT,class Storage>
void save(context_handle_t<TenT>& ctx,const TenT& a, Storage&& storage) {
    detail::verbose::call diagnostic("save", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "storage", storage);
    });
    detail::require_context(ctx);
    if constexpr(std::is_base_of_v<std::ostream,std::decay_t<Storage>>) {
        detail::write_tensor(a,storage);
    } else {
        const std::filesystem::path path(std::forward<Storage>(storage));
        std::ofstream out(path,std::ios::binary|std::ios::trunc);
        if (!out) throw std::ios_base::failure("min-tcapi: cannot open tensor file for writing");
        detail::write_tensor(a,out);
        out.close();
        if (!out) throw std::ios_base::failure("min-tcapi: tensor file close failed");
    }
}
} // namespace tcapi

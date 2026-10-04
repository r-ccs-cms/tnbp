#pragma once
#include "../types.h"
#include "verbose.h"
#include "gqten/tensor/tensor.h"
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>

namespace tcapi::detail {
inline void require_context(const gqten_handle& ctx) {
    if (!ctx.active) throw std::logic_error("min-tcapi: context is not active");
}
template<class TenT>
std::size_t checked_size(const shape_t<TenT>& shape) {
    if (shape.size() > static_cast<std::size_t>(std::numeric_limits<order_t<TenT>>::max()))
        throw std::overflow_error("min-tcapi: tensor order overflow");
    std::size_t count = 1;
    for (const auto dim : shape) {
        if (dim <= 0) throw std::invalid_argument("min-tcapi: bond dimensions must be positive; use {} for a scalar");
        const auto d = static_cast<std::size_t>(dim);
        if (count > std::numeric_limits<std::size_t>::max() / d)
            throw std::overflow_error("min-tcapi: element count overflow");
        count *= d;
    }
    if (count > std::numeric_limits<std::size_t>::max() / sizeof(elem_t<TenT>))
        throw std::overflow_error("min-tcapi: storage size overflow");
    return count;
}
template<class TenT>
void check_coordinates(const TenT& a, const elem_coors_t<TenT>& coors) {
    if (coors.size() != a.Rank()) throw std::invalid_argument("min-tcapi: coordinate rank mismatch");
    for (std::size_t i = 0; i < coors.size(); ++i)
        if (coors[i] < 0 || coors[i] >= a.Shape()[i])
            throw std::out_of_range("min-tcapi: coordinate outside tensor");
}
template<class TenT>
elem_coors_t<TenT> coordinates(std::size_t index, const shape_t<TenT>& shape) {
    elem_coors_t<TenT> coors(shape.size());
    for (std::size_t i = 0; i < shape.size(); ++i) {
        coors[i] = static_cast<elem_coor_t<TenT>>(index % static_cast<std::size_t>(shape[i]));
        index /= static_cast<std::size_t>(shape[i]);
    }
    return coors;
}
// Keep ownership until the tensor has successfully copied its shape metadata.
template<class TenT, class Fill>
TenT construct(const shape_t<TenT>& shape, Fill&& fill) {
    const auto count = checked_size<TenT>(shape);
    using E = elem_t<TenT>;
    if (shape.empty()) {
        TenT a;
        a.SetElem(elem_coors_t<TenT>{}, static_cast<E>(fill(0)));
        return a;
    }
    std::unique_ptr<E, decltype(&std::free)> data(
        static_cast<E*>(std::malloc(count * sizeof(E))), &std::free);
    if (!data) throw std::bad_alloc();
    for (std::size_t i = 0; i < count; ++i) data.get()[i] = static_cast<E>(fill(i));
    TenT a(shape, data.get());
    data.release();
    return a;
}
// gqten move assignment from a scalar does not release the old dense payload.
// Destruction + default construction is needed for clear's release contract.
template<class TenT> void reset(TenT& a) {
    a.~TenT();
    ::new (static_cast<void*>(std::addressof(a))) TenT();
}
// replacement must be a distinct, fully constructed temporary. Reset first so
// gqten's scalar move assignment cannot retain the destination's old payload.
template<class TenT> void replace(TenT& destination, TenT&& replacement) {
    reset(destination);
    destination = std::move(replacement);
}
}

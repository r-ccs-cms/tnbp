#ifndef TCI_GQTEN_IO_H
#define TCI_GQTEN_IO_H

#include <string>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <system_error>
#include <filesystem>

namespace tci {

  template <typename TenT, typename Storage,
	    std::enable_if_t<is_gqten_tensor_v<TenT>,int> =0>
  TenT load(
	    tci::context_handle_t<TenT> &ctx,
	    Storage &&strg
	    ) {
    (void)ctx;
    using S = std::decay_t<Storage>;
    TenT A;
    if constexpr (std::is_base_of_v<std::istream,S>) {
      std::istream & in = strg;
      A.StreamRead(in);
    } else if constexpr (std::is_convertible_v<S,std::filesystem::path>) {
      std::filesystem::path p{std::forward<Storage>(strg)};
      std::ifstream ifs(p,std::ios::binary);
      if (!ifs) {
	throw std::system_error(errno, std::generic_category(),
				"Failed to open file for read: " + p.string());
      }
      A.StreamRead(ifs);
    } else {
      static_assert(!std::is_same_v<S,S>, "Unsupported Storage type for load()");
    }
    return A;
  }


  template <typename TenT, typename Storage,
	    std::enable_if_t<is_gqten_tensor_v<TenT>, int> = 0>
  void save(tci::context_handle_t<TenT>& ctx,
	    const TenT& a, Storage&& strg) {
    (void)ctx;
    using S = std::decay_t<Storage>;
    
    if constexpr (std::is_base_of_v<std::ostream, S>) {
      std::ostream& out = strg;
      a.StreamWrite(out);
    } else if constexpr (std::is_convertible_v<S, std::filesystem::path>) {
      std::filesystem::path p{std::forward<Storage>(strg)};
      std::ofstream ofs(p, std::ios::binary);
      if (!ofs) {
	throw std::system_error(errno, std::generic_category(),
				"Failed to open file for write: " + p.string());
      }
      a.StreamWrite(ofs);
    }
    else {
        static_assert(!std::is_same_v<S, S>,
                      "tci::save(): Unsupported Storage type");
    }
  }

}

#endif

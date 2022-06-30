// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_split.h>  // for SplitIterator, Splitter, StrSplit
#include <fmt/format.h>              // for format

#include <boost/process/pipe.hpp>
#include <string>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/numeric_utils.hpp"
#include "modle/common/random.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::test::cmatrix {

[[nodiscard]] inline std::vector<std::vector<u32>> load_matrix_from_file(
    const std::string& path_to_file, const std::string& sep = "\t") {
  std::vector<std::vector<u32>> m;
  compressed_io::Reader r(path_to_file);
  std::string line;
  std::string buff;
  u32 n{};
  while (r.getline(line)) {
    std::vector<u32> v;

    for (const auto& tok : absl::StrSplit(line, sep)) {
      modle::utils::parse_numeric_or_throw(tok, n);
      v.push_back(n);
    }
    m.emplace_back(std::move(v));
  }
  REQUIRE(r.eof());
  return m;
}

template <class ContactMatrixT>
static void write_cmatrix_to_stream(const ContactMatrixT& m, boost::process::opstream& s) {
  std::vector<contacts_t> buff(m.ncols());
  for (usize i = 0; i < m.ncols(); ++i) {
    buff.clear();
    for (usize j = 0; j < m.ncols(); ++j) {
      buff.push_back(m.get(i, j));
    }
    const auto sbuff = fmt::format(FMT_STRING("{}\n"), fmt::join(buff, ","));
    s.write(sbuff.data(), static_cast<std::streamsize>(sbuff.size()));
  }
  s.flush();
  s.pipe().close();
}

template <class ContactMatrixT>
[[nodiscard]] static ContactMatrixT read_cmatrix_from_stream(const usize ncols, const usize nrows,
                                                             boost::process::ipstream& s) {
  std::string sbuff;
  std::getline(s, sbuff);

  std::vector<double> buff;
  for (const auto& tok : absl::StrSplit(sbuff, ',')) {
    buff.push_back(utils::parse_numeric_or_throw<double>(tok));
  }
  REQUIRE(buff.size() == nrows * ncols);
  ContactMatrixT m(nrows, ncols);
  for (usize i = 0; i < nrows; ++i) {
    for (auto j = i; j < ncols; ++j) {
      m.set(i, j, buff[(i * m.nrows()) + j]);
    }
  }
  return m;
}

template <class ContactMatrix>
inline void create_random_matrix(ContactMatrix& m, usize nnz, u64 seed = 8336046165695760686ULL) {
  using N = typename ContactMatrix::value_type;
  assert(nnz <= m.npixels());

  auto rand_eng = random::PRNG(seed);

  auto contact_gen = [&rand_eng]() {
    if constexpr (std::is_floating_point_v<N>) {
      return random::uniform_real_distribution<N>{1, 65553}(rand_eng);
    } else {
      u64 max_ = std::min(u64(65553), static_cast<u64>((std::numeric_limits<N>::max)()));
      return random::uniform_int_distribution<N>{1, static_cast<N>(max_)}(rand_eng);
    }
  };

  auto row_gen = [&]() {
    return random::uniform_int_distribution<usize>{0, m.ncols() - 1}(rand_eng);
  };

  auto col_gen = [&]() {
    return random::uniform_int_distribution<usize>{0, m.nrows() - 1}(rand_eng);
  };

  do {
    for (usize i = m.unsafe_get_nnz(); i < nnz; ++i) {
      const auto row = row_gen();
      const auto col = std::min(row + col_gen(), m.ncols() - 1);
      if (row < m.ncols()) {
        m.set(row, col, contact_gen());
        assert(m.get_n_of_missed_updates() == 0);
      }
    }
  } while (m.unsafe_get_nnz() < nnz);
}
}  // namespace modle::test::cmatrix

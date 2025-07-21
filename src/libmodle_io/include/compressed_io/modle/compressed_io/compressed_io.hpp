// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <archive.h>

#include <array>
#include <boost/iostreams/filtering_stream.hpp>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

struct archive;
struct archive_entry;

namespace modle::compressed_io {
using namespace std::literals::string_view_literals;
DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
class Reader {
  DISABLE_WARNING_POP
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_read_free)>;

 public:
  Reader() = default;
  explicit Reader(const std::filesystem::path& path, usize buff_capacity = 512 * 1024);

  bool getline(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string_view getline(char sep = '\n');
  bool readall(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string readall(char sep = '\n');
  [[nodiscard]] bool eof() const noexcept;
  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const std::filesystem::path& path);
  void reset();

  [[nodiscard]] explicit operator bool() const;
  [[nodiscard]] bool operator!() const;

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

 private:
  std::filesystem::path _path{};
  archive_ptr_t _arc{nullptr, archive_read_free};
  std::unique_ptr<archive_entry*> _arc_entry{new archive_entry* };
  std::string _buff{};
  std::string _tok_tmp_buff{};
  usize _idx{0};
  bool _eof{false};

  void handle_libarchive_errors(la_ssize_t errcode) const;
  void handle_libarchive_errors() const;
  // Returns false when reaching eof
  [[nodiscard]] bool read_next_chunk();
  // Return false when unable to find the next token occurrence
  [[nodiscard]] bool read_next_token(std::string& buff, char sep);
  [[nodiscard]] std::string_view read_next_token(char sep);
};

DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
class Writer {
  DISABLE_WARNING_POP
 public:
  enum Compression : u8f { AUTO = 0, NONE = 1, GZIP = 2, BZIP2 = 3, LZMA = 4, ZSTD = 5 };

  Writer() = default;
  explicit Writer(const std::filesystem::path& path, Compression compression = AUTO);

  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const std::filesystem::path& path);

  void write(std::string_view buff);

  [[nodiscard]] explicit operator bool() const;
  [[nodiscard]] bool operator!() const;

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

  [[nodiscard]] static Compression infer_compression_from_ext(const std::filesystem::path& p);

 private:
  std::filesystem::path _path{};
  std::ofstream _fp{};
  boost::iostreams::filtering_ostream _out{};
  Compression _compression{AUTO};

  static constexpr std::array<std::pair<std::string_view, Compression>, 6> ext_mappings{
      std::make_pair(".gz"sv, GZIP),  std::make_pair(".bz2"sv, BZIP2),
      std::make_pair(".xz"sv, LZMA),  std::make_pair(".lzma"sv, LZMA),
      std::make_pair(".zst"sv, ZSTD), std::make_pair(".zstd"sv, ZSTD)};
};

}  // namespace modle::compressed_io

// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <archive.h>

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

#include "modle/common/common.hpp"
#include "modle/common/suppress_compiler_warnings.hpp"

namespace modle::compressed_io {

namespace detail {

struct ArchiveReaderDeleter {
  void operator()(archive* arch) const noexcept;
};

struct ArchiveWriterDeleter {
  void operator()(archive* arch) const noexcept;
};

struct ArchiveEntryDeleter {
  void operator()(archive_entry* entry) const noexcept;
};

}  // namespace detail

using namespace std::literals::string_view_literals;
DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
class Reader {
  DISABLE_WARNING_POP
 public:
  using archive_ptr_t = std::unique_ptr<archive, detail::ArchiveReaderDeleter>;
  using archive_entry_ptr_t = std::unique_ptr<archive_entry, detail::ArchiveEntryDeleter>;

  Reader() = default;
  explicit Reader(const std::filesystem::path& path, std::size_t buff_capacity = 512 * 1024);

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
  archive_ptr_t _arc{};
  archive_entry_ptr_t _arc_entry{};
  std::string _buff{};
  std::string _tok_tmp_buff{};
  std::size_t _idx{0};
  bool _eof{false};

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
  enum Compression : std::uint_fast8_t { AUTO, NONE, GZIP, BZIP2, LZMA, XZ, ZSTD };

  using archive_ptr_t = std::unique_ptr<archive, detail::ArchiveWriterDeleter>;
  using archive_entry_ptr_t = std::unique_ptr<archive_entry, detail::ArchiveEntryDeleter>;

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
  archive_ptr_t _arc{};
  archive_entry_ptr_t _arc_entry{};
  Compression _compression{AUTO};
};

}  // namespace modle::compressed_io

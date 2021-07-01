#pragma once
#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <archive.h>                       // for archive
#include <archive_entry.h>                 // for archive_entry_new, archive_entry_free

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/path.hpp>  // for path
#include <fstream>
#include <iostream>
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for string_view

namespace modle::compressed_io {
using namespace std::literals::string_view_literals;
class Reader {
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_read_free)>;

 public:
  Reader() = default;  // NOLINTNEXTLINE
  Reader(const boost::filesystem::path& path, size_t buff_capacity = 512 * 1024);

  bool getline(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string_view getline(char sep = '\n');
  [[nodiscard]] bool eof() const noexcept;
  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const boost::filesystem::path& path);
  void reset();

  [[nodiscard]] const boost::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

 private:
  boost::filesystem::path _path{};
  archive_ptr_t _arc{nullptr, archive_read_free};
  std::unique_ptr<archive_entry*> _arc_entry{new archive_entry* };
  std::string _buff{};
  std::string _tok_tmp_buff{};
  size_t _idx{0};
  bool _eof{false};

  void handle_libarchive_errors(la_ssize_t errcode) const;
  void handle_libarchive_errors() const;
  // Returns false when reaching eof
  [[nodiscard]] bool read_next_chunk();
  // Return false when unable to find the next token occurrence
  [[nodiscard]] bool read_next_token(std::string& buff, char sep);
  [[nodiscard]] std::string_view read_next_token(char sep);
};

class Writer {
 public:
  enum Compression : uint8_t { AUTO = 0, NONE = 1, GZIP = 2, BZIP2 = 3, LZMA = 4, ZSTD = 5 };

  Writer() = default;
  explicit Writer(const boost::filesystem::path& path, Compression compression = AUTO);  // NOLINT

  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const boost::filesystem::path& path);

  void write(std::string_view buff);

  [[nodiscard]] const boost::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

  [[nodiscard]] static Compression infer_compression_from_ext(const boost::filesystem::path& p);

 private:
  boost::filesystem::path _path{};
  Compression _compression{AUTO};
  std::ofstream _fp{};
  boost::iostreams::filtering_ostream _out{};

  static constexpr std::array<std::pair<std::string_view, Compression>, 6> ext_mappings{
      std::make_pair(".gz"sv, GZIP),  std::make_pair(".bz2"sv, BZIP2),
      std::make_pair(".xz"sv, LZMA),  std::make_pair(".lzma"sv, LZMA),
      std::make_pair(".zst"sv, ZSTD), std::make_pair(".zstd"sv, ZSTD)};
};

}  // namespace modle::compressed_io

#pragma once
#include <archive.h>        // for archive
#include <archive_entry.h>  // for archive_entry_new, archive_entry_free

#include <filesystem>   // for path
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for string_view

namespace modle::libarchivexx {
class Reader {
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_read_free)>;

 public:
  Reader() = default;  // NOLINTNEXTLINE
  Reader(const std::filesystem::path& path, size_t buff_capacity = 512 * 1024);

  bool getline(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string_view getline(char sep = '\n');
  [[nodiscard]] bool eof() const noexcept;
  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const std::filesystem::path& path);
  void reset();

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

 private:
  std::filesystem::path _path{};
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
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_write_free)>;

 public:
  enum Compression : uint8_t {
    NONE = ARCHIVE_FILTER_NONE,
    GZIP = ARCHIVE_FILTER_GZIP,
    BZIP2 = ARCHIVE_FILTER_BZIP2,
    LZMA = ARCHIVE_FILTER_LZMA,
    XZ = ARCHIVE_FILTER_XZ,
    LZIP = ARCHIVE_FILTER_LZIP,
    LRZIP = ARCHIVE_FILTER_LRZIP,
    LZOP = ARCHIVE_FILTER_LZOP,
    GRZIP = ARCHIVE_FILTER_GRZIP,
    LZ4 = ARCHIVE_FILTER_LZ4,
    ZSTD = ARCHIVE_FILTER_ZSTD
  };

  Writer() = default;
  Writer(const std::filesystem::path& path, Compression compression = GZIP);

  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const std::filesystem::path& path);

  void write(std::string_view buff);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;
  [[nodiscard]] const char* path_c_str() const noexcept;

 private:
  std::filesystem::path _path{};
  archive_ptr_t _arc{nullptr, archive_write_free};
  std::unique_ptr<archive_entry, decltype(&archive_entry_free)> _arc_entry{archive_entry_new(),
                                                                           archive_entry_free};
  Compression _compression{GZIP};
};

}  // namespace modle::libarchivexx

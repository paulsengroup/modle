#pragma once
#include <archive.h>  // for archive

#include <filesystem>   // for path
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for string_view

namespace modle::libarchivexx {
class Reader {
 public:
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_read_free)>;
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
}  // namespace modle::libarchivexx

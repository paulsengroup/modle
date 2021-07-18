#include "modle/compressed_io.hpp"

#include <archive.h>        // for archive
#include <archive_entry.h>  // for archive_entry
#include <fmt/format.h>     // system_error
#include <fmt/ostream.h>

#include <boost/filesystem/path.hpp>  // for path
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <cassert>      // for assert
#include <string>       // for string
#include <string_view>  // for string_view

#include "modle/common/utils.hpp"  // for ndebug_defined

namespace modle::compressed_io {

Reader::Reader(const boost::filesystem::path& path, size_t buff_capacity) {
  this->_buff.reserve(buff_capacity);
  this->open(path);
}

void Reader::open(const boost::filesystem::path& path) {
  auto handle_open_errors = [&](la_ssize_t status) {
    if (status == ARCHIVE_EOF) {
      throw std::runtime_error(fmt::format(FMT_STRING("File {} appears to be empty"), this->_path));
    }
    if (status < ARCHIVE_OK) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Failed to open file {} for reading (error code {}): {}"), this->_path,
          archive_errno(this->_arc.get()), archive_error_string(this->_arc.get())));
    }
  };

  if (this->is_open()) {
    this->close();
  }

  this->_path = path;
  if (this->_path.empty()) {
    return;
  }

  this->_arc.reset(archive_read_new());
  if (!this->_arc) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to allocate a buffer of to read file {}"), this->_path));
  }

  handle_open_errors(archive_read_support_filter_all(this->_arc.get()));
  handle_open_errors(archive_read_support_format_raw(this->_arc.get()));
  handle_open_errors(
      archive_read_open_filename(this->_arc.get(), this->_path.c_str(), this->_buff.capacity()));
  handle_open_errors(archive_read_next_header(this->_arc.get(), this->_arc_entry.get()));
  this->_idx = 0;
}

bool Reader::eof() const noexcept {
  assert(this->is_open());  // NOLINT
  return this->_eof;
}

bool Reader::is_open() const noexcept { return !!this->_arc; }

void Reader::close() {
  this->_arc = nullptr;
  this->_buff.clear();
  this->_eof = false;
}

void Reader::reset() {
  this->close();
  this->open(this->_path);
  this->_idx = 0;
  this->_buff.clear();
  this->_tok_tmp_buff.clear();
}

const boost::filesystem::path& Reader::path() const noexcept { return this->_path; }
std::string Reader::path_string() const noexcept { return this->_path.string(); }
const char* Reader::path_c_str() const noexcept { return this->_path.c_str(); }

void Reader::handle_libarchive_errors(la_ssize_t errcode) const {
  if (errcode < ARCHIVE_OK) {
    this->handle_libarchive_errors();
  }
}

void Reader::handle_libarchive_errors() const {
  if (const auto status = archive_errno(this->_arc.get()); status < ARCHIVE_OK) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("The following error occurred while reading file (error code {}): {}"),
        this->_path, archive_errno(this->_arc.get()), archive_error_string(this->_arc.get())));
  }
}

bool Reader::getline(std::string& buff, char sep) {
  assert(this->is_open());  // NOLINT
  buff.clear();
  if (this->eof()) {
    return false;
  }

  while (!this->read_next_token(buff, sep)) {
    if (!this->read_next_chunk()) {
      assert(this->eof());  // NOLINT
      return false;
    }
  }
  return true;
}

std::string_view Reader::getline(char sep) {
  assert(this->is_open());  // NOLINT
  if (this->eof()) {
    return std::string_view{};
  }

  this->_tok_tmp_buff.clear();
  while (true) {
    if (const auto tok = this->read_next_token(sep); !tok.empty()) {
      return tok;
    }
    if (!this->read_next_chunk()) {
      assert(this->eof());  // NOLINT
      return std::string_view{};
    }
  }
}

bool Reader::read_next_chunk() {
  assert(!this->eof());     // NOLINT
  assert(this->is_open());  // NOLINT
  this->_buff.resize(this->_buff.capacity());
  const auto bytes_read =
      archive_read_data(this->_arc.get(), this->_buff.data(), this->_buff.capacity());
  if (bytes_read < 0) {
    handle_libarchive_errors();
  } else if (bytes_read == 0) {
    this->_eof = true;
    this->_buff.clear();
    this->_tok_tmp_buff.clear();
    return false;
  }
  this->_buff.resize(static_cast<size_t>(bytes_read));
  this->_idx = 0;
  return true;
}

bool Reader::read_next_token(std::string& buff, char sep) {
  assert(!this->eof());                      // NOLINT
  assert(this->is_open());                   // NOLINT
  assert(this->_idx <= this->_buff.size());  // NOLINT
  if (this->_idx == this->_buff.size()) {
    return false;
  }

  const auto pos = this->_buff.find(sep, this->_idx);
  const auto i = static_cast<int64_t>(this->_idx);
  if (pos == std::string::npos) {
    buff.append(this->_buff.begin() + i, this->_buff.end());
    return false;
  }

  assert(pos >= this->_idx);  // NOLINT
  buff.append(this->_buff.begin() + i, this->_buff.begin() + static_cast<int64_t>(pos));
  this->_idx = pos + 1;
  return true;
}

std::string_view Reader::read_next_token(char sep) {
  assert(!this->eof());                      // NOLINT
  assert(this->is_open());                   // NOLINT
  assert(this->_idx <= this->_buff.size());  // NOLINT
  if (this->_idx == this->_buff.size()) {
    return std::string_view{};
  }

  const auto pos = this->_buff.find(sep, this->_idx);
  const auto i = static_cast<int64_t>(this->_idx);
  if (pos == std::string::npos) {
    this->_tok_tmp_buff.append(this->_buff.begin() + i, this->_buff.end());
    return std::string_view{};
  }

  assert(pos >= this->_idx);  // NOLINT
  this->_idx = pos + 1;
  if (this->_tok_tmp_buff.empty()) {  // NOLINTNEXTLINE
    return std::string_view{this->_buff.data() + static_cast<size_t>(i),
                            pos - static_cast<size_t>(i)};
  }

  this->_tok_tmp_buff.append(this->_buff.begin() + i,
                             this->_buff.begin() + static_cast<int64_t>(pos));
  DISABLE_WARNING_PUSH
#if __GNUC__ < 8
  // The following two warnings seem to be a false positive produced by GCC 7.5
  DISABLE_WARNING_CONVERSION
#endif
  return std::string_view{this->_tok_tmp_buff};
}
DISABLE_WARNING_POP

Writer::Writer(const boost::filesystem::path& path, Compression compression)
    : _compression(compression) {
  this->open(path);
}

void Writer::open(const boost::filesystem::path& path) {
  this->_out.reset();
  if (this->is_open()) {
    this->close();
  }

  if (this->_compression == AUTO) {
    this->_compression = infer_compression_from_ext(path);
  }

  switch (this->_compression) {
    case GZIP:
      this->_out.push(boost::iostreams::gzip_compressor(
          boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
      break;
    case BZIP2:
      this->_out.push(boost::iostreams::bzip2_compressor(boost::iostreams::bzip2_params(9)));
      break;
    case LZMA:
      this->_out.push(boost::iostreams::lzma_compressor(
          boost::iostreams::lzma_params(boost::iostreams::lzma::best_compression)));
      break;
    case ZSTD:  // Disabled
      // this->_out.push(boost::iostreams::zstd_compressor(
      //     boost::iostreams::zstd_params(boost::iostreams::zstd::best_compression)));
      throw std::runtime_error("ZSTD compression not available.");
    case NONE:
      break;
    case AUTO:
      if constexpr (utils::ndebug_defined()) {
        utils::throw_with_trace(std::logic_error("Unreachable code"));
      }
  }
  this->_path = path;
  this->_fp.open(path.string(), std::ios_base::binary);
  if (!this->_fp) {
    throw fmt::system_error(errno, FMT_STRING("Failed to open file {} for writing"), this->_path);
  }
  this->_out.push(this->_fp);
}

bool Writer::is_open() const noexcept { return this->_fp.is_open(); }

void Writer::close() {
  this->_out.reset();
  boost::iostreams::close(this->_out);
  this->_fp.close();
}

const boost::filesystem::path& Writer::path() const noexcept { return this->_path; }
std::string Writer::path_string() const noexcept { return this->_path.string(); }
const char* Writer::path_c_str() const noexcept { return this->_path.c_str(); }

Writer::Compression Writer::infer_compression_from_ext(const boost::filesystem::path& p) {
  const auto ext = absl::AsciiStrToLower(p.extension().string());
  const auto* const match = std::find_if(ext_mappings.begin(), ext_mappings.end(),
                                         [&](const auto& mapping) { return mapping.first == ext; });
  if (match == ext_mappings.end()) {
    return NONE;
  }
  return match->second;
}

void Writer::write(std::string_view buff) {
  if constexpr (utils::ndebug_defined()) {
    if (!this->is_open()) {
      utils::throw_with_trace(std::runtime_error("Writer::write() was called on a closed file!"));
    }
  }
  this->_out << buff;
}

}  // namespace modle::compressed_io

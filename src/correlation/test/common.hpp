#pragma once
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_set>

namespace modle::correlation::test {

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
void write_vect_to_file(const std::string& fpath, const std::vector<N>& v) {
  auto fp = std::ofstream(fpath);
  fp << v[0];
  for (auto i = 1UL; i < v.size(); ++i) {
    fp << "," << v[i];
  }
  fp << std::endl;
  fp.close();
}

std::vector<uint32_t> generate_random_vect(std::mt19937& rnd_eng, uint32_t size, uint32_t min,
                                           uint32_t max, bool allow_duplicates = true) {
  std::uniform_int_distribution<uint32_t> dist(min, max);
  std::vector<uint32_t> v(size);
  if (allow_duplicates) {
    std::generate(v.begin(), v.end(), [&]() { return dist(rnd_eng); });
  } else {
    std::unordered_set<uint32_t> s;
    while (s.size() < size) {
      s.insert(dist(rnd_eng));
    }
    v = {s.begin(), s.end()};
  }
  return v;
}

std::pair<std::vector<uint32_t>, std::vector<uint32_t>> generate_correlated_vects(
    std::mt19937& rnd_eng, uint32_t size) {
  std::uniform_int_distribution<int32_t> dist(static_cast<int32_t>(size) / -50,  // NOLINT
                                              static_cast<int32_t>(size / 50));  // NOLINT
  std::vector<uint32_t> v1(size);
  std::vector<uint32_t> v2(size);
  std::iota(v1.begin(), v1.end(), 0);
  std::iota(v2.begin(), v2.end(), 0);
  for (auto i = 0UL; i < size; ++i) {
    int64_t n = static_cast<int64_t>(v1[i]) + dist(rnd_eng);
    v1[i] = static_cast<uint32_t>(std::max(0L, n));
    n = static_cast<int64_t>(v2[i]) + dist(rnd_eng);
    v2[i] = static_cast<uint32_t>(std::max(0L, n));
  }
  return {v1, v2};
}

template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
std::pair<double, double> corr_scipy(const std::vector<N>& v1, const std::vector<N>& v2,
                                     const std::string& method) {
  auto f1_path = std::string(std::tmpnam(nullptr));  // NOLINT
  auto f2_path = std::string(std::tmpnam(nullptr));  // NOLINT
  write_vect_to_file(f1_path, v1);
  write_vect_to_file(f2_path, v2);

  const std::string cmd =
      "python3 -c '"
      "from scipy.stats import pearsonr, spearmanr, kendalltau; from sys import argv, stderr; "
      "from numpy import genfromtxt; "
      "v1 = genfromtxt(argv[1], delimiter=\",\", dtype=int); "
      "v2 = genfromtxt(argv[2], delimiter=\",\", dtype=int); "
      "corr, pv = " +
      method +
      "(v1, v2); "
      //      "print(v1, file=stderr); print(v2, file=stderr); "
      "print(f\"{corr:.16e}\\t{pv:.16e}\", end=\"\");' " +
      f1_path + " " + f2_path;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  std::array<char, 256> buffer{};
  std::string result;
  // TODO: replace this with boost::process
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
    //    printf("%s\n", result.c_str());
  }
  std::filesystem::remove(f1_path);
  std::filesystem::remove(f2_path);

  const auto rho = std::stod(std::string(result.data(), result.find('\t')));
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto pv = std::stod(std::string(result.data() + result.find('\t')));

  return {rho, pv};
}
}  // namespace modle::correlation::test
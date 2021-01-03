#pragma once

#include <catch2/catch.hpp>

#include "modle/juicer_contacts.hpp"

namespace modle::test::juicer_contacts {
using namespace modle::juicer_contacts;
/*
TEST_CASE("Juicer Contacts run_juicer_dump_and_parse_contacts", "[parsers][juicer_contacts][short]")
{ const std::string hic_file =
      "https://encode-public.s3.amazonaws.com/2018/10/19/37bcf091-9330-46bb-9a88-379b5920a41b/"
      "ENCFF956DEY.hic";
  ContactMatrix<uint32_t>::Header header{"1", 1000, 1'700'000, 1'800'000, 10'000, 100, 10};
  const auto cmatrix = run_juicer_dump_and_parse_contacts<uint32_t>(header, hic_file);
  cmatrix.write_to_tsv("/tmp/test_cmatrix.tsv.bz2", "#chr1\t1000\t1700000\t1800000\t10000\n");
}
 */
}  // namespace modle::test::juicer_contacts
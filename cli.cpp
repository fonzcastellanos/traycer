#include "cli.hpp"

#include <cstdio>
#include <cstring>

namespace cli {

const char* const kStatusMessages[] = {
    "success", "invalid option type", "missing option value",
    "invalid option value", "unexpected option name"};

const char* StatusMessage(Status s) { return kStatusMessages[s]; }

Status ParseOpts(uint argc, char* argv[], Opt opts[], uint opts_size,
                 uint* argi) {
  uint i = 1;

  while (i < argc && std::strncmp(argv[i], "--", 2) == 0) {
    if (std::strcmp(argv[i], "--") == 0) {
      ++i;
      break;
    }

    if (i + 1 >= argc) {
      *argi = i;
      return kStatus_MissingOptVal;
    }

    uint j = 0;

    while (j < opts_size && std::strcmp(&argv[i][2], opts[j].name) != 0) {
      ++j;
    }
    if (j == opts_size) {
      *argi = i;
      return kStatus_UnexpectedOptName;
    }

    const char* format;
    switch (opts[j].type) {
      case kOptType_Int: {
        format = "%d";
        break;
      }
      case kOptType_Uint: {
        format = "%u";
        break;
      }
      case kOptType_Float: {
        format = "%f";
        break;
      }
      case kOptType_String: {
        format = "%s";
        break;
      }
      default: {
        *argi = i;
        return kStatus_InvalidOptType;
      }
    }

    int rc = std::sscanf(argv[i + 1], format, opts[j].val);
    if (rc != 1) {
      *argi = i;
      return kStatus_InvalidOptVal;
    }

    i += 2;
  }

  *argi = i;

  return kStatus_Ok;
}
}  // namespace cli
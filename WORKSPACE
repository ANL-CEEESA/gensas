workspace(name = "gensas")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Rule repository, note that it's recommended to use a pinned commit to a released version of the rules
http_archive(
    name = "rules_foreign_cc",
    sha256 = "2a4d07cd64b0719b39a7c12218a3e507672b82a97b98c6a89d38565894cf7c51",
    strip_prefix = "rules_foreign_cc-0.9.0",
    url = "https://github.com/bazelbuild/rules_foreign_cc/archive/refs/tags/0.9.0.tar.gz",
)

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")

# This sets up some common toolchains for building targets. For more details, please see
# https://github.com/bazelbuild/rules_foreign_cc/tree/main/docs#rules_foreign_cc_dependencies
rules_foreign_cc_dependencies()

_ALL_CONTENT = """\
filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
"""

# pcre source code repository
http_archive(
    name = "pcre",
    build_file_content = _ALL_CONTENT,
    strip_prefix = "pcre-8.43",
    urls = [
        "https://mirror.bazel.build/ftp.pcre.org/pub/pcre/pcre-8.43.tar.gz",
        "https://ftp.pcre.org/pub/pcre/pcre-8.43.tar.gz",
    ],
    sha256 = "0b8e7465dc5e98c757cc3650a20a7843ee4c3edf50aaf60bb33fd879690d2c73",
)

_OPENBLAS_README = """
exports_files(["README.md"])
"""
# openblas source code repository
http_archive(
    name = "openblas",
    build_file_content = _ALL_CONTENT +_OPENBLAS_README ,
    strip_prefix = "OpenBLAS-0.3.23",
    urls = [
        "https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz",
    ],
    sha256 = "5D9491D07168A5D00116CDC068A40022C3455BF9293C7CB86A65B1054D7E5114",
)

# superlu
http_archive(
    name = "superlu",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "superlu-6.0.1",
    urls = [
        "https://github.com/xiaoyeli/superlu/archive/refs/tags/v6.0.1.tar.gz",
    ],
    sha256 = "6C5A3A9A224CB2658E9DA15A6034EED44E45F6963F5A771A6B4562F7AFB8F549",
)

#armadillo
http_archive(
    name = "armadillo",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "armadillo-code-12.6.x",
    urls = [
        "https://gitlab.com/conradsnicta/armadillo-code/-/archive/12.6.x/armadillo-code-12.6.x.tar.gz"
    ],
    sha256 = "0052C111E1188737A8B9495FD163D80015D3D8F0AE3D962B16BF74BD29176A23",
)

#hdf5
http_archive(
    name = "hdf5",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "./hdf5-1.14.3",
    urls = [
        "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz"
    ],
    sha256 = "09CDB287AA7A89148C1638DD20891FDBAE08102CF433EF128FD345338AA237C7",
)

#zlib
http_archive(
    name = "zlib",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "zlib-1.3",
    urls = [
        "https://zlib.net/zlib-1.3.tar.gz"
    ],
    sha256 = "FF0BA4C292013DBC27530B3A81E1F9A813CD39DE01CA5E0F8BF355702EFA593E",
)

#libaec (szlib)
http_archive(
    name = "libaec",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "libaec-v1.0.6",
    urls = [
        "https://gitlab.dkrz.de/k202009/libaec/-/archive/v1.0.6/libaec-v1.0.6.tar.gz"
    ],
    sha256 = "2675D26B24A081CDAAA48E0D085632E568E3F85F843F62B91FC51E1B3F5FE93F",
)

#libmatio
http_archive(
    name = "libmatio",
    build_file_content = _ALL_CONTENT ,
    strip_prefix = "matio-1.5.26",
    urls = [
        "https://github.com/tbeu/matio/releases/download/v1.5.26/matio-1.5.26.tar.gz"
    ],
    sha256 = "8B47C29F58E468DBA7A5555371C6A72AD4C6AA8B15F459B2B0B65A303C063933",
)

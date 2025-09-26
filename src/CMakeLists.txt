aux_source_directory(. SRC)
include_directories(.)

add_executable(main ${SRC}
        kzg.cpp
        svc.cpp
        basic.cpp
        fc.cpp
        pst.cpp
        pst.h)

find_package(OpenSSL REQUIRED)

target_link_libraries(main mcl XKCP  OpenSSL::SSL OpenSSL::Crypto)

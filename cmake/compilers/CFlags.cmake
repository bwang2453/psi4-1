if(NOT DEFINED DEFAULT_C_FLAGS_SET OR RESET_FLAGS)

# custom flags are defined by setup --custom-cc-flags
if(DEFINED CUSTOM_C_FLAGS)
    # set custom compiler flags (for every build type)
    set(CMAKE_C_FLAGS "${CUSTOM_C_FLAGS}")
    # special flags for build types will be empty
    set(CMAKE_C_FLAGS_DEBUG "")
    set(CMAKE_C_FLAGS_RELEASE "")
    set(CMAKE_C_FLAGS_PROFILE "")
else()
    # custom flags are not defined
    if(CMAKE_C_COMPILER_ID MATCHES GNU)
        set(CMAKE_C_FLAGS         "-std=c99 -DRESTRICT=${restrict} -DFUNDERSCORE=1 -fPIC")
        # Special debug flags are set if ASan, TSan or UBSan are requested
        if(ENABLE_ASAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=address -fno-omit-frame-pointer")
        elseif(ENABLE_TSAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=thread -fno-omit-frame-pointer -pie")
        elseif(ENABLE_UBSAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=undefined -fno-omit-frame-pointer")
        else()
         set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3 -Wall -Wextra -Winit-self -Wuninitialized -Wmissing-declarations -Wwrite-strings ")
        endif()
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_LINK_FLAGS    " ")
        set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")

    if(ENABLE_VECTORIZATION)
           set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_ARCHITECTURE_FLAGS} ${DEFINITIONS}")
        endif()
        if(ENABLE_CODE_COVERAGE)
           set(CMAKE_C_FLAGS      "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
           set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -fprofile-arcs -ftest-coverage")
        endif()

        if(ENABLE_STATIC_LINKING)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES Intel)
        set(CMAKE_C_FLAGS "-restrict -DRESTRICT=${restrict} -std=c99 -fPIC")
        set(CMAKE_C_FLAGS_DEBUG   "-O0 -g -w3 -vec-report -Wall -Wuninitialized ")
        # Check if xHost flag is available and add it CMAKE_C_FLAGS_RELEASE
        set(xHost "")
        if(ENABLE_XHOST)
            check_c_compiler_flag("-xHost" has_xHost)
            if(has_xHost)
                set(xHost "-xHost")
            endif()
        endif()
        set(CMAKE_C_FLAGS_RELEASE "-O3 -ip -DNDEBUG ${xHost}")
        set(CMAKE_C_LINK_FLAGS    "-shared-intel -fPIC")
        set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -g -pg")

    if(ENABLE_VECTORIZATION)
           set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_ARCHITECTURE_FLAGS} ${DEFINITIONS}")
        endif()

        if(DEFINED MKL_FLAG)
           set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_FLAG}")
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES Clang)
        set(CMAKE_C_FLAGS         "-std=c99 -DRESTRICT=${restrict} -DFUNDERSCORE=1 -fPIC")
        # Special debug flags are set if ASan, MSan, TSan or UBSan are requested
        if(ENABLE_ASAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=address -fno-omit-frame-pointer")
        elseif(ENABLE_MSAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=memory -fno-omit-frame-pointer")
        elseif(ENABLE_TSAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=thread -fno-omit-frame-pointer")
        elseif(ENABLE_UBSAN)
           set(CMAKE_C_FLAGS_DEBUG    "-g -O1 -fsanitize=undefined -fno-omit-frame-pointer")
        else()
           set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3 -Wall -Wextra -Winit-self -Wuninitialized -Wmissing-declarations -Wwrite-strings ")
        endif()
        set(CMAKE_C_FLAGS_RELEASE "-O3")
        if(ENABLE_VECTORIZATION)
       set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_ARCHITECTURE_FLAGS} ${DEFINITIONS}")
        endif()
        # clang does not use gprof
        set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE}")

        if(ENABLE_STATIC_LINKING)
           set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Bstatic")
        endif()
    else()
        message(FATAL_ERROR "Vendor of your C compiler is not supported")
    endif()

    if(DEFINED EXTRA_C_FLAGS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
    endif()
endif()

save_compiler_flags(C)
endif()

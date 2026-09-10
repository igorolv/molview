# Сборка «Строения молекул», версия на VTK.
#
# Почему CMake, а не Makefile как во второй версии: VTK распространяется
# вместе с файлами настройки CMake, и подключается она одной строкой
# find_package(VTK). Собирать её вручную, выписывая полсотни -l, бессмысленно.
#
#   cmake -B build -G "Ninja"      настроить
#   cmake --build build            собрать
#   ctest --test-dir build         прогнать проверки ядра
#
# Запускать из оболочки MINGW64, которую ставит MSYS2: только там в PATH
# лежит тот g++, которым собрана VTK из пакета pacman. Смешивать её
# с w64devkit нельзя — разные сборки libstdc++ не совместимы по ABI.

cmake_minimum_required(VERSION 3.20)
project(molviewvtk LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Без явного указания CMake собирает без оптимизации и без отладочной
# информации разом — то есть в самом бесполезном режиме.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Исходники в UTF-8, и строки в них должны остаться UTF-8 после компиляции:
# иначе русские сообщения превратятся в мусор ещё до Windows API.
add_compile_options(-Wall -Wextra -finput-charset=UTF-8 -fexec-charset=UTF-8)

# ── ядро расчёта ────────────────────────────────────────────────────────
# Про Windows не знает ничего. Поэтому один и тот же объектный код идёт
# и в оконную программу, и в консольные тесты.
file(GLOB CHEM_SRC CONFIGURE_DEPENDS src/chem/*.cpp)
add_library(chem STATIC ${CHEM_SRC})
target_include_directories(chem PUBLIC src)

# ── VTK ─────────────────────────────────────────────────────────────────
# Просить только нужные модули, а не всё сразу: полный список тянет на
# компоновку сотню библиотек и превращает сборку в минуту ожидания.
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    DomainsChemistry
    RenderingCore
    RenderingOpenGL2
    RenderingFreeType
)

# ── оконная программа ───────────────────────────────────────────────────
file(GLOB GUI_SRC CONFIGURE_DEPENDS src/render/*.cpp src/scene/*.cpp src/ui/*.cpp)
add_executable(molview WIN32 ${GUI_SRC} src/main.cpp)
# msimg32 — ради AlphaBlend, которым кадр VTK ложится на фон от GDI+.
target_link_libraries(molview PRIVATE chem ${VTK_LIBRARIES}
                      gdiplus comctl32 gdi32 msimg32 user32 ole32)

# WIN32 выше даёт -mwindows (программа без консоли за спиной), а точка входа
# у неё wWinMain с широкими строками — за это отвечает -municode.
target_link_options(molview PRIVATE -municode)

# Единого статического .exe больше нет: VTK собрана динамически, и
# статически с ней не слинковаться. Рядом с программой лягут DLL.

# БЕЗ ЭТОЙ СТРОКИ БУДЕТ ЧЁРНЫЙ ЭКРАН. Она регистрирует фабрики реализаций
# (какой рендер-бэкенд подставить вместо абстрактного vtkRenderWindow).
# Ошибок при этом нет ни при сборке, ни при запуске — просто ничего не рисуется.
vtk_module_autoinit(TARGETS molview MODULES ${VTK_LIBRARIES})

# ── проверки ядра ───────────────────────────────────────────────────────
# Консольные, про Windows не знают, поэтому WIN32 у них нет.
enable_testing()
foreach(name selftest smilestest)
    add_executable(${name} tests/${name}.cpp)
    target_link_libraries(${name} PRIVATE chem)
    target_include_directories(${name} PRIVATE tests)
    target_link_options(${name} PRIVATE -static -static-libgcc -static-libstdc++)
    # Рабочий каталог — корень проекта: findDataFile ищет data/ относительно него.
    add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()

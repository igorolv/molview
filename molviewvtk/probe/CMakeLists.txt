# Проба собирается отдельно от основного проекта: она нужна один раз,
# до начала работы, и тащить её в общую сборку незачем.
#
#   cmake -B build -G "Ninja"
#   cmake --build build
#   ./build/probe.exe

cmake_minimum_required(VERSION 3.20)
project(probe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Просить у VTK только нужные модули, а не всё сразу: полный список тянет
# на компоновку сотню библиотек и превращает сборку в минуту ожидания.
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    DomainsChemistry
    RenderingCore
    RenderingOpenGL2
    RenderingFreeType
)

add_executable(probe WIN32 probe.cpp)
target_compile_options(probe PRIVATE -Wall -Wextra -finput-charset=UTF-8 -fexec-charset=UTF-8)
target_link_options(probe PRIVATE -municode)
target_link_libraries(probe PRIVATE ${VTK_LIBRARIES} gdi32 user32)

# БЕЗ ЭТОЙ СТРОКИ БУДЕТ ЧЁРНЫЙ ЭКРАН. Она разворачивается в определения,
# которые регистрируют фабрики реализаций (какой именно рендер-бэкенд
# подставить вместо абстрактного vtkRenderWindow). Ошибок при этом нет
# ни при сборке, ни при запуске — просто ничего не рисуется.
vtk_module_autoinit(TARGETS probe MODULES ${VTK_LIBRARIES})

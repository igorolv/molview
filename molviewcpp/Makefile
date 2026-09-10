# Сборка проекта «Строение молекул».
# Тулчейн — w64devkit (g++ + make), никаких внешних библиотек.
#
#   make            собрать программу и тесты
#   make test       прогнать обе проверки ядра
#   make probe      проба из раздела 3 задания: собирается ли GDI+
#   make clean      убрать промежуточные файлы

CXX      := g++
# -MMD -MP просят компилятор выписать зависимости от заголовков: без них
# правка .h не приводит к пересборке и в компоновку попадают устаревшие объекты
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Isrc -Itests -MMD -MP \
            -finput-charset=UTF-8 -fexec-charset=UTF-8
LDFLAGS  := -static -static-libgcc -static-libstdc++
# GDI+ идёт ПОСЛЕ объектных файлов, иначе компоновщик её не подхватит
GUILIBS  := -lgdiplus -lcomctl32 -lgdi32 -luser32 -lole32

BUILD    := build

# Вычислительное ядро. Про Windows не знает ничего — именно поэтому
# тесты собираются как обычные консольные программы.
CHEM_SRC := $(wildcard src/chem/*.cpp)
CHEM_OBJ := $(patsubst src/%.cpp,$(BUILD)/%.o,$(CHEM_SRC))

# Рендер и интерфейс — только для оконной программы.
GUI_SRC  := $(wildcard src/render/*.cpp) $(wildcard src/ui/*.cpp) src/main.cpp
GUI_OBJ  := $(patsubst src/%.cpp,$(BUILD)/%.o,$(GUI_SRC))

APP      := $(BUILD)/molview.exe
TESTS    := $(BUILD)/selftest.exe $(BUILD)/smilestest.exe

.PHONY: all test probe clean

all: $(APP) $(TESTS)

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD)/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# -mwindows убирает консольное окно за спиной программы,
# -municode выбирает точкой входа wWinMain с широкими строками
$(APP): $(GUI_OBJ) $(CHEM_OBJ)
	$(CXX) $^ -o $@ -mwindows -municode $(LDFLAGS) $(GUILIBS)

$(BUILD)/selftest.exe: $(BUILD)/selftest.o $(CHEM_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(BUILD)/smilestest.exe: $(BUILD)/smilestest.o $(CHEM_OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

test: $(TESTS)
	@$(BUILD)/selftest.exe && $(BUILD)/smilestest.exe

# Проба GDI+. Отдельная цель: ядро от Windows не зависит и собирается без неё.
probe: $(BUILD)/probe.exe
	@$(BUILD)/probe.exe

$(BUILD)/probe.exe: probe.cpp
	@mkdir -p $(BUILD)
	$(CXX) -std=c++17 $< -o $@ -lgdiplus $(LDFLAGS)

clean:
	rm -rf $(BUILD)

# зависимости от заголовков, выписанные компилятором рядом с объектными файлами
-include $(GUI_OBJ:.o=.d) $(CHEM_OBJ:.o=.d) $(BUILD)/selftest.d $(BUILD)/smilestest.d

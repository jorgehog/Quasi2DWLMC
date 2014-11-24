CONFIG -= qt

SOURCES += main.cpp

LIBS += -L$$PWD/../WLMC/lib -lWLMC -larmadillo

QMAKE_CXXFLAGS += -std=c++11

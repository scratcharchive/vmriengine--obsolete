TEMPLATE =	app

OBJECTS_DIR = build
DESTDIR = bin
TARGET = vmriengine

INCLUDEPATH += .
DEPENDPATH += .
SOURCES += vmriengine.cpp
HEADERS += vmrisimulator.h vmriinterface.h
SOURCES += vmrisimulator.cpp vmriinterface.cpp
HEADERS += phantom1.h
SOURCES += phantom1.cpp
HEADERS += abstractphantom.h
SOURCES += abstractphantom.cpp

INCLUDEPATH += utils
DEPENDPATH += utils
HEADERS += textfile.h qjson.h
SOURCES += textfile.cpp qjson.cpp


LIBS += -lqjson
